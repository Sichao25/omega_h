#include "Omega_h_egads_lite.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_timer.hpp"
#include "Omega_h_for.hpp"
#include <Omega_h_file.hpp> //vtk

#include <Omega_h_fail.hpp>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#endif

#include <egads.h>

enum EgadsObjectClass {
  EGADS_CONTXT = CONTXT,
  EGADS_TRANSFORM = TRANSFORM,
  EGADS_TESSELATION = TESSELLATION,
  EGADS_NIL = NIL,
  /*EGADS_EMPTY = EMPTY, not doing this one
   * because an EGADS error exists by the same name
   */
  EGADS_REFERENCE = REFERENCE,
  EGADS_PCURVE = PCURVE,
  EGADS_CURVE = CURVE,
  EGADS_SURFACE = SURFACE,
  EGADS_NODE = NODE,
  EGADS_EDGE = EDGE,
  EGADS_LOOP = LOOP,
  EGADS_FACE = FACE,
  EGADS_SHELL = SHELL,
  EGADS_BODY = BODY,
  EGADS_MODEL = MODEL
};

#undef CONTXT
#undef TRANSFORM
#undef TESSELLATION
#undef NIL
#undef EMPTY
#undef REFERENCE
#undef PCURVE
#undef CURVE
#undef SURFACE
#undef NODE
#undef EDGE
#undef LOOP
#undef FACE
#undef SHELL
#undef BODY
#undef MODEL

#ifdef __clang__
#pragma clang diagnostic pop
#endif

namespace Omega_h {

OMEGA_H_INLINE void call_egads(
    int result, char const* code, char const* file, int line) {
  if (EGADS_SUCCESS == result) return;
  OMEGA_H_CHECK_PRINTF(false,
      "EGADS call %s returned %d at %s +%d\n", code, result, file, line);
}

#define CALL(f) call_egads((f), #f, __FILE__, __LINE__)

static int const dims2oclass[4] = {
    EGADS_NODE, EGADS_EDGE, EGADS_FACE, EGADS_BODY};

struct Egads {
  ego context;
  ego model;
  ego body;
  int counts[3];
  ego* entities[3];
  std::map<std::set<ego>, ego> classifier;
  LOs counts_d;
  ego* entities_d[3]; //HACK stored as arrays of GOs... see conversion fns below
};

OMEGA_H_INLINE Omega_h::GO egoToGo(ego obj) {
  return (Omega_h::GO) obj;
}

OMEGA_H_INLINE Omega_h::GO egoPtrToGo(ego* obj) {
  return (Omega_h::GO) obj;
}

OMEGA_H_INLINE ego* goToEgoPtr(Omega_h::GO obj) {
  return (ego*) obj;
}

OMEGA_H_INLINE ego goToEgo(Omega_h::GO obj) {
  return (ego) obj;
}

Omega_h::Write<Omega_h::GO> OhWriteEgo(int n) {
  assert(sizeof(Omega_h::GO) == sizeof(ego));
  return Omega_h::Write<Omega_h::GO>(n);
}

Egads* egads_lite_load(std::string const& filename) {
  auto eg = new Egads;
  CALL(EG_open(&eg->context));
  CALL(EG_loadModel(eg->context, 0, filename.c_str(), &eg->model));
  int nbodies;
  ego* bodies;

  for (int i = 0; i < 3; ++i)
    printf("dims2oclass[%d] %d\n", i, dims2oclass[i]);
  const auto egModel = eg->model;
  Omega_h::LOs d2oc = {dims2oclass[0],
                       dims2oclass[1],
                       dims2oclass[2],
                       dims2oclass[3]};
  Omega_h::Write<int> egCounts_d(3);
  auto egEnts_d = OhWriteEgo(3);
  auto egBody_d = OhWriteEgo(1);
  auto getTopo = OMEGA_H_LAMBDA(int) {
    printf("cuda eg_getTopo\n");
    ego model_geom;
    int model_oclass;
    int model_mtype;
    int* body_senses;
    //needed outputs
    int nbodies_local;
    ego* bodies_local;
    printf("eg_getTopo 0.1\n");
    CALL(EG_getTopology(egModel, &model_geom, &model_oclass, &model_mtype,
        nullptr, &nbodies_local, &bodies_local, &body_senses));
    printf("nbodies_local %d\n", nbodies_local);
    assert(nbodies_local == 1);
    egBody_d[0] = egoToGo(bodies_local[0]);
    printf("device body %p\n", bodies_local[0]);
    for (int i = 0; i < 3; ++i) {
      printf("d2oc[%d] %d\n", i, d2oc[i]);
      int counts;
      ego* ents;
      CALL(EG_getBodyTopos(bodies_local[0], nullptr, d2oc[i], &counts, &ents));
      egCounts_d[i] = counts;
      egEnts_d[i] = egoPtrToGo(ents);
      printf("device %d count %d ents %p\n",
          i, egCounts_d[i], egEnts_d[i]);
    }
    printf("eg_getTopo 0.3\n");
  };
  parallel_for(1, getTopo, "getEgadsTopo");
  assert(cudaSuccess == cudaDeviceSynchronize());
  const auto egEnts = Omega_h::HostRead<Omega_h::GO>(egEnts_d);
  const auto egCounts = Omega_h::HostRead<int>(egCounts_d);
  printf("created reads\n");
  for (int i = 0; i < 3; ++i) {
    eg->counts[i] = egCounts[i];
    eg->entities[i] = goToEgoPtr(egEnts[i]);
    //store the pointer to the array of egads faces
    eg->entities_d[i] = goToEgoPtr(egEnts[i]);
    printf("host %d count %d ents %p\n", i, eg->counts[i], eg->entities[i]);
  }
  printf("3.0\n");
  const auto egBody = Omega_h::HostRead<Omega_h::GO>(egBody_d);
  printf("host body %p\n", egBody[0]);
  eg->body = goToEgo(egBody[0]);
  printf("3.1\n");

  // preprocess edge and vertex adjacency to faces
  for (int i = 0; i < 2; ++i) {
    printf("3.11\n");
    Omega_h::Write<int> setSizes_d(eg->counts[i]);
    auto egBody = eg->body;
    auto egCounts = eg->counts[2];
    auto egEnts = eg->entities[2];
    //count the set sizes
    auto countIndexBody = OMEGA_H_LAMBDA(int) {
      for (int j = 0; j < egCounts; ++j) {
        auto face = egEnts[j];
        int nadj_ents;
        ego* adj_ents;
        CALL(EG_getBodyTopos(egBody, face, d2oc[i], &nadj_ents, &adj_ents));
        for (int k = 0; k < nadj_ents; ++k) {
          auto adj_ent = adj_ents[k];
          auto egIdx = EG_indexBodyTopo(egBody, adj_ent);
          assert(egIdx > 0); //egads error codes are <=0
          auto idx = egIdx-1;
          setSizes_d[idx]++;
        }
      }
    };
    parallel_for(1, countIndexBody, "getIndexBody");
    assert(cudaSuccess == cudaDeviceSynchronize());
    printf("3.12\n");

    const auto setSizes = Omega_h::HostRead<int>(setSizes_d);
    int totSize = 0;
    for(int j=0; j<setSizes.size(); j++) {
      totSize += setSizes[j];
    }
    auto idxs2adj_faces_d = OhWriteEgo(totSize);
    printf("3.13\n");

    Omega_h::Write<int> setCounts_d(eg->counts[i]);
    printf("3.2\n");
    //fill the sets
    auto getIndexBody = OMEGA_H_LAMBDA(int) {
      for (int j = 0; j < egCounts; ++j) {
        auto face = egEnts[j];
        int nadj_ents;
        ego* adj_ents;
        CALL(EG_getBodyTopos(egBody, face, d2oc[i], &nadj_ents, &adj_ents));
        for (int k = 0; k < nadj_ents; ++k) {
          auto adj_ent = adj_ents[k];
          auto egIdx = EG_indexBodyTopo(egBody, adj_ent);
          assert(egIdx > 0);
          auto idx = egIdx-1;
          auto ohIdx = setCounts_d[idx]++;
          idxs2adj_faces_d[ohIdx] = egoToGo(face);
        }
      }
    };
    parallel_for(1, getIndexBody, "getIndexBody");
    assert(cudaSuccess == cudaDeviceSynchronize());

    printf("3.3\n");
    std::vector<std::set<ego>> idxs2adj_faces(eg->counts[i]);
    //copy array into vector of sets
    const auto idxs2adj_faces_h = Omega_h::HostRead<Omega_h::GO>(idxs2adj_faces_d);
    int idx = 0;
    for(int j = 0; j < setSizes.size(); j++) {
      for(int k = 0; k < setSizes[j]; k++) {
        auto face = goToEgo(idxs2adj_faces_h[idx++]);
        idxs2adj_faces[j].insert(face);
      }
    }

    printf("3.4\n");
    //copy each device entity pointer to the host
    auto entPtrs_d = OhWriteEgo(eg->counts[i]);
    auto copyDevPtrs = OMEGA_H_LAMBDA(int j) {
      ego* ents = goToEgoPtr(egEnts_d[i]);
      auto entPtr = ents[j];
      entPtrs_d[j] = egoToGo(entPtr);
    };
    parallel_for(eg->counts[i], copyDevPtrs, "copyDevPtrs");
    assert(cudaSuccess == cudaDeviceSynchronize());
    auto entPtrs_h = Omega_h::HostRead<Omega_h::GO>(entPtrs_d);
    eg->entities[i] = new ego[eg->counts[i]];
    for(int j=0; j < eg->counts[i]; j++) {
      eg->entities[i][j] = goToEgo(entPtrs_h[j]);
    }
    
    printf("3.5\n");
    for (int j = 0; j < eg->counts[i]; ++j) {
      auto adj_faces = idxs2adj_faces[j];
      // HACK!: we have a really insane CAD model with nonsensical topology.
      // this essentially manifests as edges that are adjacent to only one
      // model face.
      // we actually want to just ignore these edges, so we won't create
      // classifier entries for them.
      if (adj_faces.size() == 1) continue;
      eg->classifier[adj_faces] = eg->entities[i][j];
    }
    printf("3.6\n");
  } //done loop over vertices and edges

  //copy each device face pointer to the host
  const int fdim = 2;
  auto entPtrs_d = OhWriteEgo(eg->counts[fdim]);
  auto copyDevPtrs = OMEGA_H_LAMBDA(int j) {
    ego* ents = goToEgoPtr(egEnts_d[fdim]);
    auto entPtr = ents[j];
    entPtrs_d[j] = egoToGo(entPtr);
  };
  parallel_for(eg->counts[fdim], copyDevPtrs, "copyDevPtrs");
  assert(cudaSuccess == cudaDeviceSynchronize());
  auto entPtrs_h = Omega_h::HostRead<Omega_h::GO>(entPtrs_d);
  eg->entities[fdim] = new ego[eg->counts[fdim]];
  for(int j=0; j < eg->counts[fdim]; j++) {
    eg->entities[fdim][j] = goToEgo(entPtrs_h[j]);
  }
  printf("3.7\n");

  //set struct device pointer for entity count
  eg->counts_d = LOs{eg->counts[0], eg->counts[1], eg->counts[2]};

  return eg;
}

static int get_dim(ego e) {
  ego ref;
  int oclass;
  int mtype;
  int nchild;
  ego* children;
  int* senses;
  CALL(EG_getTopology(
      e, &ref, &oclass, &mtype, nullptr, &nchild, &children, &senses));
  for (int i = 0; i <= 3; ++i)
    if (dims2oclass[i] == oclass) return i;
  return -1;
}

void egads_lite_classify(Egads* eg, int nadj_faces, int const adj_face_ids[],
    int* class_dim, int* class_id) {
  std::set<ego> uniq_adj_faces;
  for (int i = 0; i < nadj_faces; ++i) {
    const auto adjFaceId = adj_face_ids[i]-1;
    auto adj_face = eg->entities[2][adjFaceId];
    uniq_adj_faces.insert(adj_face);
  }
  auto it = eg->classifier.find(uniq_adj_faces);
  if (it != eg->classifier.end()) {
    auto ent = it->second;
    Omega_h::LOs d2oc = {dims2oclass[0],
                         dims2oclass[1],
                         dims2oclass[2],
                         dims2oclass[3]};
    Omega_h::Write<int> classDimAndId_d(2);
    auto egBody = eg->body;
    auto getEntClass = OMEGA_H_LAMBDA(int) {
      //get model entity dimension
      ego ref;
      int oclass;
      int mtype;
      int nchild;
      ego* children;
      int* senses;
      CALL(EG_getTopology(ent, &ref, &oclass, &mtype, nullptr, &nchild, &children, &senses));
      classDimAndId_d[0] = -1;
      for (int i = 0; i <= 3; ++i) {
        if (d2oc[i] == oclass) {
          classDimAndId_d[0] = i;
          break;
        }
      }
      //get model entity id
      auto egIdx = EG_indexBodyTopo(egBody, ent);
      assert(egIdx > 0);
      classDimAndId_d[1] = egIdx;
    };
    parallel_for(1, getEntClass, "getEntClass");
    assert(cudaSuccess == cudaDeviceSynchronize());
    auto classDimAndId = Omega_h::HostRead<int>(classDimAndId_d);
    //set host vars
    *class_dim = classDimAndId[0];
    *class_id = classDimAndId[1];
  }
}

void egads_lite_free(Egads* eg) {
  for (int i = 0; i < 3; ++i) {
    EG_free(eg->entities[i]); //TODO will this work with new ego[...] ?
  }
  CALL(EG_deleteObject(eg->model));
  CALL(EG_close(eg->context));
  delete eg;
}

void egads_lite_reclassify(Mesh* mesh, Egads* eg) {
  OMEGA_H_CHECK(mesh->dim() == 3);
  auto face_class_dims = mesh->get_array<I8>(FACE, "class_dim");
  auto face_class_ids = mesh->get_array<ClassId>(FACE, "class_id");
  for (Int dim = 0; dim < 2; ++dim) {
    auto ents2faces = mesh->ask_up(dim, FACE);
    auto adj_class_dims = read(unmap(ents2faces.ab2b, face_class_dims, 1));
    auto keep_edges = each_eq_to(adj_class_dims, I8(2));
    auto ents2eq_faces = filter_graph_edges(ents2faces, keep_edges);
    auto adj_eq_face_ids = unmap(ents2eq_faces.ab2b, face_class_ids, 1);
    auto host_a2ab = HostRead<LO>(ents2eq_faces.a2ab);
    auto host_face_ids = HostRead<LO>(adj_eq_face_ids);
    auto class_dims = mesh->get_array<I8>(dim, "class_dim");
    auto class_ids = mesh->get_array<ClassId>(dim, "class_id");
    auto host_class_dims = HostWrite<I8>(deep_copy(class_dims));
    auto host_class_ids = HostWrite<LO>(deep_copy(class_ids));
    for (LO i = 0; i < mesh->nents(dim); ++i) {
      auto b = host_a2ab[i];
      auto e = host_a2ab[i + 1];
      Int class_dim = host_class_dims[i];
      LO class_id = host_class_ids[i];
      egads_lite_classify(
          eg, e - b, host_face_ids.data() + b, &class_dim, &class_id);
      host_class_dims[i] = I8(class_dim);
      host_class_ids[i] = class_id;
    }
    class_dims = Read<I8>(host_class_dims.write());
    class_ids = Read<LO>(host_class_ids.write());
    mesh->set_tag(dim, "class_id", class_ids);
    mesh->set_tag(dim, "class_dim", class_dims);
  }
}

OMEGA_H_INLINE Vector<3> get_closest_point(ego g, Vector<3> in) {
  Vector<2> ignored;
  Vector<3> out = in;
  CALL(EG_invEvaluate(g, in.data(), ignored.data(), out.data()));
  return out;
}

Reals egads_lite_get_snap_warp(Mesh* mesh, Egads* eg, bool verbose) {
  fprintf(stderr, "numverts %d\n", mesh->nverts());
  //Omega_h::vtk::write_parallel("preWarp", mesh, mesh->dim());
  OMEGA_H_CHECK(mesh->dim() == 3);
  if (verbose) std::cout << "Querying closest points for surface vertices...\n";
  auto t0 = now();
  auto class_dims = mesh->get_array<I8>(VERT, "class_dim");
  auto class_ids = mesh->get_array<ClassId>(VERT, "class_id");
  auto coords = mesh->coords();
  auto closePts = Write<Real>(mesh->nverts() * 3);
  auto warp = Write<Real>(mesh->nverts() * 3);
  GOs egEnts_d{egoPtrToGo(eg->entities_d[0]),
               egoPtrToGo(eg->entities_d[1]),
               egoPtrToGo(eg->entities_d[2])};
  auto egCounts_d = eg->counts_d;
  auto egBody_d = eg->body;
  printf("vtx,class_id,class_dim,isEdge,isFace,ptX,ptY,ptZ,clPtX,clPtY,clPtZ,warpX,warpY,warpZ\n");
  auto calc_warp = OMEGA_H_LAMBDA(LO i) {
    auto a = get_vector<3>(coords, i);
    Int class_dim = class_dims[i];
    OMEGA_H_CHECK(class_dim >= 0);
    OMEGA_H_CHECK(class_dim <= 3);
    auto d = vector_3(0, 0, 0);
    auto clPt = vector_3(0, 0, 0);
    if (0 < class_dim && class_dim < 3) { //edges and faces only
      auto index = class_ids[i] - 1;
      OMEGA_H_CHECK(index >= 0);
      OMEGA_H_CHECK(index < egCounts_d[class_dim]);
      auto ents = goToEgoPtr(egEnts_d[class_dim]);
      auto g = ents[index];
      auto index2 = EG_indexBodyTopo(egBody_d, g);
      assert(index2 > 0);
      OMEGA_H_CHECK(index2 == index + 1);
      auto b = get_closest_point(g, a);
      clPt = b;
      d = b - a;
      {
      int isEdge = (g->oclass == EGADS_EDGE);
      int isFace = (g->oclass == EGADS_FACE);
      printf("%d,%d,%d,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",
          i, index2, class_dim, isEdge, isFace, a[0], a[1], a[2],
          clPt[0], clPt[1], clPt[2],
          d[0],d[1],d[2]);
      }
    }
    set_vector(warp, i, d);
    set_vector(closePts, i, clPt);
  };
  parallel_for(mesh->nverts(), std::move(calc_warp), "calc_warp"); 
  assert(cudaSuccess == cudaDeviceSynchronize());
  auto t1 = now();
  mesh->add_tag(0, "warpVec", 3, read(warp));
  mesh->add_tag(0, "closePts", 3, read(closePts));
  //Omega_h::vtk::write_parallel("warpVec", mesh, mesh->dim());
  if (verbose) {
    std::cout << "Querying closest points for surface vertices took "
              << (t1 - t0) << " seconds\n";
  }
  return warp;
}

void checkCudaError(int line) {
#ifdef __NVCC__
  cudaError_t code = cudaDeviceSynchronize();
  const char * errorMessage = cudaGetErrorString(code);
  if( code != cudaSuccess ) {
    fprintf(stderr, "CUDA error on line %d Error code: %d (%s)\n", line, code, errorMessage);
  }
  assert(code == cudaSuccess);
#endif
}

void hackClassification(Omega_h::Mesh* mesh) {
  fprintf(stderr, "hacking classification\n");
  OMEGA_H_CHECK(mesh->dim() == 3);
  auto vtx_class_dims = mesh->get_array<Omega_h::I8>(Omega_h::VERT, "class_dim");
  auto vtx_class_ids_r = mesh->get_array<Omega_h::ClassId>(Omega_h::VERT, "class_id");
  auto vtx_class_ids_w = Omega_h::deep_copy(vtx_class_ids_r, "vtxClassIds_w");
  auto setVtxClass = OMEGA_H_LAMBDA(int i) {
    if(vtx_class_dims[i] == 1 && vtx_class_ids_w[i] == 1) {
      printf("vtx %i reclassified\n",i);
      vtx_class_ids_w[i] = 7;
    }
  };
  Omega_h::parallel_for(mesh->nents(0), setVtxClass, "setVtxClass");
  fprintf(stderr, "done hacking vtx classification\n");
  mesh->set_tag(0, "class_id", Omega_h::read(vtx_class_ids_w));

  auto edge_class_dims = mesh->get_array<Omega_h::I8>(Omega_h::EDGE, "class_dim");
  auto edge_class_ids_r = mesh->get_array<Omega_h::ClassId>(Omega_h::EDGE, "class_id");
  auto edge_class_ids_w = Omega_h::deep_copy(edge_class_ids_r, "edgeClassIds_w");
  auto setEdgeClass = OMEGA_H_LAMBDA(int i) {
    if(edge_class_dims[i] == 1 && edge_class_ids_w[i] == 1) {
      printf("edge %i reclassified\n",i);
      edge_class_ids_w[i] = 7;
    }
  };
  Omega_h::parallel_for(mesh->nents(1), setEdgeClass, "setEdgeClass");
  fprintf(stderr, "done hacking edge classification\n");
  mesh->set_tag(1, "class_id", Omega_h::read(edge_class_ids_w));

}

void setCudaStackSz() {
  size_t stackLimit;
  cuCtxGetLimit(&stackLimit, CU_LIMIT_STACK_SIZE);
  checkCudaError(__LINE__);
  printf("original stack limit %d\n", stackLimit);
  stackLimit=8*1024;
  cuCtxSetLimit(CU_LIMIT_STACK_SIZE,stackLimit);
  checkCudaError(__LINE__);
  cuCtxGetLimit(&stackLimit, CU_LIMIT_STACK_SIZE);
  checkCudaError(__LINE__);
  printf("new stack limit %d\n", stackLimit);
  printf("stack limit %d\n", stackLimit);
}

}  // namespace Omega_h
