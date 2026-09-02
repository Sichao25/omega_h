#include <Omega_h_adapt.hpp>
#include <Omega_h_cmdline.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_for.hpp>

#ifdef OMEGA_H_USE_EGADSLITE
#include <Omega_h_egads_lite.hpp>
#endif

#include <cuda.h>

#include <iostream>

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

static void compute_implied_metric(Omega_h::Mesh* mesh) {
  auto metrics = Omega_h::get_implied_metrics(mesh);
  metrics = Omega_h::limit_metric_gradation(mesh, metrics, 1.0);
  mesh->add_tag(
      Omega_h::VERT, "metric", Omega_h::symm_ncomps(mesh->dim()), metrics);
}

static void compute_target_metric(Omega_h::Mesh* mesh) {
  auto metric = Omega_h::diagonal(Omega_h::metric_eigenvalues_from_lengths(
      Omega_h::vector_3(0.1, 0.1, 0.1)));
  auto metrics = Omega_h::repeat_symm(mesh->nverts(), metric);
  mesh->add_tag(Omega_h::VERT, "target_metric",
      Omega_h::symm_ncomps(mesh->dim()), metrics);
}

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  setCudaStackSz();
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("mesh_in.meshb");
  cmdline.add_arg<std::string>("mesh_out.meshb");
#ifdef OMEGA_H_USE_EGADSLITE
  auto& model_flag = cmdline.add_flag("--model", "optional EGADS model");
  model_flag.add_arg<std::string>("model.step");
#endif
  auto& viz_flag = cmdline.add_flag("--viz", "optional VTK progress log");
  viz_flag.add_arg<std::string>("path_vtk");
  if (!cmdline.parse_final(lib.world(), &argc, argv)) return -1;
  auto path_in = cmdline.get<std::string>("mesh_in.meshb");
  auto path_out = cmdline.get<std::string>("mesh_out.meshb");
  Omega_h::Mesh mesh(&lib);
  std::cout << "reading in " << path_in << '\n';
  Omega_h::meshb::read(&mesh, path_in);
  std::cout << "computing metric tags\n";
  compute_implied_metric(&mesh);
  compute_target_metric(&mesh);
  std::cout << "computing minimum quality\n";
  Omega_h::AdaptOpts opts(&mesh);
#ifdef OMEGA_H_USE_EGADSLITE
  auto has_model = cmdline.parsed("--model");
  if (has_model) {
    auto model_path = cmdline.get<std::string>("--model", "model.step");
    std::cout << "reading in " << model_path << '\n';
    auto eg = Omega_h::egads_lite_load(model_path);
    Omega_h::egads_lite_reclassify(&mesh, eg);
    opts.egads_model = eg;
    hackClassification(&mesh); //there are problems...
    Omega_h::vtk::write_parallel("postHack", &mesh, mesh.dim());
  }
#endif
  fprintf(stderr, "numverts %d\n", mesh.nents(0));
  auto has_viz = cmdline.parsed("--viz");
  Omega_h::vtk::Writer writer;
  if (has_viz) {
    auto viz_path = cmdline.get<std::string>("--viz", "path_vtk");
    writer = Omega_h::vtk::Writer(viz_path, &mesh);
    writer.write();
  }
  opts.verbosity = Omega_h::EXTRA_STATS;
  opts.max_length_allowed = opts.max_length_desired * 2.0;
  while (Omega_h::approach_metric(&mesh, opts)) {
    Omega_h::adapt(&mesh, opts);
    if (has_viz) writer.write();
  }
  std::cout << "writing out " << path_out << '\n';
  mesh.remove_tag(Omega_h::VERT, "metric");
  Omega_h::meshb::write(&mesh, path_out);
#ifdef OMEGA_H_USE_EGADSLITE
  if (has_model) {
    Omega_h::egads_lite_free(opts.egads_model);
  }
#endif
}
