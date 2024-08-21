#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_for.hpp> //parallel_for
#include <Omega_h_array_ops.hpp> //get_min
#include <Omega_h_int_scan.hpp> //offset_scan
#include <Omega_h_simplex.hpp> //simplex_degree
#include <Kokkos_NestedSort.hpp> //sort_team
#include <Omega_h_build.hpp> //build_box

#include <cstdlib>

using namespace Omega_h;

int verbose = 0;

void render(Mesh& m, Graph patches, std::string suffix) {
  const auto num_patches = patches.nnodes();
  auto offsets = patches.a2ab;
  Write<LO> degree(num_patches);
  parallel_for(num_patches, OMEGA_H_LAMBDA(LO i) {
    degree[i] = offsets[i+1]-offsets[i];
  });
  m.add_tag(VERT, "patchDegree", 1, read(degree));
  std::string name = "patch_" + suffix + ".vtk";
  vtk::write_parallel(name, &m, m.dim());
}

void writeGraph(Graph g, std::string name="") {
  if(verbose<2) return;
  HostRead offsets(g.a2ab);
  HostRead values(g.ab2b);
  std::cout << "== " << name << " ==\n";
  for(int node = 0; node < g.nnodes(); node++) {
    std::cout << node << ": ";
     for(int edge = offsets[node]; edge < offsets[node+1]; edge++) {
       std::cout << values[edge] << " ";
     }
     std::cout << "\n";
  }
  if(verbose<3) return;
  std::cout << "offsets = {";
  for(int i = 0; i < offsets.size(); i++) {
    std::cout << offsets[i];
    if(i!=offsets.size()-1) std::cout << ",";
  }
  std::cout << "};\n";
  std::cout << "values = {";
  for(int i = 0; i < values.size(); i++) {
    std::cout << values[i];
    if(i!=values.size()-1) std::cout << ",";
  }
  std::cout << "};\n";
}

[[nodiscard]] Graph adj_segment_sort(Graph& g) {
  using ExecSpace = Kokkos::DefaultExecutionSpace;
  using TeamPol = Kokkos::TeamPolicy<ExecSpace>;
  using TeamMem = typename TeamPol::member_type;
  auto offsets = g.a2ab;
  auto elms_r = g.ab2b.view(); //read only
  Kokkos::View<LO*, ExecSpace> elms("elms", elms_r.size());
  Kokkos::deep_copy(elms, elms_r);
  auto segment_sort = KOKKOS_LAMBDA(const TeamMem& t) {
    auto i = t.league_rank();
    auto patch = Kokkos::subview(elms, Kokkos::make_pair(offsets[i], offsets[i+1]));
    Kokkos::Experimental::sort_team(t, patch);
  };
  Kokkos::parallel_for(TeamPol(g.nnodes(), Kokkos::AUTO()), segment_sort);
  return Graph(offsets,Write<LO>(elms));
}

[[nodiscard]] Graph remove_duplicate_edges(Graph g) {
  auto offsets = g.a2ab;
  auto values = g.ab2b;
  Write<I8> keep(g.nedges());
  auto markDups = OMEGA_H_LAMBDA(LO i) {
    keep[offsets[i]] = 1; //always keep the first edge in the segment
    for(int j=offsets[i]+1; j<offsets[i+1]; j++) {
      keep[j] = (values[j-1] != values[j]);
    }
  };
  parallel_for(g.nnodes(), markDups);
  auto filtered = filter_graph_edges(g,keep);
  return filtered;
}

[[nodiscard]] Read<I8> patchSufficient(Graph patches, Int minPatchSize) {
  const auto num_patches = patches.nnodes();
  auto offsets = patches.a2ab;
  Write<I8> done(num_patches);
  parallel_for(num_patches, OMEGA_H_LAMBDA(LO i) {
    done[i] = ((offsets[i+1]-offsets[i]) >= minPatchSize);
  });
  return read(done);
}

/**
 * \brief expand the patches
 * \param m (in) mesh of simplices
 * \param patches (in) graph of key entities to elements
 * \param adjElms (in) second order element-to-element adjacencies
 *                     used for expansion
 * \return an expanded graph from key entities to elements
*/
//FIXME ideally, this used the Omega_h_map and Omega_h_graph functions
[[nodiscard]] Graph expandPatches(Mesh& m, Graph patches, Graph adjElms, Read<I8> patchDone) {
  auto adjElms_offsets = adjElms.a2ab;
  auto adjElms_elms = adjElms.ab2b;
  const auto num_patches = patches.nnodes();
  auto patch_offsets = patches.a2ab;
  auto patch_elms = patches.ab2b;
  Write<LO> degree(num_patches);
  parallel_for(num_patches, OMEGA_H_LAMBDA(LO patch) {
    degree[patch] = patch_offsets[patch+1] - patch_offsets[patch];
    if(!patchDone[patch]) {
      for(int j=patch_offsets[patch]; j<patch_offsets[patch+1]; j++) {
        auto elm = patch_elms[j];
        degree[patch] += adjElms_offsets[elm+1]-adjElms_offsets[elm]; //counts duplicates
      }
    }
  });
  auto patchExpDup_offsets = offset_scan(read(degree));
  Write<LO> patchExpDup_elms(patchExpDup_offsets.last());
  parallel_for(num_patches, OMEGA_H_LAMBDA(LO patch) {
    auto idx = patchExpDup_offsets[patch];
    for(int j=patch_offsets[patch]; j<patch_offsets[patch+1]; j++) {
      patchExpDup_elms[idx++] = patch_elms[j];
      if(!patchDone[patch]) {
        auto elm = patch_elms[j];
        for(int k=adjElms_offsets[elm]; k<adjElms_offsets[elm+1]; k++) {
          patchExpDup_elms[idx++] = adjElms_elms[k];
        }
      }
    }
  });
  Graph patchExpDup(patchExpDup_offsets,patchExpDup_elms);
  auto sorted = adj_segment_sort(patchExpDup);
  auto dedup = remove_duplicate_edges(sorted);
  writeGraph(dedup, "dedup");
  return dedup;
}

/**
 * \brief form a patch of at least minPatchSize elements surrounding each entity
 *        of dimension keyDim
 * \remark the patch is expanded via 2nd order adjacencies using meshDim-1 as
 *         the bridge entity (e.g., faces for 3d, edges for 2d)
 * \param m (in) mesh of simplices
 * \param keyDim (in) the dimension of mesh entities that the patches are
 *        created around
 * \param minPatchSize (in) the minimum number of elements in each patch
 * \return a graph whose source nodes are the entities of keyDim dimension, and
 *         edges are connecting to elements in the patch
 */
[[nodiscard]] Graph formPatches(Mesh& m, LO keyDim, Int minPatchSize) {
  OMEGA_H_CHECK(keyDim >= 0 && keyDim < m.dim());
  OMEGA_H_CHECK(minPatchSize > 0);
  auto patches = m.ask_up(VERT,m.dim());
  writeGraph(patches);
  render(m,patches,"init");
  auto patchDone = patchSufficient(patches, minPatchSize);
  if( get_min(patchDone) == 1 )
    return patches;
  auto adjElms = m.ask_dual();
  writeGraph(adjElms, "adjElms");
  for(Int iter = 0; iter < 10; iter++) {
    if(verbose>=2) std::cout << iter << " expanding patch\n";
    patches = expandPatches(m, patches, adjElms, patchDone);
    render(m,patches,std::to_string(iter));
    patchDone = patchSufficient(patches, minPatchSize);
    if( get_min(patchDone) == 1 ) {
      if(verbose>=1) std::cout << "iterations: " << iter << "\n";
      return patches;
    }
  }
  assert(false);
  return Graph();
}

void testGraphSort() {
  {
    Graph g({0,6},{1,3,2,4,4,3});
    auto res = adj_segment_sort(g);
    Graph expected({0,6},{1,2,3,3,4,4});
    OMEGA_H_CHECK(res == expected);
  }
  {
    Graph g({0,4,6},{4,3,2,2,1,0});
    auto res = adj_segment_sort(g);
    Graph expected({0,4,6},{2,2,3,4,0,1});
    OMEGA_H_CHECK(res == expected);
  }
}

void testGraphDuplicateRemoval() {
  {
     Graph g({0,2,6},{1,1,2,3,7,7});
     auto res = remove_duplicate_edges(g);
     Graph expected({0,1,4},{1,2,3,7});
     OMEGA_H_CHECK(res == expected);
  }
  {
     Graph g({0,6},{1,1,1,1,1,1});
     auto res = remove_duplicate_edges(g);
     Graph expected({0,1},{1});
     OMEGA_H_CHECK(res == expected);
  }
}

void test2x2(Omega_h::Library& lib) {
  auto world = lib.world();
  const auto x = 2.0;
  const auto y = 2.0;
  const auto z = 0.0;
  const auto nx = 2;
  const auto ny = 2;
  const auto nz = 0;
  const auto symmetric = false;
  auto mesh = Omega_h::build_box(world, OMEGA_H_SIMPLEX, x, y, z, nx, ny, nz, symmetric);
  const auto minPatchSize = 3;
  auto patches = formPatches(mesh, VERT, minPatchSize);
  Graph expected(
   {0,4,7,11,17,20,24,27,30,34},
   {0,1,2,6,1,2,4,1,2,4,5,0,1,2,3,5,6,1,4,5,1,3,5,6,3,6,7,0,6,7,0,3,6,7});
  OMEGA_H_CHECK(patches == expected);
}

void test1x5(Omega_h::Library& lib) {
  auto world = lib.world();
  const auto x = 1.0;
  const auto y = 5.0;
  const auto z = 0.0;
  const auto nx = 1;
  const auto ny = 5;
  const auto nz = 0;
  const auto symmetric = false;
  auto mesh = Omega_h::build_box(world, OMEGA_H_SIMPLEX, x, y, z, nx, ny, nz, symmetric);
  {
    const auto minPatchSize = 3;
    auto patches = formPatches(mesh, VERT, minPatchSize);
    Graph expected(
      {0,3,6,9,12,15,18,21,24,27,30,33,36},
      {0,1,2,1,2,3,0,1,2,0,1,2,1,3,4,3,4,6,5,6,9,4,5,6,7,8,9,7,8,9,7,8,9,5,7,9});
    OMEGA_H_CHECK(patches == expected);
  }
  {
    const auto minPatchSize = 4;
    auto patches = formPatches(mesh, VERT, minPatchSize);
    Graph expected(
       {0,4,9,13,17,22,27,32,37,41,45,49,54},
       {0,1,2,3,0,1,2,3,4,0,1,2,3,0,1,2,3,1,2,3,4,6,1,3,
        4,5,6,4,5,6,7,9,3,4,5,6,9,5,7,8,9,5,7,8,9,5,7,8,9,5,6,7,8,9});
    OMEGA_H_CHECK(patches == expected);
  }
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  testGraphSort();
  testGraphDuplicateRemoval();
  OMEGA_H_CHECK(argc == 4);
  verbose = std::stoi(argv[3]);
  test2x2(lib);
  test1x5(lib);
  Mesh mesh(&lib);
  binary::read(argv[1], lib.world(), &mesh);
  auto patches = formPatches(mesh, VERT, 3);
  vtk::write_parallel(argv[2], &mesh, mesh.dim());
}
