#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_for.hpp> //parallel_for
#include <Omega_h_array_ops.hpp> //get_min
#include <Kokkos_NestedSort.hpp> //sort_team

#include <cstdlib>

using namespace Omega_h;

[[nodiscard]] Graph getElmToElm2ndOrderAdj(Mesh& m, Int bridgeDim) {
  const auto elmDim = m.dim();
  OMEGA_H_CHECK(bridgeDim >= 0 && bridgeDim < elmDim);
  auto elm2Bridge = m.ask_down(elmDim, bridgeDim);
  auto e2bOffsets = elm2Bridge.a2ab;
  auto e2bValues = elm2Bridge.ab2b;
  //get bridge-to-elm
  auto bridge2Elm = m.ask_up(bridgeDim, elmDim);
  //traverse elm-to-bridge and bridge-to-elm to form elm-to-elm
  //are there map and graph functions to compute this?
  //first count how many elms there are per elm
  Write<LO> e2eDegree(m.nelems());
  auto b2eOffsets = bridge2Elm.a2ab;
  parallel_for(elm, OMEGA_H_LAMBDA(LO i) {
    //HERE loop over e2bOffsets
      //get bridge
      e2eDegree[i] = b2eOffsets[bridge+1]-b2eOffsets[bridge]; //counts duplicates!
  });

  return Graph();
}

[[nodiscard]] bool patchSufficient(Graph patches, Int minPatchSize) {
  const auto num_patches = patches.a2ab.size()-1;
  auto offsets = patches.a2ab;
  Write<LO> degree(num_patches);
  parallel_for(num_patches, OMEGA_H_LAMBDA(LO i) {
    degree[i] = offsets[i+1]-offsets[i];
  });
  auto minDegree = get_min(read(degree));
  if( minDegree < minPatchSize)
    return false;
  else 
    return true;
}

[[nodiscard]] Graph adj_segment_sort(Graph& g) {
  using ExecSpace = Kokkos::DefaultExecutionSpace;
  using TeamPol = Kokkos::TeamPolicy<ExecSpace>;
  using TeamMem = typename TeamPol::member_type;
  auto offsets = g.a2ab;
  auto elms_r = g.ab2b.view(); //read only
  Kokkos::View<LO*, ExecSpace> elms("elms", elms.size());
  Kokkos::deep_copy(elms, elms_r);
  auto segment_sort = KOKKOS_LAMBDA(const TeamMem& t) {
    //Sort a row of A using the whole team.
    auto i = t.league_rank();
    auto patch = Kokkos::subview(elms, Kokkos::make_pair(offsets[i], offsets[i+1]));
    Kokkos::Experimental::sort_team(t, patch);
  };
  Kokkos::parallel_for(TeamPol(offsets.size()-1, Kokkos::AUTO()), segment_sort);
  return Graph(offsets,Write<LO>(elms));
}

[[nodiscard]] Graph remove_duplicate_edges(Graph g) {
  return Graph();
}

/**
 * \brief expand the patches
 * \param m (in) mesh of simplices
 * \param patches (in) graph of key entities to elements
 * \param bridgeDim (in) the entity dimension used for expansion via second
 *        order element-to-element adjacencies
 * \return an expanded graph from key entities to elements
*/ 
[[nodiscard]] Graph expandPatches(Mesh& m, Graph patches, Int bridgeDim) {
  OMEGA_H_CHECK(bridgeDim >= 0 && bridgeDim < m.dim());
  auto patch_elms = patches.ab2b;
  auto adjElms = getElmToElm2ndOrderAdj(m, bridgeDim);
  auto expanded = unmap_graph(patch_elms, adjElms);
  auto sorted = adj_segment_sort(expanded);
  auto dedup = remove_duplicate_edges(sorted);
  return dedup;
}

/**
 * \brief form a patch of at least minPatchSize elements surrounding each entity
 *        of dimension keyDim
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
  auto patches = Graph();
  for(Int bridgeDim = m.dim()-1; bridgeDim >= 0; bridgeDim--) {
    auto bridgePatches = expandPatches(m, patches, bridgeDim);
    if( patchSufficient(bridgePatches, minPatchSize) ) {
      return bridgePatches;
    }
  }
  assert(false);
  return Graph();
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  OMEGA_H_CHECK(argc == 3);
  Mesh mesh(&lib);
  binary::read(argv[1], lib.world(), &mesh);
  auto patches = formPatches(mesh, VERT, 3);
  vtk::write_parallel(argv[2], &mesh, mesh.dim());
}
