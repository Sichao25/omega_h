#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_for.hpp> //parallel_for
#include <Omega_h_array_ops.hpp> //get_min
#include <Omega_h_int_scan.hpp> //offset_scan
#include <Omega_h_simplex.hpp> //simplex_degree
#include <Kokkos_NestedSort.hpp> //sort_team

#include <cstdlib>

using namespace Omega_h;

[[nodiscard]] Graph adj_segment_sort(Graph& g) { //not tested
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

[[nodiscard]] Graph remove_duplicate_edges(Graph g) { //not tested
  auto offsets = g.a2ab;
  auto values = g.ab2b;
  Write<I8> keep(g.nedges());
  auto markDups = OMEGA_H_LAMBDA(LO i) {
    keep[offsets[i]] = 1;
    for(int j=offsets[i]+1; j<offsets[i+1]; j++) {
      keep[j] = (values[j-1] != values[j]);
    }
  };
  parallel_for(g.nnodes(), markDups);
  auto filtered = filter_graph_edges(g,keep);
  return filtered;
}

[[nodiscard]] Graph getElmToElm2ndOrderAdj(Mesh& m, Int bridgeDim) {
  const auto elmDim = m.dim();
  OMEGA_H_CHECK(bridgeDim >= 0 && bridgeDim < elmDim);
  if(bridgeDim == elmDim-1)
    return m.ask_dual();
  auto elm2Bridge = m.ask_down(elmDim, bridgeDim);
  auto e2bDegree = simplex_degree(elmDim, bridgeDim);
  auto e2bValues = elm2Bridge.ab2b;
  //get bridge-to-elm
  auto bridge2Elm = m.ask_up(bridgeDim, elmDim);
  auto b2eOffsets = bridge2Elm.a2ab;
  auto b2eValues = bridge2Elm.ab2b;
  //traverse elm-to-bridge and bridge-to-elm to form elm-to-elm
  //first count how many adj elms there are per elm
  Write<LO> e2eDegree(m.nelems(),0);
  parallel_for(m.nelems(), OMEGA_H_LAMBDA(LO i) {
    for(int j=i*e2bDegree; j<(i+1)*e2bDegree; j++) {
      auto bridge = e2bValues[j];
      e2eDegree[i] += b2eOffsets[bridge+1]-b2eOffsets[bridge];
    }
  });
  auto e2eOffsets = offset_scan(read(e2eDegree));
  Write<LO> e2eValues(e2eOffsets.last());
  parallel_for(m.nelems(), OMEGA_H_LAMBDA(LO i) {
    auto pos = e2eOffsets[i];
    for(int j=i*e2bDegree; j<(i+1)*e2bDegree; j++) {
      auto bridge = e2bValues[j];
      for(int k=b2eOffsets[bridge]; k<b2eOffsets[bridge+1]; k++) {
        assert(pos < e2eOffsets[i+1]);
        e2eValues[pos++] = b2eValues[k];
      }
    }
  });
  auto e2e_unsorted = Graph(e2eOffsets, read(e2eValues));
  auto e2e_dups = adj_segment_sort(e2e_unsorted);
  return remove_duplicate_edges(e2e_dups);
}

[[nodiscard]] bool patchSufficient(Graph patches, Int minPatchSize) {
  const auto num_patches = patches.nnodes();
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


/**
 * \brief expand the patches
 * \param m (in) mesh of simplices
 * \param patches (in) graph of key entities to elements
 * \param bridgeDim (in) the entity dimension used for expansion via second
 *        order element-to-element adjacencies
 * \return an expanded graph from key entities to elements
*/
//FIXME ideally, this used the Omega_h_map and Omega_h_graph functions
[[nodiscard]] Graph expandPatches(Mesh& m, Graph patches, Int bridgeDim) {
  OMEGA_H_CHECK(bridgeDim >= 0 && bridgeDim < m.dim());
  auto adjElms = getElmToElm2ndOrderAdj(m, bridgeDim);
  auto adjElms_offsets = adjElms.a2ab;
  auto adjElms_elms = adjElms.ab2b;
  const auto num_patches = patches.nnodes();
  auto patch_offsets = patches.a2ab;
  auto patch_elms = patches.ab2b;
  Write<LO> degree(num_patches);
  parallel_for(num_patches, OMEGA_H_LAMBDA(LO patch) {
    degree[patch] = 0;
    for(int j=patch_offsets[patch]; j<patch_offsets[patch+1]; j++) {
      auto elm = patch_elms[j];
      degree[patch] += adjElms_offsets[elm+1]-adjElms_offsets[elm]; //counts duplicates
    }
  });
  auto patchExpDup_offsets = offset_scan(read(degree));
  Write<LO> patchExpDup_elms(patchExpDup_offsets.last());
  parallel_for(num_patches, OMEGA_H_LAMBDA(LO patch) {
    auto idx = 0;
    for(int j=patch_offsets[patch]; j<patch_offsets[patch+1]; j++) {
      auto elm = patch_elms[j];
      for(int k=adjElms_offsets[elm]; k<adjElms_offsets[elm+1]; k++) {
        patchExpDup_elms[idx++] = adjElms_elms[k];
      }
    }
  });
  Graph patchExpDup(patchExpDup_offsets,patchExpDup_elms);
  auto sorted = adj_segment_sort(patchExpDup);
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
  auto patches = m.ask_up(VERT,m.dim());
  if( patchSufficient(patches, minPatchSize) ) {
    return patches;
  }
  for(Int bridgeDim = m.dim()-1; bridgeDim >= 0; bridgeDim--) {
    auto bridgePatches = expandPatches(m, patches, bridgeDim);
    if( patchSufficient(bridgePatches, minPatchSize) ) {
      return bridgePatches;
    }
  }
  assert(false); //fails here
  return Graph();
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  OMEGA_H_CHECK(argc == 3);
  Mesh mesh(&lib);
  binary::read(argv[1], lib.world(), &mesh);
  {
     Graph g({0,2,6},{1,1,2,3,7,7});
     auto res = remove_duplicate_edges(g);
     Graph expected({0,1,4},{1,2,3,7});
     OMEGA_H_CHECK(res == expected);
  }
     Graph g({0,6},{1,1,1,1,1,1});
     auto res = remove_duplicate_edges(g);
     Graph expected({0,1},{1});
     OMEGA_H_CHECK(res == expected);
  }
  auto patches = formPatches(mesh, VERT, 3);
  vtk::write_parallel(argv[2], &mesh, mesh.dim());
}
