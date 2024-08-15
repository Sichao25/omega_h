#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>

#include <cstdlib>

using namespace Omega_h;

/**
 * \brief expand the patches
 * \param m (in) mesh of simplices
 * \param patches (in) graph of key entities to elements
 * \param bridgeDim (in) the entity dimension used for expansion via second
 *        order element-to-element adjacencies
 * \return an expanded graph from key entities to elements
*/ 
Graph expandPatches(Mesh& m, Graph patches, Int bridgeDim) {
  OMEGA_H_CHECK(bridgeDim >= 0 && bridgeDim < m.dim());
  return Graph();
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
Graph formPatches(Mesh& m, LO keyDim, Int minPatchSize) {
  OMEGA_H_CHECK(keyDim >= 0 && keyDim < m.dim());
  OMEGA_H_CHECK(minPatchSize > 0);
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
