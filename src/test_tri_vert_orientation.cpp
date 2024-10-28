#include <Omega_h_macros.h>

#include <Omega_h_adj.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_mesh.hpp>
#include <cstdio>
#include <Omega_h_file.hpp>

bool is_ccw_oriented(Omega_h::Mesh& mesh) {
  OMEGA_H_CHECK(mesh.dim() == 2);
  bool ccw;
  // find mesh orientation: node number orientation of faces
  // clockwise or counter clockwise: 1 for clockwise, -1 for counter clockwise
  const auto& face2nodes = mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;
  const auto& coords = mesh.coords();
  Omega_h::Write<Omega_h::LO> ccw_face(mesh.nfaces(), 0, "ccw_face");

  auto find_ccw_face = OMEGA_H_LAMBDA(const Omega_h::LO face) {
    const auto faceVerts = Omega_h::gather_verts<3>(face2nodes, face);
    const auto faceCoords = Omega_h::gather_vectors<3, 2>(coords, faceVerts);
    // calculate using (y2 - y1) * (x3 - x2) - (y3 - y2) * (x2 - x1)
    const bool cw = (faceCoords[1][1] - faceCoords[0][1]) *
                            (faceCoords[2][0] - faceCoords[1][0]) -
                        (faceCoords[2][1] - faceCoords[1][1]) *
                            (faceCoords[1][0] - faceCoords[0][0]) >
                    0;
    ccw_face[face] = (!cw) ? 1 : -1;
  };
  Omega_h::parallel_for(mesh.nfaces(), find_ccw_face, "find_ccw_face");
  auto ccw_check = OMEGA_H_LAMBDA(const Omega_h::LO face) {
    // check all the faces are either clockwise or counter clockwise but same
    // orientation
    OMEGA_H_CHECK((ccw_face[face] == ccw_face[0]) && (ccw_face[face] != 0));
  };
  Omega_h::parallel_for(mesh.nfaces(), ccw_check, "ccw_check");

  auto ccw_check_host = Omega_h::HostRead<Omega_h::LO>(ccw_face);
  ccw = ccw_check_host[0] == 1;
  printf(
      "**** ORIENTATION: %s ****\n", (ccw) ? "COUNTER CLOCKWISE" : "CLOCKWISE");

  return ccw;
}

int main(int argc, char** argv) {
    Omega_h::Library library(&argc, &argv);
    Omega_h::Mesh mesh(&library);
    mesh = Omega_h::binary::read(argv[1], library.self());
    bool ccw = is_ccw_oriented(mesh);

    if (!ccw) {
        std::fprintf(stderr, "ERROR: Mesh is not counter clockwise oriented\n");
        return 1;
    }

    return 0;
}