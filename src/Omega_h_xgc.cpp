#include "Omega_h_file.hpp"

#include <fstream>
#include <string>
#include <vector>

#include "Omega_h_build.hpp"
#include "Omega_h_mark.hpp"

namespace Omega_h {

// XGC mesh format based on https://xgc.pppl.gov/html/mesh_file_format.html
namespace xgc {

namespace {

void read_node_file(std::istream& stream, std::vector<Real>& coords,
    std::vector<Int>& boundary_markers) {
  Int nverts, dim, nattributes, nboundary_markers;
  stream >> nverts >> dim >> nattributes >> nboundary_markers;
  OMEGA_H_CHECK(nverts > 0);
  OMEGA_H_CHECK(dim == 2);
  OMEGA_H_CHECK(nattributes == 0);
  OMEGA_H_CHECK(nboundary_markers == 0 || nboundary_markers == 1);

  coords.resize(static_cast<std::size_t>(nverts) * 2);
  boundary_markers.resize(static_cast<std::size_t>(nverts), 0);

  for (Int i = 0; i < nverts; ++i) {
    Int index;
    Real r, z;
    stream >> index >> r >> z;
    OMEGA_H_CHECK(index == i + 1);
    coords[static_cast<std::size_t>(i * 2)] = r;
    coords[static_cast<std::size_t>(i * 2 + 1)] = z;
    if (nboundary_markers == 1) {
      Int bm;
      stream >> bm;
      boundary_markers[static_cast<std::size_t>(i)] = bm;
    }
  }
}

void read_ele_file(std::istream& stream, std::vector<LO>& conn) {
  Int ntris, nodes_per_tri, nattributes;
  stream >> ntris >> nodes_per_tri >> nattributes;
  OMEGA_H_CHECK(ntris > 0);
  OMEGA_H_CHECK(nodes_per_tri == 3);
  OMEGA_H_CHECK(nattributes == 0);

  conn.resize(static_cast<std::size_t>(ntris) * 3);

  for (Int i = 0; i < ntris; ++i) {
    Int index;
    stream >> index;
    for (Int j = 0; j < 3; ++j) {
      Int node;
      stream >> node;
      conn[static_cast<std::size_t>(i * 3 + j)] = node - 1;
    }
  }
}

}  // end anonymous namespace

Mesh read(filesystem::path const& basename, CommPtr comm) {
  auto node_path = basename.string() + ".node";
  auto ele_path = basename.string() + ".ele";

  auto mesh = Mesh(comm->library());
  if (comm->rank() == 0) {
    std::ifstream node_stream(node_path);
    if (!node_stream.is_open()) {
      Omega_h_fail("couldn't open \"%s\"\n", node_path.c_str());
    }
    std::ifstream ele_stream(ele_path);
    if (!ele_stream.is_open()) {
      Omega_h_fail("couldn't open \"%s\"\n", ele_path.c_str());
    }

    std::vector<Real> coords;
    std::vector<Int> boundary_markers;
    read_node_file(node_stream, coords, boundary_markers);

    std::vector<LO> conn;
    read_ele_file(ele_stream, conn);

    auto nverts = static_cast<LO>(coords.size() / 2);
    auto ntris = static_cast<LO>(conn.size() / 3);

    HostWrite<Real> host_coords(nverts * 2);
    for (LO i = 0; i < nverts * 2; ++i) {
      host_coords[i] = coords[static_cast<std::size_t>(i)];
    }

    HostWrite<LO> host_conn(ntris * 3);
    for (LO i = 0; i < ntris * 3; ++i) {
      host_conn[i] = conn[static_cast<std::size_t>(i)];
    }

    auto ev2v_dev = LOs(host_conn);
    auto coords_dev = Reals(host_coords);
    build_from_elems_and_coords(&mesh, OMEGA_H_SIMPLEX, FACE, ev2v_dev, coords_dev);

    if (!boundary_markers.empty()) {
      HostWrite<I8> host_bm(nverts);
      for (LO i = 0; i < nverts; ++i) {
        host_bm[i] = static_cast<I8>(
            boundary_markers[static_cast<std::size_t>(i)]);
      }
      mesh.add_tag(VERT, "boundary_marker", 1, Read<I8>(host_bm));
    }
  }
  mesh.set_comm(comm);
  mesh.balance();
  return mesh;
}

void write(filesystem::path const& basename, Mesh* mesh) {
  OMEGA_H_CHECK(mesh->comm()->size() == 1);
  OMEGA_H_CHECK(mesh->dim() == 2);

  auto node_path = basename.string() + ".node";
  auto ele_path = basename.string() + ".ele";

  auto nverts = mesh->nverts();
  auto nelems = mesh->nelems();
  auto dim = mesh->dim();
  auto coords = mesh->coords();
  auto e2v = mesh->ask_elem_verts();

  auto exposed_sides = mark_exposed_sides(mesh);
  auto exposed_nodes = mark_down(mesh, dim - 1, 0, exposed_sides);

  HostRead<Real> h_coords(coords);
  HostRead<LO> h_e2v(e2v);
  HostRead<I8> h_exposed(exposed_nodes);

  std::ofstream node_stream(node_path);
  if (!node_stream.is_open()) {
    Omega_h_fail("couldn't open \"%s\"\n", node_path.c_str());
  }
  node_stream << nverts << " 2 0 1\n";
  for (LO i = 0; i < nverts; ++i) {
    node_stream << i + 1 << " " << h_coords[i * dim] << " "
                << h_coords[i * dim + 1] << " "
                << static_cast<int>(h_exposed[i]) << "\n";
  }
  node_stream.close();

  std::ofstream ele_stream(ele_path);
  if (!ele_stream.is_open()) {
    Omega_h_fail("couldn't open \"%s\"\n", ele_path.c_str());
  }
  ele_stream << nelems << " 3 0\n";
  for (LO i = 0; i < nelems; ++i) {
    ele_stream << i + 1;
    for (LO j = 0; j < dim + 1; ++j) {
      ele_stream << " " << h_e2v[i * (dim + 1) + j] + 1;
    }
    ele_stream << "\n";
  }
  ele_stream.close();
}

}  // namespace xgc

}  // end namespace Omega_h
