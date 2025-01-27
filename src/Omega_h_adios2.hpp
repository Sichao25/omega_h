#ifndef OMEGA_H_ADIOS2_HPP
#define OMEGA_H_ADIOS2_HPP
#include <adios2.h>
#include <Omega_h_mesh.hpp>
//#if ADIOS2_USE_MPI
//#include <mpi.h>
//#endif
#include "Omega_h_filesystem.hpp" // filesystem
#include "Omega_h_library.hpp"

namespace Omega_h {
void write_adios2(filesystem::path const& path,
                  std::map<Mesh*, std::string>& mesh_map);
void write_adios2(filesystem::path const& path, Mesh *mesh, std::string pref="");

Mesh read_adios2(filesystem::path const& path, Library* lib, std::string pref="");
}  // namespace Omega_h
#endif
