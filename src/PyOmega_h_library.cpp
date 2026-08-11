#include <Omega_h_library.hpp>
#include <PyOmega_h.hpp>
#include <iostream>

namespace Omega_h {

/* The lifetime of the Library object is quite important (it must contain the
   lifetime of pretty much all other Omega_h objects). We use a shared_ptr
   holder so that py::keep_alive policies on Mesh and Comm objects can keep
   the Library alive until all dependent objects have been destroyed.
   This avoids the "Attempting to use an MPI routine after finalizing MPICH"
   error that occurs when Python GC destroys the Library before Meshes/Comms.
 */

void pybind11_library(py::module& m) {
  // Bind Omega_h::Library with shared_ptr holder so that keep_alive
  // policies on dependent objects (Mesh, Comm) can extend its lifetime.
  py::class_<Omega_h::Library, std::shared_ptr<Omega_h::Library>>(
    m, "OmegaHLibrary")
    .def(py::init([]() { return std::make_shared<Omega_h::Library>(); }),
         "Default constructor")
    .def("world",
         [](std::shared_ptr<Omega_h::Library> self) { return self->world(); },
         py::keep_alive<0, 1>(),
         "Get the world communicator (keeps library alive)");
}

}  // namespace Omega_h
