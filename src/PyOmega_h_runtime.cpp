#include <PyOmega_h.hpp>

#include <stdexcept>

#ifdef OMEGA_H_USE_MPI
#include <Omega_h_mpi.h>
#endif
#ifdef OMEGA_H_USE_KOKKOS
#include <Omega_h_pool_kokkos.hpp>
#endif

namespace Omega_h {

namespace {
#ifdef OMEGA_H_USE_MPI
bool g_py_owns_mpi = false;
#endif
#ifdef OMEGA_H_USE_KOKKOS
bool g_py_owns_kokkos = false;
#endif
}  // namespace


void initialize() {
#ifdef OMEGA_H_USE_MPI
  int is_initialized = 0;
  MPI_Initialized(&is_initialized);
  if (!is_initialized) {
    if (MPI_Init(nullptr, nullptr) != MPI_SUCCESS) {
      throw std::runtime_error("PyOmega_h: MPI_Init failed");
    }
    g_py_owns_mpi = true;
  }
#endif
#ifdef OMEGA_H_USE_KOKKOS
  if (!Kokkos::is_initialized()) {
    Kokkos::initialize();
    g_py_owns_kokkos = true;
  }
#endif
}

void finalize() {
#ifdef OMEGA_H_USE_KOKKOS
  if (g_py_owns_kokkos && !Kokkos::is_finalized()) {
    KokkosPool::destroyGlobalPool();
    Kokkos::finalize();
    g_py_owns_kokkos = false;
  }
#endif
#ifdef OMEGA_H_USE_MPI
  int is_finalized = 0;
  MPI_Finalized(&is_finalized);
  if (g_py_owns_mpi && !is_finalized) {
    MPI_Finalize();
    g_py_owns_mpi = false;
  }
#endif
}

}  // namespace Omega_h
