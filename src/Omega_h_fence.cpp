#include <Omega_h_fence.hpp>
#include <Omega_h_fail.hpp>

#ifdef OMEGA_H_USE_KOKKOS
#include <Omega_h_kokkos.hpp>
#endif

namespace Omega_h {

void fence() {
#if defined(OMEGA_H_USE_KOKKOS)
  Kokkos::fence();
#endif
}

}  // namespace Omega_h
