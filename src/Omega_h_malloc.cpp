#include <Omega_h_fail.hpp>
#include <Omega_h_malloc.hpp>
#include <Omega_h_pool.hpp>
#include <Omega_h_profile.hpp>
#include <cstdlib>

namespace Omega_h {

void* host_malloc(std::size_t size) {
  OMEGA_H_TIME_FUNCTION;
  return ::std::malloc(size);
}

void host_free(void* ptr, std::size_t) {
  OMEGA_H_TIME_FUNCTION;
  ::std::free(ptr);
}

static Pool* host_pool = nullptr;

//if true then either host pooling for a non-kokkos build is enabled or
//kokkos based pooling is enabled for the default memory space
static bool pooling_enabled = false;

void enable_pooling() {
  host_pool = new Pool(host_malloc, host_free);
  pooling_enabled = true;
}

void disable_pooling() {
  delete host_pool;
  host_pool = nullptr;
  pooling_enabled = false;
}

bool is_pooling_enabled() { return pooling_enabled; }

void* maybe_pooled_host_malloc(std::size_t size) {
  #ifdef OMEGA_H_USE_KOKKOS
  Omega_h::fail("This function should not be called when Kokkos is enabled. "
      "Use SharedRef<...> if reference counting is needed or directly call "
      "KokkosPool::getGlobalPool().allocate(...) or Kokkos::kokkos_malloc(...)\n");
  #else
  if (host_pool) return allocate(*host_pool, size);
  return host_malloc(size);
  #endif
}

void maybe_pooled_host_free(void* ptr, std::size_t size) {
  #ifdef OMEGA_H_USE_KOKKOS
  Omega_h::fail("This function should not be called when Kokkos is enabled. "
      "Use SharedRef<...> if reference counting is needed or directly call "
      "KokkosPool::getGlobalPool().deallocate(...) or Kokkos::kokkos_free(...)\n");
  #else
  if (host_pool) deallocate(*host_pool, ptr, size);
  else host_free(ptr, size);
  #endif
}
}  // namespace Omega_h
