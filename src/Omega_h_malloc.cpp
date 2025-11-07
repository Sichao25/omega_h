#include <Omega_h_fail.hpp>
#include <Omega_h_malloc.hpp>
#include <Omega_h_pool.hpp>
#ifdef OMEGA_H_USE_KOKKOS
#include <Omega_h_pool_kokkos.hpp>
#endif
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

void* maybe_pooled_device_malloc(std::size_t size) {
  #ifdef OMEGA_H_USE_KOKKOS
  if (pooling_enabled) return KokkosPool::getGlobalPool().allocate(size);
  return Kokkos::kokkos_malloc(size);
  #else
  Omega_h_fail("Non-Kokkos pooled device malloc is not implemented\n");
  #endif
}

void maybe_pooled_device_free(void* ptr, std::size_t size) {
  #ifdef OMEGA_H_USE_KOKKOS
  if (pooling_enabled)
    KokkosPool::getGlobalPool().deallocate(ptr);
  return Kokkos::kokkos_free(ptr);
  #else
  Omega_h_fail("Non-Kokkos pooled device free is not implemented\n");
  #endif
}

void* maybe_pooled_host_malloc(std::size_t size) {
  #ifdef OMEGA_H_USE_KOKKOS
  if (pooling_enabled) return KokkosPool::getGlobalPool().allocate(size);
  return Kokkos::kokkos_malloc(size);
  #else
  if (pooling_enabled) return allocate(*host_pool, size);
  return host_malloc(size);
  #endif
}

void maybe_pooled_host_free(void* ptr, std::size_t size) {
  #ifdef OMEGA_H_USE_KOKKOS
  if (pooling_enabled)
    KokkosPool::getGlobalPool().deallocate(ptr);
  return Kokkos::kokkos_free(ptr);
  #else
  if (pooling_enabled) deallocate(*host_pool, ptr, size);
  else
    host_free(ptr, size);
  #endif
}
}  // namespace Omega_h
