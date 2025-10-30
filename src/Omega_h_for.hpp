#ifndef OMEGA_H_FOR_HPP
#define OMEGA_H_FOR_HPP

#include <Omega_h_defines.hpp>
#include <Omega_h_profile.hpp>

#include <Omega_h_int_iterator.hpp>
#include <Omega_h_shared_alloc.hpp>

#ifdef OMEGA_H_USE_KOKKOS
#include <Omega_h_kokkos.hpp>
#endif

namespace Omega_h {

#if defined(OMEGA_H_USE_KOKKOS) || !defined(OMEGA_H_USE_CUDA)

template <typename InputIterator, typename UnaryFunction>
void for_each(InputIterator first, InputIterator last, UnaryFunction&& f) {
  if (first >= last) return;
  Omega_h::entering_parallel = true;
  auto const f2 = std::move(f);
  Omega_h::entering_parallel = false;
#if defined(OMEGA_H_USE_KOKKOS)
  LO const n = last - first;
  if (n > 0) Kokkos::parallel_for(policy(n), f2);
#elif defined(OMEGA_H_USE_OPENMP)
  LO const n = last - first;
#pragma omp parallel for
  for (LO i = 0; i < n; ++i) {
    f2(first[i]);
  }
#else
  for (; first != last; ++first) {
    f2(*first);
  }
#endif
}

#endif
template <typename UnaryFunction>
void parallel_for(LO n, UnaryFunction&& f) {
  OMEGA_H_TIME_FUNCTION;
  auto const first = IntIterator(0);
  auto const last = IntIterator(n);
  ::Omega_h::for_each(first, last, f);
}

template <typename T>
void parallel_for(LO n, T const& f, char const* name = "") {
#if defined(OMEGA_H_USE_KOKKOS)
  if (n > 0) {
    Kokkos::parallel_for(name, policy(n), f);
    Kokkos::fence();
  }
#else
  (void)name;
  auto const first = IntIterator(0);
  auto const last = IntIterator(n);
  auto f2 = f;
  ::Omega_h::for_each(first, last, std::move(f2));
#endif
}

template <typename T>
void parallel_for(char const* name, LO n, T const& f) {
  parallel_for(n, f, name);
}

}  // end namespace Omega_h

#endif
