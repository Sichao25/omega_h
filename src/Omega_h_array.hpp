#ifndef OMEGA_H_ARRAY_HPP
#define OMEGA_H_ARRAY_HPP

#include <Omega_h_defines.hpp>
#include <Omega_h_fail.hpp>
#include <initializer_list>
#ifdef OMEGA_H_USE_KOKKOS
#include <Omega_h_kokkos.hpp>
#include <Omega_h_pool_kokkos.hpp>
#include <Omega_h_memory.hpp>
#else
#include <Omega_h_shared_alloc.hpp>
#include <memory> //shared_ptr
#include <string>
#endif

namespace Omega_h {

template <typename T>
T* nonnull(T* p);

template <typename T>
class HostWrite;

#ifdef OMEGA_H_USE_KOKKOS
template <typename T>
class KokkosViewWrapper {
 public:
  KokkosViewWrapper(size_t n, const std::string& name_in)
      : view_(KokkosPool::getGlobalPool().allocateView<T>(n)),
        label_(name_in) {}

  [[nodiscard]] auto label() const -> std::string { return label_; }

  [[nodiscard]] auto getView() const -> const View<T*>& {
    return view_;
  }

  ~KokkosViewWrapper() {
    KokkosPool::getGlobalPool().deallocateView<T>(view_);
  }

  View<T*> view_;
  std::string label_;
};

#endif

template <typename T>
class Write {
#ifdef OMEGA_H_USE_KOKKOS
  View<T*> view_;
  SharedRef<KokkosViewWrapper<T>> manager_;  // reference counting
#else
  SharedAlloc shared_alloc_;
#endif

 public:
  using value_type = T;
  OMEGA_H_INLINE Write();
#ifdef OMEGA_H_USE_KOKKOS
  Write(View<T*> view_in);
#endif
  Write(LO size_in, std::string const& name = "");
  Write(LO size_in, T value, std::string const& name = "");
  Write(LO size_in, T offset, T stride, std::string const& name = "");
  Write(std::initializer_list<T> l, std::string const& name = "");
  Write(HostWrite<T> host_write);
  OMEGA_H_INLINE LO size() const OMEGA_H_NOEXCEPT;
  OMEGA_H_DEVICE T& operator[](LO i) const OMEGA_H_NOEXCEPT;
  OMEGA_H_INLINE T* data() const noexcept;
#ifdef OMEGA_H_USE_KOKKOS
  OMEGA_H_INLINE View<T*> const& view() const { return view_; }
#endif
  void set(LO i, T value) const;
  T get(LO i) const;
  long use_count() const;
  OMEGA_H_INLINE bool exists() const noexcept;
#ifdef OMEGA_H_USE_KOKKOS
  std::string name() const;
#else
  std::string const& name() const;
#endif
  OMEGA_H_INLINE T* begin() const noexcept { return data(); }
  OMEGA_H_INLINE T* end() const OMEGA_H_NOEXCEPT { return data() + size(); }
};

template <typename T>
class Read {
  Write<T> write_;

 public:
  using value_type = T;
  OMEGA_H_INLINE Read() {}
  Read(Write<T> write);
  Read(LO size, T value, std::string const& name = "");
  Read(LO size, T offset, T stride, std::string const& name = "");
  Read(std::initializer_list<T> l, std::string const& name = "");
  OMEGA_H_INLINE LO size() const OMEGA_H_NOEXCEPT { return write_.size(); }
  OMEGA_H_DEVICE T const& operator[](LO i) const OMEGA_H_NOEXCEPT {
#ifdef OMEGA_H_CHECK_BOUNDS
    OMEGA_H_CHECK_OP(0, <=, i);
    OMEGA_H_CHECK_OP(i, <, size());
#endif
    return write_[i];
  }
  OMEGA_H_INLINE T const* data() const OMEGA_H_NOEXCEPT {
    return write_.data();
  }
#ifdef OMEGA_H_USE_KOKKOS
  View<const T*> view() const;
#endif
  T get(LO i) const;
  T first() const;
  T last() const;
  OMEGA_H_INLINE bool exists() const OMEGA_H_NOEXCEPT {
    return write_.exists();
  }
  std::string name() const { return write_.name(); }
  OMEGA_H_INLINE T const* begin() const noexcept { return data(); }
  OMEGA_H_INLINE T const* end() const noexcept { return data() + size(); }
};

template <typename T>
Read<T> read(Write<T> a) {
  return Read<T>(a);
}

using Bytes = Read<Byte>;
using LOs = Read<LO>;
using GOs = Read<GO>;
using Reals = Read<Real>;

template <typename T>
class HostRead {
  Read<T> read_;
#if defined(OMEGA_H_USE_KOKKOS)
  typename Kokkos::View<const T*, Kokkos::HostSpace> mirror_;
#elif defined(OMEGA_H_USE_CUDA)
  std::shared_ptr<T[]> mirror_;
#endif
 public:
  using value_type = T;
  HostRead() = default;
  HostRead(Read<T> read);
  LO size() const;
  inline T const& operator[](LO i) const OMEGA_H_NOEXCEPT;
  T const* data() const;
  T get(LO i) const;
  T last() const;
};

template <typename T>
class HostWrite {
  Write<T> write_;
#ifdef OMEGA_H_USE_KOKKOS
  typename View<T*>::host_mirror_type mirror_;
#elif defined(OMEGA_H_USE_CUDA)
  std::shared_ptr<T[]> mirror_;
#endif
 public:
  using value_type = T;
  HostWrite() = default;
  HostWrite(LO size_in, std::string const& name = "");
  
  /**
 * \brief Constructs a HostWrite object with specified size, offset, and stride.
 *
 * \param size_in The number of entries in the array.
 * \param offset The initial value to set for the first entry of the array.
 * \param stride The difference between the values of consecutive entries in the array.
 * \param name The name of the array for identification purposes (default is an empty string).
 *
 * This constructor initializes a HostWrite array with the given size, setting the first entry to the specified offset,
 * and each subsequent entry to the previous entry's value plus the stride. The array is given a name for identification.
 * For example, `HostWrite<Real> h_write(10, 7.0, 0.0);` will create a write array of size 10 and all filled with 7.0.
 */
  HostWrite(LO size_in, T offset, T stride, std::string const& name = "");

  HostWrite(Write<T> write_in);
  HostWrite(std::initializer_list<T> l, std::string const& name = "");
  Write<T> write() const;
  LO size() const OMEGA_H_NOEXCEPT;
  inline T& operator[](LO i) const OMEGA_H_NOEXCEPT;
  T* data() const;
  OMEGA_H_INLINE bool exists() const OMEGA_H_NOEXCEPT {
    return write_.exists();
  }
  void set(LO i, T value);
  T get(LO i) const;
};

template <typename T>
void fill(Write<T> a, T val);
template <typename T>
void fill_linear(Write<T> a, T offset, T stride);
template <class T>
void copy_into(Read<T> a, Write<T> b);
template <class T>
Write<T> deep_copy(Read<T> a, std::string const& name = "");

/* begin explicit instantiation declarations */
#define OMEGA_H_EXPL_INST_DECL(T)                                              \
  extern template T* nonnull(T*);                                              \
  extern template T const* nonnull(T const*);                                  \
  extern template class Read<T>;                                               \
  extern template class Write<T>;                                              \
  extern template class HostRead<T>;                                           \
  extern template class HostWrite<T>;                                          \
  extern template void fill(Write<T> a, T val);                                \
  extern template void fill_linear(Write<T> a, T, T);                          \
  extern template void copy_into(Read<T> a, Write<T> b);                       \
  extern template Write<T> deep_copy(Read<T> a, std::string const&);
OMEGA_H_EXPL_INST_DECL(I8)
OMEGA_H_EXPL_INST_DECL(I32)
OMEGA_H_EXPL_INST_DECL(I64)
OMEGA_H_EXPL_INST_DECL(Real)
#undef OMEGA_H_EXPL_INST_DECL
/* end explicit instantiation declarations */

}  // end namespace Omega_h

#ifdef OMEGA_H_USE_KOKKOS
#include "Omega_h_array_kokkos.hpp"
#else
#include "Omega_h_array_default.hpp"
#endif

#endif
