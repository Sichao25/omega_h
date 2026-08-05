#ifndef OMEGA_H_MIXEDMESH_HPP
#define OMEGA_H_MIXEDMESH_HPP

#include <Omega_h_adj.hpp>
#include <Omega_h_comm.hpp>
#include <Omega_h_dist.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_tag.hpp>
#include <array>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace Omega_h {

class MixedMesh {
 public:
  MixedMesh();
  MixedMesh(Library* library);
  void set_library(Library* library);
  void set_comm(CommPtr const& new_comm);
  void set_dim(Int dim_in);
  void set_verts(LO nverts_in);
  void set_ents(Int ent_dim, Adj down);
  Library* library() const;
  CommPtr comm() const;
  inline Int dim() const {
    OMEGA_H_CHECK(0 <= dim_ && dim_ <= 3);
    return dim_;
  }
  inline Omega_h_Family family() const { return family_; }

  void set_verts_type(LO nverts_in);
  void set_ents(Topo_type high_type, Topo_type low_type, Adj h2l);
  LO nents(Topo_type ent_type) const;
  Int ent_dim(Topo_type ent_type) const;

  LO npyrams() const;
  LO nwedges() const;
  LO nhexs() const;
  LO ntets() const;
  LO nquads() const;
  LO ntris() const;
  LO nregions_mix() const;
  LO nfaces_mix() const;
  LO nedges_mix() const;
  LO nverts_mix() const;

  template <typename T>
  void add_tag(Topo_type ent_type, std::string const& name, Int ncomps);
  template <typename T>
  void add_tag(Topo_type ent_type, std::string const& name, Int ncomps,
      Read<T> array, bool internal = false);
  template <typename T>
  void set_tag(Topo_type ent_type, std::string const& name, Read<T> array,
      bool internal = false);
  TagBase const* get_tagbase(Topo_type ent_type, std::string const& name) const;
  template <typename T>
  Tag<T> const* get_tag(Topo_type ent_type, std::string const& name) const;
  template <typename T>
  Read<T> get_array(Topo_type ent_type, std::string const& name) const;
  void remove_tag(Topo_type ent_type, std::string const& name);
  bool has_tag(Topo_type ent_type, std::string const& name) const;
  [[nodiscard]] Int ntags(Topo_type ent_type) const;
  TagBase const* get_tag(Topo_type ent_type, Int i) const;
  bool has_ents(Topo_type ent_type) const;
  bool has_adj(Topo_type from_type, Topo_type to_type) const;
  Adj get_adj(Topo_type from_type, Topo_type to_type) const;
  Adj ask_down(Topo_type from_type, Topo_type to_type);
  LOs ask_verts_of(Topo_type ent_type);
  Adj ask_up(Topo_type from_type, Topo_type to_type);

  typedef std::shared_ptr<const TagBase> TagPtr;
  typedef std::shared_ptr<const Adj> AdjPtr;

 private:
  typedef std::vector<TagPtr> TagVector;
  typedef TagVector::iterator TagIter;
  typedef TagVector::const_iterator TagCIter;
  TagIter tag_iter(Topo_type ent_type, std::string const& name);
  TagCIter tag_iter(Topo_type ent_type, std::string const& name) const;
  void check_type(Topo_type ent_type) const;
  void check_type2(Topo_type ent_type) const;
  void add_adj(Topo_type from_type, Topo_type to_type, Adj adj);
  Adj derive_adj(Topo_type from_type, Topo_type to_type);
  Adj ask_adj(Topo_type from_type, Topo_type to_type);
  void react_to_set_tag(Topo_type ent_type, std::string const& name);
  Library* library_;
  Omega_h_Family family_;
  Int dim_;
  CommPtr comm_;
  LO nents_[DIMS];
  TagVector tags_[DIMS];
  LO nents_type_[TOPO_TYPES];
  TagVector tags_type_[TOPO_TYPES];
  AdjPtr adjs_type_[TOPO_TYPES][TOPO_TYPES];

 public:
  void add_coords_mix(Reals array);
  Reals coords_mix() const;
};

using TagSet = std::array<std::set<std::string>, DIMS>;
  
void get_all_type_tags(MixedMesh* mesh, Int dim, Topo_type ent_type, TagSet* tags);

#define OMEGA_H_EXPL_INST_DECL(T)                                              \
  extern template Read<T> MixedMesh::get_array<T>(                                  \
      Topo_type ent_type, std::string const& name) const;                      \
  extern template void MixedMesh::add_tag<T>(                                       \
      Topo_type ent_type, std::string const& name, Int ncomps);                \
  extern template void MixedMesh::add_tag<T>(Topo_type ent_type,                    \
      std::string const& name, Int ncomps, Read<T> array, bool internal);      \
  extern template void MixedMesh::set_tag(Topo_type ent_type,                       \
      std::string const& name, Read<T> array, bool internal);
OMEGA_H_EXPL_INST_DECL(I8)
OMEGA_H_EXPL_INST_DECL(I32)
OMEGA_H_EXPL_INST_DECL(I64)
OMEGA_H_EXPL_INST_DECL(Real)
#undef OMEGA_H_EXPL_INST_DECL

}  // namespace Omega_h

#endif
