#ifndef OMEGA_H_EGADS_LITE_HPP
#define OMEGA_H_EGADS_LITE_HPP

#include <Omega_h_array.hpp>
#include <string>

namespace Omega_h {

class Mesh;
struct Egads;

Egads* egads_lite_load(std::string const& filename);
void egads_lite_classify(Egads* eg, int nadj_faces, int const adj_face_ids[],
    int* class_dim, int* class_id);
void egads_lite_free(Egads* eg);
void egads_lite_reclassify(Mesh* mesh, Egads* eg);
Reals egads_lite_get_snap_warp(Mesh* mesh, Egads* eg, bool verbose);

}  // namespace Omega_h

#endif
