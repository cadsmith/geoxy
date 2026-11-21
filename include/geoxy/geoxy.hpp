/***************************************************************************************************
 *  @file      geoxy/geoxy.hpp
 *  @brief     Every public header of geoxy, in one include.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_GEOXY_HPP
#define GEOXY_GEOXY_HPP

#include "geoxy/accuracy.hpp"         // IWYU pragma: export
#include "geoxy/angle.hpp"            // IWYU pragma: export
#include "geoxy/arc.hpp"              // IWYU pragma: export
#include "geoxy/area.hpp"             // IWYU pragma: export
#include "geoxy/arrangement.hpp"      // IWYU pragma: export
#include "geoxy/at.hpp"               // IWYU pragma: export
#include "geoxy/bezier/cubic.hpp"     // IWYU pragma: export
#include "geoxy/bezier/quadratic.hpp" // IWYU pragma: export
#include "geoxy/box.hpp"              // IWYU pragma: export
#include "geoxy/chamfer.hpp"          // IWYU pragma: export
#include "geoxy/clip.hpp"             // IWYU pragma: export
#include "geoxy/cloud.hpp"            // IWYU pragma: export
#include "geoxy/construct.hpp"        // IWYU pragma: export
#include "geoxy/contain.hpp"          // IWYU pragma: export
#include "geoxy/curvature.hpp"        // IWYU pragma: export
#include "geoxy/curve.hpp"            // IWYU pragma: export
#include "geoxy/dash.hpp"             // IWYU pragma: export
#include "geoxy/derivative.hpp"       // IWYU pragma: export
#include "geoxy/divide.hpp"           // IWYU pragma: export
#include "geoxy/domain.hpp"           // IWYU pragma: export
#include "geoxy/ellipse.hpp"          // IWYU pragma: export
#include "geoxy/ends.hpp"             // IWYU pragma: export
#include "geoxy/error.hpp"            // IWYU pragma: export
#include "geoxy/extents.hpp"          // IWYU pragma: export
#include "geoxy/fillet.hpp"           // IWYU pragma: export
#include "geoxy/fit.hpp"              // IWYU pragma: export
#include "geoxy/flatten.hpp"          // IWYU pragma: export
#include "geoxy/grid.hpp"             // IWYU pragma: export
#include "geoxy/inertia.hpp"          // IWYU pragma: export
#include "geoxy/intersect.hpp"        // IWYU pragma: export
#include "geoxy/join.hpp"             // IWYU pragma: export
#include "geoxy/length.hpp"           // IWYU pragma: export
#include "geoxy/line.hpp"             // IWYU pragma: export
#include "geoxy/loop.hpp"             // IWYU pragma: export
#include "geoxy/matrix.hpp"           // IWYU pragma: export
#include "geoxy/nearest.hpp"          // IWYU pragma: export
#include "geoxy/offset.hpp"           // IWYU pragma: export
#include "geoxy/orient.hpp"           // IWYU pragma: export
#include "geoxy/outline.hpp"          // IWYU pragma: export
#include "geoxy/part.hpp"             // IWYU pragma: export
#include "geoxy/place.hpp"            // IWYU pragma: export
#include "geoxy/point.hpp"            // IWYU pragma: export
#include "geoxy/polyline.hpp"         // IWYU pragma: export
#include "geoxy/prepared.hpp"         // IWYU pragma: export
#include "geoxy/ray.hpp"              // IWYU pragma: export
#include "geoxy/region.hpp"           // IWYU pragma: export
#include "geoxy/reverse.hpp"          // IWYU pragma: export
#include "geoxy/segment.hpp"          // IWYU pragma: export
#include "geoxy/space/frame.hpp"      // IWYU pragma: export
#include "geoxy/space/vector.hpp"     // IWYU pragma: export
#include "geoxy/spline.hpp"           // IWYU pragma: export
#include "geoxy/tangents.hpp"         // IWYU pragma: export
#include "geoxy/tolerance.hpp"        // IWYU pragma: export
#include "geoxy/transform.hpp"        // IWYU pragma: export
#include "geoxy/tree.hpp"             // IWYU pragma: export
#include "geoxy/triangulate.hpp"      // IWYU pragma: export
#include "geoxy/vector.hpp"           // IWYU pragma: export
#include "geoxy/version.hpp"          // IWYU pragma: export
#include "geoxy/workers.hpp"          // IWYU pragma: export
#include "geoxy/xline.hpp"            // IWYU pragma: export

/**
 *  Points, vectors, maps, angles, boxes and tolerances, and the curves a drawing is made of with
 *  what is worked out from them. Nothing in it names an entity, reads a file or holds state: every
 *  function is a value in, a value out.
 *
 *  @namespace Geoxy The geometry and math kernel of a 2D CAD application.
 */

/**
 *  @namespace Geoxy::Space
 *  What leaves the plane: vectors in space, and the frames of planes by their normals.
 */

#endif /* GEOXY_GEOXY_HPP */
