/***************************************************************************************************
 *  @file      geoxy/area.hpp
 *  @brief     The area a curve encloses and its centroid, an open curve closed by the line from its
 *             end to its start, as AutoCAD's AREA closes it.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_AREA_HPP
#define GEOXY_AREA_HPP

#include "geoxy/accuracy.hpp"
#include "geoxy/api.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/xline.hpp"

#include <optional>

namespace Geoxy
{
  /**
   *  The area a line encloses.
   *
   *  @param  line The line.
   *  @return      0.
   */
  [[nodiscard]] GEOXY_API auto Area(const Line &line) -> double;

  /**
   *  The area a ray encloses.
   *
   *  @param  ray The ray.
   *  @return     0: it never ends, and closes nothing.
   */
  [[nodiscard]] GEOXY_API auto Area(const Ray &ray) -> double;

  /**
   *  The area an xline encloses.
   *
   *  @param  xline The xline.
   *  @return       0: it never ends, and closes nothing.
   */
  [[nodiscard]] GEOXY_API auto Area(const Xline &xline) -> double;

  /**
   *  The area an arc encloses with its chord, in closed form: r^2 (s - sin s) / 2 for the radius r
   *  and the sweep s, whatever the arc's place, s - sin s taken from its series where the two
   *  nearly cancel. Within a step of a double of the area its sweep has.
   *
   *  @param  arc The arc.
   *  @return     The signed area, positive for an arc that runs counterclockwise; every turn of it
   *              counted.
   */
  [[nodiscard]] GEOXY_API auto Area(const Arc &arc) -> double;

  /**
   *  The area an ellipse encloses with its chord, in closed form: det M (s - sin s) / 2 for its map
   *  M and its sweep s, whatever its place, as an arc's.
   *
   *  @param  ellipse The ellipse.
   *  @return         The signed area, positive for an ellipse that runs counterclockwise; every
   *                  turn of it counted.
   */
  [[nodiscard]] GEOXY_API auto Area(const Ellipse &ellipse) -> double;

  /**
   *  The area a polyline encloses, closed by the line from its last vertex to its first when it is
   *  open.
   *
   *  The straight segments' share is summed by the shoelace formula about the first vertex, each
   *  product and each sum carried with its rounding error, as Ogita, Rump and Oishi do: as accurate
   *  as working in twice the precision. A bulged segment adds what its arc encloses with its chord,
   *  worked out from the chord and the bulge alone: c^2 (b / 3 + b^3 / 15 - ...) for a bulge b up
   *  to 1/2 and a chord of length c, and r^2 (s - sin s) / 2 past it, the radius and the sweep from
   *  the bulge.
   *
   *  @param  polyline The polyline.
   *  @return          The signed area, positive for a polyline that runs counterclockwise; a part
   *                   that runs the other way subtracts; 0 for fewer than two vertices.
   */
  [[nodiscard]] GEOXY_API auto Area(const Polyline &polyline) -> double;

  /**
   *  The area a spline encloses with its chord: for a spline that is not rational, of degree 2 to
   *  26, the double nearest it, its Bézier pieces' control points and their closed form worked out
   *  in double-double with every rounding bounded, and where that bound does not settle it, as for
   *  one of a higher degree, the sum of its pieces' Bernstein forms; for one of degree 1, the
   *  polygon of its poles; for a rational one, strip by strip as Length measures it, X x X' / W^2
   *  integrated by Gauss-Legendre rules whose error is bounded from the weight's roots: within a
   *  few steps of a double of the area, or of the integral of its rate's length where its parts
   *  cancel. A spline and its reverse are measured along the same way, whichever of the two reads
   *  first, so that the reverse's area is the area negated exactly; one that reads the same either
   *  way runs out and back along itself and encloses nothing.
   *
   *  @param  spline   The spline.
   *  @param  accuracy How far the area may be off, for a rational spline.
   *  @return          The signed area, positive for a spline that runs counterclockwise; 0 for one
   *                   that is not valid.
   */
  [[nodiscard]] GEOXY_API auto Area(const Spline &spline, Accuracy accuracy) -> double;

  /**
   *  The area a quadratic Bézier curve encloses with its chord, in closed form: p1 x p2 / 3, p_i
   *  its control points less its start, each held exactly in two doubles and the product in
   *  double-double, so that it is the double nearest the area but where the product cancels far
   *  below its terms.
   *
   *  @param  quadratic The curve.
   *  @return           The signed area, positive for a curve that runs counterclockwise.
   */
  [[nodiscard]] GEOXY_API auto Area(const Bezier::Quadratic &quadratic) -> double;

  /**
   *  The area a cubic Bézier curve encloses with its chord, in closed form: 3 (p1 x p2 + p1 x p3 +
   *  2 p2 x p3) / 20, p_i its control points less its start, each held exactly in two doubles and
   *  the sum in double-double, so that it is the double nearest the area but where the products
   *  cancel far below themselves.
   *
   *  @param  cubic The curve.
   *  @return       The signed area, positive for a curve that runs counterclockwise.
   */
  [[nodiscard]] GEOXY_API auto Area(const Bezier::Cubic &cubic) -> double;

  /**
   *  The area any curve encloses.
   *
   *  @param  curve    The curve.
   *  @param  accuracy How far the area may be off, for a kind that is approximated.
   *  @return          The signed area, as its own kind gives it.
   */
  [[nodiscard]] GEOXY_API auto Area(const Curve &curve, Accuracy accuracy) -> double;

  /**
   *  The centroid of the area a line encloses.
   *
   *  @param  line The line.
   *  @return      None: it encloses nothing.
   */
  [[nodiscard]] GEOXY_API auto Centroid(const Line &line) -> std::optional<Point>;

  /**
   *  The centroid of the area a ray encloses.
   *
   *  @param  ray The ray.
   *  @return     None: it encloses nothing.
   */
  [[nodiscard]] GEOXY_API auto Centroid(const Ray &ray) -> std::optional<Point>;

  /**
   *  The centroid of the area an xline encloses.
   *
   *  @param  xline The xline.
   *  @return       None: it encloses nothing.
   */
  [[nodiscard]] GEOXY_API auto Centroid(const Xline &xline) -> std::optional<Point>;

  /**
   *  The centroid of the area an arc encloses with its chord, in closed form: 4 r sin^3(s / 2) / (3
   *  (s - sin s)) from the center along the bisector of the sweep.
   *
   *  @param  arc The arc.
   *  @return     The centroid, the moments' over the signed area; none for no area.
   */
  [[nodiscard]] GEOXY_API auto Centroid(const Arc &arc) -> std::optional<Point>;

  /**
   *  The centroid of the area an ellipse encloses with its chord, in closed form: the arc's of the
   *  unit circle, mapped.
   *
   *  @param  ellipse The ellipse.
   *  @return         The centroid, the moments' over the signed area; none for no area.
   */
  [[nodiscard]] GEOXY_API auto Centroid(const Ellipse &ellipse) -> std::optional<Point>;

  /**
   *  The centroid of the area a polyline encloses, closed when it is open.
   *
   *  @param  polyline The polyline.
   *  @return          The centroid, the moments' over the signed area, parts that run the other way
   *                   subtracting; none for no area.
   */
  [[nodiscard]] GEOXY_API auto Centroid(const Polyline &polyline) -> std::optional<Point>;

  /**
   *  The centroid of the area a spline encloses with its chord, its moments worked out as Area
   *  works out the area, X / W times its rate: the same, to the bit, for a spline and its reverse.
   *
   *  @param  spline   The spline.
   *  @param  accuracy How far the area and each moment may be off, for a rational spline.
   *  @return          The centroid, the moments' over the signed area; none for no area, or for a
   *                   spline that is not valid.
   */
  [[nodiscard]] GEOXY_API auto Centroid(const Spline &spline, Accuracy accuracy)
      -> std::optional<Point>;

  /**
   *  The centroid of the area a quadratic Bézier curve encloses with its chord.
   *
   *  @param  quadratic The curve.
   *  @return           The centroid; none for no area.
   */
  [[nodiscard]] GEOXY_API auto Centroid(const Bezier::Quadratic &quadratic) -> std::optional<Point>;

  /**
   *  The centroid of the area a cubic Bézier curve encloses with its chord.
   *
   *  @param  cubic The curve.
   *  @return       The centroid; none for no area.
   */
  [[nodiscard]] GEOXY_API auto Centroid(const Bezier::Cubic &cubic) -> std::optional<Point>;

  /**
   *  The centroid of the area any curve encloses.
   *
   *  @param  curve    The curve.
   *  @param  accuracy How far the area and each moment may be off, for a kind that is approximated.
   *  @return          The centroid, as its own kind gives it.
   */
  [[nodiscard]] GEOXY_API auto Centroid(const Curve &curve, Accuracy accuracy)
      -> std::optional<Point>;
} // namespace Geoxy

#endif /* GEOXY_AREA_HPP */
