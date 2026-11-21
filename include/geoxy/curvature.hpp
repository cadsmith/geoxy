/***************************************************************************************************
 *  @file      geoxy/curvature.hpp
 *  @brief     The unit tangent, the normal and the signed curvature of each curve at a place.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_CURVATURE_HPP
#define GEOXY_CURVATURE_HPP

#include "geoxy/api.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <optional>

namespace Geoxy
{
  /**
   *  The unit tangent of a line.
   *
   *  @param  line  The line.
   *  @param  along Where, 0 at its start and 1 at its end; the same everywhere.
   *  @return       The way from its start to its end, of length 1; none for a line of no length.
   */
  [[nodiscard]] GEOXY_API auto Tangent(const Line &line, double along) -> std::optional<Vector>;

  /**
   *  The unit tangent of a ray.
   *
   *  @param  ray   The ray.
   *  @param  along Where, in its ways from its start; the same everywhere.
   *  @return       Its way, of length 1; none for a way of no length.
   */
  [[nodiscard]] GEOXY_API auto Tangent(const Ray &ray, double along) -> std::optional<Vector>;

  /**
   *  The unit tangent of an xline.
   *
   *  @param  xline The xline.
   *  @param  along Where, in its ways from its point; the same everywhere.
   *  @return       Its way, of length 1; none for a way of no length.
   */
  [[nodiscard]] GEOXY_API auto Tangent(const Xline &xline, double along) -> std::optional<Vector>;

  /**
   *  The unit tangent of an arc, in closed form.
   *
   *  @param  arc   The arc.
   *  @param  along The fraction of its sweep.
   *  @return       The direction of its angle there turned a quarter turn, the way its sweep runs:
   *                each component the double nearest the exact one; none for an arc of no radius or
   *                of no sweep.
   */
  [[nodiscard]] GEOXY_API auto Tangent(const Arc &arc, double along) -> std::optional<Vector>;

  /**
   *  The unit tangent of an ellipse.
   *
   *  @param  ellipse The ellipse.
   *  @param  along   The fraction of its sweep.
   *  @return         Its map's image of the unit circle's tangent, of length 1, the way its sweep
   *                  runs; none for no sweep, or at a tip of an ellipse seen edge on.
   */
  [[nodiscard]] GEOXY_API auto Tangent(const Ellipse &ellipse, double along)
      -> std::optional<Vector>;

  /**
   *  The unit tangent of a polyline.
   *
   *  @param  polyline The polyline.
   *  @param  along    A segment's index and the fraction along it.
   *  @return          The tangent of the segment the value falls in, as Derivative takes it: at a
   *                   vertex the segment leaving it, at the end the last; none for a polyline of no
   *                   segment or a segment of no length.
   */
  [[nodiscard]] GEOXY_API auto Tangent(const Polyline &polyline, double along)
      -> std::optional<Vector>;

  /**
   *  The unit tangent of a spline.
   *
   *  @param  spline The spline.
   *  @param  along  A knot value; at a knot, the span that starts there, as Derivative takes it.
   *  @return        The tangent, of length 1; none where it stops, or where the spline cannot be
   *                 read, as At tells it.
   */
  [[nodiscard]] GEOXY_API auto Tangent(const Spline &spline, double along) -> std::optional<Vector>;

  /**
   *  The unit tangent of a quadratic Bézier curve.
   *
   *  @param  quadratic The curve.
   *  @param  along     0 at its start, 1 at its end.
   *  @return           The tangent, of length 1; none where it stops.
   */
  [[nodiscard]] GEOXY_API auto Tangent(const Bezier::Quadratic &quadratic, double along)
      -> std::optional<Vector>;

  /**
   *  The unit tangent of a cubic Bézier curve.
   *
   *  @param  cubic The curve.
   *  @param  along 0 at its start, 1 at its end.
   *  @return       The tangent, of length 1; none where it stops, as at an end whose control point
   *                it starts on.
   */
  [[nodiscard]] GEOXY_API auto Tangent(const Bezier::Cubic &cubic, double along)
      -> std::optional<Vector>;

  /**
   *  The unit tangent of any curve.
   *
   *  @param  curve The curve.
   *  @param  along The value, as its own kind reads one.
   *  @return       The tangent, as its own kind gives it.
   */
  [[nodiscard]] GEOXY_API auto Tangent(const Curve &curve, double along) -> std::optional<Vector>;

  /**
   *  The unit normal of a line.
   *
   *  @param  line  The line.
   *  @param  along Where.
   *  @return       Its tangent turned a quarter turn counterclockwise, to its left; none where it
   *                has no tangent.
   */
  [[nodiscard]] GEOXY_API auto Normal(const Line &line, double along) -> std::optional<Vector>;

  /**
   *  The unit normal of a ray.
   *
   *  @param  ray   The ray.
   *  @param  along Where.
   *  @return       Its tangent turned a quarter turn counterclockwise; none where it has none.
   */
  [[nodiscard]] GEOXY_API auto Normal(const Ray &ray, double along) -> std::optional<Vector>;

  /**
   *  The unit normal of an xline.
   *
   *  @param  xline The xline.
   *  @param  along Where.
   *  @return       Its tangent turned a quarter turn counterclockwise; none where it has none.
   */
  [[nodiscard]] GEOXY_API auto Normal(const Xline &xline, double along) -> std::optional<Vector>;

  /**
   *  The unit normal of an arc.
   *
   *  @param  arc   The arc.
   *  @param  along The fraction of its sweep.
   *  @return       Its tangent turned a quarter turn counterclockwise: toward its center where it
   *                runs counterclockwise; none where it has no tangent.
   */
  [[nodiscard]] GEOXY_API auto Normal(const Arc &arc, double along) -> std::optional<Vector>;

  /**
   *  The unit normal of an ellipse.
   *
   *  @param  ellipse The ellipse.
   *  @param  along   The fraction of its sweep.
   *  @return         Its tangent turned a quarter turn counterclockwise; none where it has none.
   */
  [[nodiscard]] GEOXY_API auto Normal(const Ellipse &ellipse, double along)
      -> std::optional<Vector>;

  /**
   *  The unit normal of a polyline.
   *
   *  @param  polyline The polyline.
   *  @param  along    A segment's index and the fraction along it.
   *  @return          Its tangent turned a quarter turn counterclockwise; none where it has none.
   */
  [[nodiscard]] GEOXY_API auto Normal(const Polyline &polyline, double along)
      -> std::optional<Vector>;

  /**
   *  The unit normal of a spline.
   *
   *  @param  spline The spline.
   *  @param  along  A knot value.
   *  @return        Its tangent turned a quarter turn counterclockwise; none where it has none.
   */
  [[nodiscard]] GEOXY_API auto Normal(const Spline &spline, double along) -> std::optional<Vector>;

  /**
   *  The unit normal of a quadratic Bézier curve.
   *
   *  @param  quadratic The curve.
   *  @param  along     0 at its start, 1 at its end.
   *  @return           Its tangent turned a quarter turn counterclockwise; none where it has none.
   */
  [[nodiscard]] GEOXY_API auto Normal(const Bezier::Quadratic &quadratic, double along)
      -> std::optional<Vector>;

  /**
   *  The unit normal of a cubic Bézier curve.
   *
   *  @param  cubic The curve.
   *  @param  along 0 at its start, 1 at its end.
   *  @return       Its tangent turned a quarter turn counterclockwise; none where it has none.
   */
  [[nodiscard]] GEOXY_API auto Normal(const Bezier::Cubic &cubic, double along)
      -> std::optional<Vector>;

  /**
   *  The unit normal of any curve.
   *
   *  @param  curve The curve.
   *  @param  along The value, as its own kind reads one.
   *  @return       Its tangent turned a quarter turn counterclockwise; none where it has none.
   */
  [[nodiscard]] GEOXY_API auto Normal(const Curve &curve, double along) -> std::optional<Vector>;

  /**
   *  The signed curvature of a line.
   *
   *  @param  line  The line.
   *  @param  along Where.
   *  @return       0; none for a line of no length, which has no tangent.
   */
  [[nodiscard]] GEOXY_API auto Curvature(const Line &line, double along) -> std::optional<double>;

  /**
   *  The signed curvature of a ray.
   *
   *  @param  ray   The ray.
   *  @param  along Where.
   *  @return       0; none for a way of no length.
   */
  [[nodiscard]] GEOXY_API auto Curvature(const Ray &ray, double along) -> std::optional<double>;

  /**
   *  The signed curvature of an xline.
   *
   *  @param  xline The xline.
   *  @param  along Where.
   *  @return       0; none for a way of no length.
   */
  [[nodiscard]] GEOXY_API auto Curvature(const Xline &xline, double along) -> std::optional<double>;

  /**
   *  The signed curvature of an arc: 1 over its radius, exactly.
   *
   *  @param  arc   The arc.
   *  @param  along The fraction of its sweep; the same everywhere.
   *  @return       The double nearest 1 / |r|, positive where it runs counterclockwise and negative
   *                where it runs clockwise; none for no radius or no sweep.
   */
  [[nodiscard]] GEOXY_API auto Curvature(const Arc &arc, double along) -> std::optional<double>;

  /**
   *  The signed curvature of an ellipse, in closed form.
   *
   *  The unit circle's tangent t mapped by M, the curvature is det M / |M t|^3, signed by its
   *  sweep.
   *
   *  @param  ellipse The ellipse.
   *  @param  along   The fraction of its sweep.
   *  @return         The curvature, positive turning counterclockwise; none where it has no
   *                  tangent.
   */
  [[nodiscard]] GEOXY_API auto Curvature(const Ellipse &ellipse, double along)
      -> std::optional<double>;

  /**
   *  The signed curvature of a polyline: of the segment the value falls in.
   *
   *  @param  polyline The polyline.
   *  @param  along    A segment's index and the fraction along it.
   *  @return          0 on a line; 4 b / (c (1 + b^2)) on a bulge b of chord c, 1 over its radius
   *                   signed by its bulge; none where it has no tangent.
   */
  [[nodiscard]] GEOXY_API auto Curvature(const Polyline &polyline, double along)
      -> std::optional<double>;

  /**
   *  The signed curvature of a spline.
   *
   *  From the first and second derivatives of its span's rational Bézier piece in homogeneous
   *  coordinates, by the quotient rule: x' × x'' / |x'|^3.
   *
   *  @param  spline The spline.
   *  @param  along  A knot value; at a knot, the span that starts there.
   *  @return        The curvature, positive turning counterclockwise; none where it has no tangent.
   */
  [[nodiscard]] GEOXY_API auto Curvature(const Spline &spline, double along)
      -> std::optional<double>;

  /**
   *  The signed curvature of a quadratic Bézier curve.
   *
   *  @param  quadratic The curve.
   *  @param  along     0 at its start, 1 at its end.
   *  @return           x' × x'' / |x'|^3, positive turning counterclockwise; none where it stops.
   */
  [[nodiscard]] GEOXY_API auto Curvature(const Bezier::Quadratic &quadratic, double along)
      -> std::optional<double>;

  /**
   *  The signed curvature of a cubic Bézier curve.
   *
   *  @param  cubic The curve.
   *  @param  along 0 at its start, 1 at its end.
   *  @return       x' × x'' / |x'|^3, positive turning counterclockwise; none where it stops.
   */
  [[nodiscard]] GEOXY_API auto Curvature(const Bezier::Cubic &cubic, double along)
      -> std::optional<double>;

  /**
   *  The signed curvature of any curve.
   *
   *  @param  curve The curve.
   *  @param  along The value, as its own kind reads one.
   *  @return       The curvature, as its own kind gives it.
   */
  [[nodiscard]] GEOXY_API auto Curvature(const Curve &curve, double along) -> std::optional<double>;
} // namespace Geoxy

#endif /* GEOXY_CURVATURE_HPP */
