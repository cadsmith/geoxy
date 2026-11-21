/***************************************************************************************************
 *  @file      geoxy/chamfer.hpp
 *  @brief     CHAMFER: two curves cut back from their corner by two distances, or by a distance and
 *             an angle, and joined by a line; and every corner of a polyline cut so.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_CHAMFER_HPP
#define GEOXY_CHAMFER_HPP

#include "geoxy/angle.hpp"
#include "geoxy/api.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/fillet.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/tolerance.hpp"

#include <array>
#include <expected>
#include <optional>

namespace Geoxy
{
  /**
   *  How far a chamfer cuts back from a corner: CHAMFER's two distances, or its distance and angle.
   */
  struct Bevel
  {
    double               First  = 0.0; /**< How far from the corner along the first curve.    */
    double               Second = 0.0; /**< How far along the second; not read with an angle. */
    std::optional<Angle> Slant  = {};  /**< CHAMFER's angle: the chamfer's from the first curve,
                                            toward the corner, in (0, 180) degrees; where given,
                                            the chamfer runs from the first cut at it until it
                                            meets the second.                                   */
  };

  /**
   *  What CHAMFER asks of two curves: how far it cuts back, and where each was picked.
   */
  struct Chamfering
  {
    Bevel                Cut   = {}; /**< How far it cuts back. */
    std::array<Point, 2> Picks = {}; /**< Where each curve was picked, on the part of it to
                                          keep.                                                  */
  };

  /**
   *  Two curves chamfered: each cut back from their corner, and the line joining the cuts.
   */
  struct Bevelled
  {
    Curve               First  = {}; /**< The first curve, cut back. */
    Curve               Second = {}; /**< The second, likewise.      */
    std::optional<Line> Edge   = {}; /**< The chamfer, from the first's cut to the second's;
                                          none where both distances are 0.                     */
  };

  /**
   *  Two curves cut back from their corner and joined by a line, as CHAMFER draws them.
   *
   *  The corner and the parts kept are FILLET's: where the curves, taken as what they lie on, meet
   *  nearest the picks, each keeping the side its pick lies on. Each is cut where it is the
   *  distance from the corner, measured straight: along a line, and across to an arc, an ellipse, a
   *  spline or a Bézier curve; with an angle, the second is cut where the line from the first's
   *  cut, at that angle from the first toward the corner, meets it. Between two lines the cuts are
   *  the doubles nearest the exact points, the angle's sine and cosine taken as the doubles
   *  Direction gives them; on curves, where Intersect rounds the crossing of a circle or a ray
   *  about the corner as found.
   *
   *  @param  first      One curve.
   *  @param  second     The other.
   *  @param  chamfering How far it cuts back, and the picks.
   *  @return            The curves and the chamfer; Parallel for two parallel lines, Apart for
   *                     curves that never meet, TooLarge where a cut falls past the far end of the
   *                     part kept or the angle's line never meets the second's part, Degenerate for
   *                     a distance below nothing or an angle outside (0, 180) degrees, NotFinite
   *                     for a number that is not finite.
   */
  [[nodiscard]] GEOXY_API auto Chamfer(const Curve &first, const Curve &second,
                                       const Chamfering &chamfering)
      -> std::expected<Bevelled, Misfit>;

  /**
   *  A polyline with each corner cut back and joined by a line where it fits, as CHAMFER's Polyline
   *  option cuts one.
   *
   *  Each corner where two segments meet at an angle is cut as two curves are chamfered, the
   *  segment before it taking the first distance and the one after the second; the corners are cut
   *  in order, each within what the one before left of their segment. A corner whose cut does not
   *  fit, past what is left of a segment by more than the tolerance, is left, and reported; a
   *  segment left shorter than the tolerance is dropped; a corner where the two segments run on
   *  smoothly is passed by.
   *
   *  @param  polyline  The polyline.
   *  @param  bevel     How far each corner is cut back, a distance above nothing.
   *  @param  tolerance Within which a cut still fits and a segment is of no length.
   *  @return           The polyline and the corners left; NotFinite for a number that is not
   *                    finite, Degenerate for no distance above nothing, a distance below nothing
   *                    or an angle outside (0, 180) degrees.
   */
  [[nodiscard]] GEOXY_API auto Chamfer(const Polyline &polyline, const Bevel &bevel,
                                       const Tolerance &tolerance)
      -> std::expected<Reshaped, Misfit>;
} // namespace Geoxy

#endif /* GEOXY_CHAMFER_HPP */
