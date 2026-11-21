/***************************************************************************************************
 *  @file      geoxy/clip.hpp
 *  @brief     Curves cut to a box: the stretches of a curve's parameter that lie within it.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_CLIP_HPP
#define GEOXY_CLIP_HPP

#include "geoxy/api.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/xline.hpp"

#include <vector>

namespace Geoxy
{
  /**
   *  The stretches of a line within a box, its edges included.
   *
   *  @param  line The line.
   *  @param  box  The box.
   *  @return      At most one stretch of its parameter, longer than nothing; none for an empty box.
   */
  [[nodiscard]] GEOXY_API auto Clip(const Line &line, const Box &box) -> std::vector<Range>;

  /**
   *  The stretches of a ray within a box, its edges included.
   *
   *  @param  ray The ray.
   *  @param  box The box.
   *  @return     At most one stretch of its parameter, longer than nothing; none for an empty box.
   */
  [[nodiscard]] GEOXY_API auto Clip(const Ray &ray, const Box &box) -> std::vector<Range>;

  /**
   *  The stretches of an xline within a box, its edges included.
   *
   *  @param  xline The xline.
   *  @param  box   The box.
   *  @return       At most one stretch of its parameter, longer than nothing; none for an empty
   *                box.
   */
  [[nodiscard]] GEOXY_API auto Clip(const Xline &xline, const Box &box) -> std::vector<Range>;

  /**
   *  The stretches of an arc within a box, its edges included.
   *
   *  @param  arc The arc.
   *  @param  box The box.
   *  @return     The stretches of its parameter, in order, apart and longer than nothing; none for
   *              an empty box. Past 65,536 turns, its crossings are not sought.
   */
  [[nodiscard]] GEOXY_API auto Clip(const Arc &arc, const Box &box) -> std::vector<Range>;

  /**
   *  The stretches of an ellipse within a box, its edges included.
   *
   *  @param  ellipse The ellipse.
   *  @param  box     The box.
   *  @return         The stretches of its parameter, in order, apart and longer than nothing; none
   *                  for an empty box. Past 65,536 turns, its crossings are not sought.
   */
  [[nodiscard]] GEOXY_API auto Clip(const Ellipse &ellipse, const Box &box) -> std::vector<Range>;

  /**
   *  The stretches of a polyline within a box, its edges included.
   *
   *  @param  polyline The polyline.
   *  @param  box      The box.
   *  @return          The stretches of its parameter, in order and apart, a stretch running on
   *                   across vertices; one of nothing for a lone vertex within; none for no
   *                   vertices or an empty box.
   */
  [[nodiscard]] GEOXY_API auto Clip(const Polyline &polyline, const Box &box) -> std::vector<Range>;

  /**
   *  The stretches of a spline within a box, its edges included.
   *
   *  Beside a weight 2^52 times its neighbours' or more, a spline's parameter tells places apart on
   *  the curve as one double. A stretch whose ends it tells as one is nothing it can tell, and is
   *  not given; stretches on either side of a gap it cannot tell are given as one.
   *
   *  @param  spline The spline.
   *  @param  box    The box.
   *  @return        The stretches of its parameter, in order, apart and longer than nothing; none
   *                 for an empty box, or for a spline that is not valid.
   */
  [[nodiscard]] GEOXY_API auto Clip(const Spline &spline, const Box &box) -> std::vector<Range>;

  /**
   *  The stretches of a quadratic Bézier curve within a box, its edges included.
   *
   *  @param  quadratic The curve.
   *  @param  box       The box.
   *  @return           The stretches of its parameter, in order, apart and longer than nothing;
   *                    none for an empty box.
   */
  [[nodiscard]] GEOXY_API auto Clip(const Bezier::Quadratic &quadratic, const Box &box)
      -> std::vector<Range>;

  /**
   *  The stretches of a cubic Bézier curve within a box, its edges included.
   *
   *  @param  cubic The curve.
   *  @param  box   The box.
   *  @return       The stretches of its parameter, in order, apart and longer than nothing; none
   *                for an empty box.
   */
  [[nodiscard]] GEOXY_API auto Clip(const Bezier::Cubic &cubic, const Box &box)
      -> std::vector<Range>;

  /**
   *  The stretches of any curve within a box, its edges included.
   *
   *  @param  curve The curve.
   *  @param  box   The box.
   *  @return       The stretches of its parameter, as its own kind gives them.
   */
  [[nodiscard]] GEOXY_API auto Clip(const Curve &curve, const Box &box) -> std::vector<Range>;
} // namespace Geoxy

#endif /* GEOXY_CLIP_HPP */
