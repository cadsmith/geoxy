/***************************************************************************************************
 *  @file      geoxy/ends.hpp
 *  @brief     Where a curve starts and where it ends.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_ENDS_HPP
#define GEOXY_ENDS_HPP

#include "geoxy/api.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/place.hpp"

#include <optional>

namespace Geoxy
{
  /**
   *  Where a curve starts and where it ends, each as a place.
   */
  struct Endpoints
  {
    std::optional<Place> Start = {}; /**< Where it starts: none for an xline.        */
    std::optional<Place> End   = {}; /**< Where it ends: none for a ray or an xline. */

    /**
     *  Whether two pairs of ends are equal exactly.
     *
     *  @param  lhs One pair.
     *  @param  rhs The other.
     *  @return     Whether they are.
     */
    [[nodiscard]] friend auto operator==(const Endpoints &lhs, const Endpoints &rhs) -> bool
    {
      return lhs.Start == rhs.Start && lhs.End == rhs.End;
    }
  };

  /**
   *  Where a curve starts and ends: the points At gives at the ends of its domain, ENDPOINT's
   *  snaps.
   *
   *  @param  curve The curve.
   *  @return       Both ends of a bounded curve, one point twice for a closed one or a polyline of
   *                one vertex; the start alone of a ray; none of an xline, of a polyline of no
   *                vertex or of a spline that is not valid.
   */
  [[nodiscard]] GEOXY_API auto Ends(const Curve &curve) -> Endpoints;
} // namespace Geoxy

#endif /* GEOXY_ENDS_HPP */
