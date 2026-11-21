/***************************************************************************************************
 *  @file      geoxy/box.hpp
 *  @brief     Axis-aligned boxes, empty until grown.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_BOX_HPP
#define GEOXY_BOX_HPP

#include "geoxy/api.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <limits>

namespace Geoxy
{
  /**
   *  An axis-aligned box, its edges included.
   *
   *  Empty until grown: its least corner at +infinity and its greatest at -infinity, so that
   *  growing it by a point gives that point.
   */
  struct Box
  {
    /**
     *  The corner of the least coordinates.
     */
    Point Min = {.X = std::numeric_limits<double>::infinity(),
                 .Y = std::numeric_limits<double>::infinity()};

    /**
     *  The corner of the greatest coordinates.
     */
    Point Max = {.X = -std::numeric_limits<double>::infinity(),
                 .Y = -std::numeric_limits<double>::infinity()};

    /**
     *  Whether both corners are equal exactly.
     *
     *  @param  lhs One box.
     *  @param  rhs The other.
     *  @return     Whether they are equal.
     */
    [[nodiscard]] friend constexpr auto operator==(const Box &lhs, const Box &rhs) -> bool
    {
      return lhs.Min == rhs.Min && lhs.Max == rhs.Max;
    }
  };

  /**
   *  Whether a box holds no point.
   *
   *  @param  box The box.
   *  @return     Whether it is empty, which a box with a NaN corner is.
   */
  [[nodiscard]] constexpr auto IsEmpty(const Box &box) -> bool
  {
    return !(box.Min.X <= box.Max.X) || !(box.Min.Y <= box.Max.Y);
  }

  /**
   *  The smallest box holding a box and a point.
   *
   *  @param  box   The box.
   *  @param  point The point.
   *  @return       The box grown to hold the point.
   */
  [[nodiscard]] constexpr auto Grown(const Box &box, Point point) -> Box
  {
    return {.Min = {.X = std::min(box.Min.X, point.X), .Y = std::min(box.Min.Y, point.Y)},
            .Max = {.X = std::max(box.Max.X, point.X), .Y = std::max(box.Max.Y, point.Y)}};
  }

  /**
   *  The smallest box holding two.
   *
   *  @param  first  One box.
   *  @param  second The other.
   *  @return        The box holding both.
   */
  [[nodiscard]] constexpr auto Joined(const Box &first, const Box &second) -> Box
  {
    const Point least = {.X = std::min(first.Min.X, second.Min.X),
                         .Y = std::min(first.Min.Y, second.Min.Y)};
    const Point most  = {.X = std::max(first.Max.X, second.Max.X),
                         .Y = std::max(first.Max.Y, second.Max.Y)};
    return {.Min = least, .Max = most};
  }

  /**
   *  Whether a box holds a point, its edges included.
   *
   *  @param  box   The box.
   *  @param  point The point.
   *  @return       Whether the point is inside or on an edge.
   */
  [[nodiscard]] constexpr auto Contains(const Box &box, Point point) -> bool
  {
    return box.Min.X <= point.X
        && point.X   <= box.Max.X
        && box.Min.Y <= point.Y
        && point.Y   <= box.Max.Y;
  }

  /**
   *  Whether two boxes share a point, an edge or a corner being enough.
   *
   *  @param  first  One box.
   *  @param  second The other.
   *  @return        Whether they overlap.
   */
  [[nodiscard]] constexpr auto Overlaps(const Box &first, const Box &second) -> bool
  {
    return first.Min.X  <= second.Max.X
        && second.Min.X <= first.Max.X
        && first.Min.Y  <= second.Max.Y
        && second.Min.Y <= first.Max.Y;
  }

  /**
   *  The width and height of a box.
   *
   *  @pre        The box is not empty.
   *  @param  box The box.
   *  @return     Its width along x and its height along y.
   */
  [[nodiscard]] constexpr auto Size(const Box &box) -> Vector
  {
    return box.Max - box.Min;
  }

  /**
   *  The middle of a box.
   *
   *  Worked from halves of its corners, so that no finite box overflows.
   *
   *  @pre        The box is not empty.
   *  @param  box The box.
   *  @return     The point halfway between its corners.
   */
  [[nodiscard]] GEOXY_API auto Center(const Box &box) -> Point;

  /**
   *  A box grown by a margin on every side.
   *
   *  @pre           The box is not empty.
   *  @param  box    The box.
   *  @param  margin The margin; a negative one shrinks it, and past nothing it is empty.
   *  @return        The box padded.
   */
  [[nodiscard]] GEOXY_API auto Padded(const Box &box, double margin) -> Box;
} // namespace Geoxy

#endif /* GEOXY_BOX_HPP */
