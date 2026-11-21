/***************************************************************************************************
 *  @file      tests/geoxy/loops.hpp
 *  @brief     Loops and regions to test with: squares, rectangles, circles, slots and polygons,
 *             moved and turned about.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TESTS_LOOPS_HPP
#define GEOXY_TESTS_LOOPS_HPP

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/line.hpp"
#include "geoxy/loop.hpp"
#include "geoxy/point.hpp"
#include "geoxy/region.hpp"
#include "geoxy/vector.hpp"

#include <cstddef>
#include <ranges>
#include <variant>
#include <vector>

namespace Geoxy::Tests
{
  /**
   *  A polygon's loop: lines from each point to the next, the last to the first.
   *
   *  @param  points The points, in order.
   *  @return        The loop.
   */
  inline auto Polygon(const std::vector<Point> &points) -> Loop
  {
    Loop loop = {};
    for (std::size_t index = 0; index < points.size(); ++index)
    {
      loop.Edges.emplace_back(
          Line{.Start = points[index], .End = points[(index + 1) % points.size()]});
    }
    return loop;
  }

  /**
   *  A rectangle's loop, counterclockwise from its least corner.
   *
   *  @param  least    Its least corner.
   *  @param  greatest Its greatest corner.
   *  @return          The loop.
   */
  inline auto Rectangle(Point least, Point greatest) -> Loop
  {
    return Polygon(
        {least, {.X = greatest.X, .Y = least.Y}, greatest, {.X = least.X, .Y = greatest.Y}});
  }

  /**
   *  A square's loop, counterclockwise from its least corner.
   *
   *  @param  least Its least corner.
   *  @param  side  Its side.
   *  @return       The loop.
   */
  inline auto Square(Point least, double side) -> Loop
  {
    return Rectangle(least, {.X = least.X + side, .Y = least.Y + side});
  }

  /**
   *  A circle's loop: one arc of a whole turn, counterclockwise.
   *
   *  @param  center The center.
   *  @param  radius The radius.
   *  @return        The loop.
   */
  inline auto Disc(Point center, double radius) -> Loop
  {
    return {.Edges = {Circle(center, radius)}};
  }

  /**
   *  A slot's loop, counterclockwise: two sides along x and a half circle at each end.
   *
   *  @param  center The center of its first half circle.
   *  @param  length How far the second's center is along x.
   *  @param  radius The half circles' radius.
   *  @return        The loop.
   */
  inline auto Slot(Point center, double length, double radius) -> Loop
  {
    const Point far = {.X = center.X + length, .Y = center.Y};
    return {
        .Edges = {
            Line{.Start = {.X = center.X, .Y = center.Y - radius},
                 .End   = {.X = far.X, .Y = far.Y - radius}},
            Arc{.Center = far, .Radius = radius, .Start = Degrees(-90), .Sweep = Degrees(180)},
            Line{.Start = {.X = far.X, .Y = far.Y + radius},
                 .End   = {.X = center.X, .Y = center.Y + radius}},
            Arc{.Center = center, .Radius = radius, .Start = Degrees(90), .Sweep = Degrees(180)},
        }};
  }

  /**
   *  A curve run the other way.
   *
   *  @param  curve A line or an arc.
   *  @return       It reversed.
   */
  inline auto Reversed(const Curve &curve) -> Curve
  {
    if (const auto *arc = std::get_if<Arc>(&curve))
    {
      return Arc{.Center = arc->Center,
                 .Radius = arc->Radius,
                 .Start  = arc->Start + arc->Sweep,
                 .Sweep  = -arc->Sweep};
    }
    const Line &line = std::get<Line>(curve);
    return Line{.Start = line.End, .End = line.Start};
  }

  /**
   *  A loop of lines and arcs run the other way.
   *
   *  @param  loop The loop.
   *  @return      It reversed.
   */
  inline auto Reversed(const Loop &loop) -> Loop
  {
    Loop reversed = {};
    for (const Curve &edge : std::views::reverse(loop.Edges))
    {
      reversed.Edges.push_back(Reversed(edge));
    }
    return reversed;
  }

  /**
   *  A curve moved.
   *
   *  @param  curve  A line or an arc.
   *  @param  offset How far.
   *  @return        It moved.
   */
  inline auto Moved(const Curve &curve, Vector offset) -> Curve
  {
    if (const auto *arc = std::get_if<Arc>(&curve))
    {
      Arc moved    = *arc;
      moved.Center = moved.Center + offset;
      return moved;
    }
    const Line &line = std::get<Line>(curve);
    return Line{.Start = line.Start + offset, .End = line.End + offset};
  }

  /**
   *  A region moved.
   *
   *  @param  region The region, of lines and arcs.
   *  @param  offset How far.
   *  @return        It moved.
   */
  inline auto Moved(const Region &region, Vector offset) -> Region
  {
    Region moved = {};
    for (const Loop &loop : region.Loops)
    {
      Loop next = {};
      for (const Curve &edge : loop.Edges)
      {
        next.Edges.push_back(Moved(edge, offset));
      }
      moved.Loops.push_back(next);
    }
    return moved;
  }
} // namespace Geoxy::Tests

#endif /* GEOXY_TESTS_LOOPS_HPP */
