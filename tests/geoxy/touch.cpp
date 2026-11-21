/***************************************************************************************************
 *  @file      tests/geoxy/touch.cpp
 *  @brief     Tests for where two curves come nearest about a touch, polished in double-double and
 *             rounded once.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "touch.hpp"

#include "crossing.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <variant>
#include <vector>

namespace Geoxy::Touches
{
  namespace
  {
    /**
     *  A map of the plane no double's rounding touches: x and y swapped, a quarter turn, a mirror
     *  across the y axis, or a scaling by 128.
     */
    enum class Map : std::uint8_t
    {
      Swap,   /**< x and y swapped.  */
      Turn,   /**< A quarter turn.   */
      Mirror, /**< x made -x.        */
      Scale,  /**< Everything * 128. */
    };

    /**
     *  A point mapped.
     *
     *  @param  map   The map.
     *  @param  point The point.
     *  @return       The point mapped.
     */
    auto Transformed(Map map, Point point) -> Point
    {
      switch (map)
      {
      case Map::Swap:
        return {.X = point.Y, .Y = point.X};
      case Map::Turn:
        return {.X = -point.Y, .Y = point.X};
      case Map::Mirror:
        return {.X = -point.X, .Y = point.Y};
      default:
        return {.X = point.X * 128.0, .Y = point.Y * 128.0};
      }
    }

    /**
     *  A point mapped back.
     *
     *  @param  map   The map.
     *  @param  point The point.
     *  @return       The point it was.
     */
    auto Back(Map map, Point point) -> Point
    {
      switch (map)
      {
      case Map::Swap:
        return {.X = point.Y, .Y = point.X};
      case Map::Turn:
        return {.X = point.Y, .Y = -point.X};
      case Map::Mirror:
        return {.X = -point.X, .Y = point.Y};
      default:
        return {.X = point.X / 128.0, .Y = point.Y / 128.0};
      }
    }

    /**
     *  An arc's start mapped, in whole degrees, exactly: turned a quarter, or reflected across the
     *  line x = y or the y axis.
     *
     *  @param  map   The map.
     *  @param  start The start, in degrees.
     *  @return       The start mapped.
     */
    auto Started(Map map, double start) -> double
    {
      switch (map)
      {
      case Map::Swap:
        return 90.0 - start;
      case Map::Turn:
        return start + 90.0;
      case Map::Mirror:
        return 180.0 - start;
      default:
        return start;
      }
    }

    /**
     *  Whether a map reflects the plane, so that a sweep or a bulge turns the other way.
     *
     *  @param  map The map.
     *  @return     Whether it does.
     */
    auto Flips(Map map) -> bool
    {
      return map == Map::Swap || map == Map::Mirror;
    }

    /**
     *  A way mapped.
     *
     *  @param  map The map.
     *  @param  way The way.
     *  @return     The way mapped.
     */
    auto Transformed(Map map, Vector way) -> Vector
    {
      return Transformed(map, Point{} + way) - Point{};
    }

    /**
     *  A line mapped.
     *
     *  @param  map   The map.
     *  @param  shape The line.
     *  @return       It mapped.
     */
    auto Shape(Map map, const Line &shape) -> Curve
    {
      return Line{.Start = Transformed(map, shape.Start), .End = Transformed(map, shape.End)};
    }

    /**
     *  A ray mapped.
     *
     *  @param  map   The map.
     *  @param  shape The ray.
     *  @return       It mapped.
     */
    auto Shape(Map map, const Ray &shape) -> Curve
    {
      return Ray{.Start = Transformed(map, shape.Start), .Way = Transformed(map, shape.Way)};
    }

    /**
     *  An endless line mapped.
     *
     *  @param  map   The map.
     *  @param  shape The endless line.
     *  @return       It mapped.
     */
    auto Shape(Map map, const Xline &shape) -> Curve
    {
      return Xline{.Through = Transformed(map, shape.Through), .Way = Transformed(map, shape.Way)};
    }

    /**
     *  An arc mapped, its angles in whole degrees.
     *
     *  @param  map   The map.
     *  @param  shape The arc.
     *  @return       It mapped.
     */
    auto Shape(Map map, const Arc &shape) -> Curve
    {
      const double sweep = shape.Sweep.InDegrees();
      return Arc{.Center = Transformed(map, shape.Center),
                 .Radius = map == Map::Scale ? shape.Radius * 128.0 : shape.Radius,
                 .Start  = Degrees(Started(map, shape.Start.InDegrees())),
                 .Sweep  = Degrees(Flips(map) ? -sweep : sweep)};
    }

    /**
     *  An ellipse mapped.
     *
     *  @param  map   The map.
     *  @param  shape The ellipse.
     *  @return       It mapped.
     */
    auto Shape(Map map, const Ellipse &shape) -> Curve
    {
      return Ellipse{.Map   = {.X      = Transformed(map, shape.Map.X),
                               .Y      = Transformed(map, shape.Map.Y),
                               .Origin = Transformed(map, shape.Map.Origin)},
                     .Start = shape.Start,
                     .Sweep = shape.Sweep};
    }

    /**
     *  A polyline mapped, its bulges turning the other way where the map reflects.
     *
     *  @param  map   The map.
     *  @param  shape The polyline.
     *  @return       It mapped.
     */
    auto Shape(Map map, const Polyline &shape) -> Curve
    {
      Polyline mapped = shape;
      for (Vertex &vertex : mapped.Vertices)
      {
        vertex.Position = Transformed(map, vertex.Position);
        vertex.Bulge    = Flips(map) ? -vertex.Bulge : vertex.Bulge;
      }
      return mapped;
    }

    /**
     *  A spline mapped.
     *
     *  @param  map   The map.
     *  @param  shape The spline.
     *  @return       It mapped.
     */
    auto Shape(Map map, const Spline &shape) -> Curve
    {
      Spline mapped = shape;
      for (Point &pole : mapped.Poles)
      {
        pole = Transformed(map, pole);
      }
      return mapped;
    }

    /**
     *  A quadratic Bézier curve mapped.
     *
     *  @param  map   The map.
     *  @param  shape The curve.
     *  @return       It mapped.
     */
    auto Shape(Map map, const Bezier::Quadratic &shape) -> Curve
    {
      return Bezier::Quadratic{.Start   = Transformed(map, shape.Start),
                               .Control = Transformed(map, shape.Control),
                               .End     = Transformed(map, shape.End)};
    }

    /**
     *  A cubic Bézier curve mapped.
     *
     *  @param  map   The map.
     *  @param  shape The curve.
     *  @return       It mapped.
     */
    auto Shape(Map map, const Bezier::Cubic &shape) -> Curve
    {
      return Bezier::Cubic{.Start = Transformed(map, shape.Start),
                           .Out   = Transformed(map, shape.Out),
                           .In    = Transformed(map, shape.In),
                           .End   = Transformed(map, shape.End)};
    }

    /**
     *  A curve mapped exactly: its points, its ways and its axes mapped, an arc's angles turned or
     *  reflected in whole degrees, a bulge's sign flipped where the map reflects.
     *
     *  @param  map   The map.
     *  @param  curve The curve, an arc's angles whole degrees.
     *  @return       The curve mapped.
     */
    auto Transformed(Map map, const Curve &curve) -> Curve
    {
      return std::visit([map](const auto &shape) -> Curve { return Shape(map, shape); }, curve);
    }

    /**
     *  The points of a meeting's touches, mapped back, in order of x then y.
     *
     *  @param  meeting The meeting.
     *  @param  back    How a point is mapped back.
     *  @return         The points.
     */
    template <typename BackT>
    auto Touches(const Intersections &meeting, const BackT &back) -> std::vector<Point>
    {
      std::vector<Point> points;
      for (const Hit &hit : meeting.Hits)
      {
        if (hit.Kind == Contact::Touch)
        {
          points.push_back(back(hit.At));
        }
      }
      std::ranges::sort(points, [](Point lhs, Point rhs) -> bool
                        { return lhs.X < rhs.X || (!(rhs.X < lhs.X) && lhs.Y < rhs.Y); });
      return points;
    }

    /**
     *  Two curves near tangent, and the tolerance they touch within.
     */
    struct Pair
    {
      Curve  First  = Line{}; /**< One curve.              */
      Curve  Second = Line{}; /**< The other.              */
      double Length = 0.0;    /**< The tolerance's length. */
    };

    /**
     *  Pairs that touch, of every kind a touch is polished on, each a hair off tangent, through it
     *  or clear of it: lines, rays and endless lines against arcs and ellipses, Bézier curves and
     *  splines against lines and arcs, a polyline's bulge against a line and its segment against an
     *  arc, ellipses against arcs.
     *
     *  @param  hair How far off tangent.
     *  @return      The pairs.
     */
    auto Touching(double hair) -> std::vector<Pair>
    {
      constexpr double Length = 1e-9;
      const double     top    = 1.0 + hair;
      const Curve      circle =
          Arc{.Center = {.X = 0.25, .Y = 0.0}, .Radius = 1.0, .Start = Degrees(0), .Sweep = Turn};
      const Curve hump    = Bezier::Quadratic{.Start   = {.X = -1.0, .Y = 0.5},
                                              .Control = {.X = 0.125, .Y = 1.5 - hair},
                                              .End     = {.X = 1.0, .Y = 0.5}};
      const Curve bump    = Bezier::Cubic{.Start = {.X = -1.0, .Y = 0.0},
                                          .Out   = {.X = -0.25, .Y = 4.0 / 3.0},
                                          .In    = {.X = 0.75, .Y = 4.0 / 3.0},
                                          .End   = {.X = 1.0, .Y = 0.0}};
      const Curve oval    = Ellipse{.Map   = {.X      = {.X = 2.0, .Y = 0.0},
                                              .Y      = {.X = 0.5, .Y = 1.0},
                                              .Origin = {.X = 0.0, .Y = -2.0}},
                                    .Start = Degrees(0),
                                    .Sweep = Turn};
      const Curve rounded = Spline{
          .Degree  = 2,
          .Knots   = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0},
          .Poles   = {{.X = -1.0, .Y = 0.0}, {.X = 0.0, .Y = 3.0 - hair}, {.X = 1.25, .Y = 0.0}},
          .Weights = {1.0, 0.5, 1.0}};
      const Curve bowed = Polyline{.Vertices = {{.Position = {.X = -2.0, .Y = 0.0}, .Bulge = 0.0},
                                                {.Position = {.X = -1.0, .Y = 0.0}, .Bulge = -1.0},
                                                {.Position = {.X = 1.0, .Y = 0.0}, .Bulge = 0.0},
                                                {.Position = {.X = 2.0, .Y = 0.0}, .Bulge = 0.0}},
                                   .Closed   = false};
      const Curve ledge = Polyline{.Vertices = {{.Position = {.X = -2.0, .Y = top}, .Bulge = 0.0},
                                                {.Position = {.X = 2.0, .Y = top}, .Bulge = 0.0},
                                                {.Position = {.X = 2.0, .Y = 3.0}, .Bulge = 0.0}},
                                   .Closed   = false};
      const Curve line  = Line{.Start = {.X = -3.0, .Y = top}, .End = {.X = 3.0, .Y = top}};
      const Curve level = Line{.Start = {.X = -2.0, .Y = 1.0}, .End = {.X = 2.0, .Y = 1.0}};
      const Curve ray   = Ray{.Start = {.X = -3.0, .Y = top}, .Way = {.X = 1.0, .Y = 0.0}};
      const Curve xline = Xline{.Through = {.X = 0.0, .Y = -top}, .Way = {.X = 0.5, .Y = 0.0}};
      const Curve under = Arc{.Center = {.X = -0.5, .Y = -3.0 - top},
                              .Radius = 1.0,
                              .Start  = Degrees(0),
                              .Sweep  = Degrees(180)};
      const Curve above = Arc{.Center = {.X = 0.25, .Y = 1.0 + top},
                              .Radius = 1.0,
                              .Start  = Degrees(180),
                              .Sweep  = Degrees(180)};
      return {
          {.First = line   , .Second = circle, .Length = Length},
          {.First = circle , .Second = ray   , .Length = Length},
          {.First = xline  , .Second = oval  , .Length = Length},
          {.First = hump   , .Second = level , .Length = Length},
          {.First = bump   , .Second = ray   , .Length = Length},
          {.First = rounded, .Second = level , .Length = Length},
          {.First = bowed  , .Second = line  , .Length = Length},
          {.First = circle , .Second = ledge , .Length = Length},
          {.First = oval   , .Second = under , .Length = Length},
          {.First = above  , .Second = circle, .Length = Length},
      };
    }
  } // namespace

  TEST_CASE("touch: the same point whichever curve is first and however the plane is mapped")
  {
    int touched = 0;
    for (const double hair : {0.0, 1e-13, -1e-13, 1e-11})
    {
      for (const Pair &pair : Touching(hair))
      {
        const Tolerance     tolerance = {.Length = pair.Length, .Radians = 1e-9};
        const Intersections met       = Intersect(pair.First, pair.Second, tolerance);
        const auto          same      = [](Point point) -> Point
        {
          return point;
        };
        const std::vector<Point> base  = Touches(met, same);
        touched                       += base.empty() ? 0 : 1;
        CHECK(Touches(Intersect(pair.Second, pair.First, tolerance), same) == base);
        for (const Map map : {Map::Swap, Map::Turn, Map::Mirror, Map::Scale})
        {
          const double        scale = map == Map::Scale ? 128.0 : 1.0;
          const Intersections mapped =
              Intersect(Transformed(map, pair.First), Transformed(map, pair.Second),
                        {.Length = pair.Length * scale, .Radians = 1e-9});
          CHECK(Touches(mapped, [map](Point point) -> Point { return Back(map, point); }) == base);
        }
      }
    }
    CHECK(touched == 40);
  }

  TEST_CASE("touch: a line and a circle touch where the circle's point nearest the line is")
  {
    /**
     *  The circle about (0.25, 0) of radius 1 comes nearest y = 1 + h at (0.25, 1), its foot on the
     *  line (0.25, 1 + h): the touch is halfway, to the bit.
     */
    for (const double hair : {0.0, 0x1p-40, -0x1p-40})
    {
      const Curve line =
          Line{.Start = {.X = -3.0, .Y = 1.0 + hair}, .End = {.X = 3.0, .Y = 1.0 + hair}};
      const Curve circle =
          Arc{.Center = {.X = 0.25, .Y = 0.0}, .Radius = 1.0, .Start = Degrees(0), .Sweep = Turn};
      const std::optional<Crossings::Crossed> touch = Nearest({&line, &circle}, {0.54, 0.26}, 1e-9);
      REQUIRE(touch.has_value());
      CHECK(touch->At.X == 0.25);
      CHECK(touch->At.Y == 1.0 + (hair / 2.0));
    }
  }

  TEST_CASE("touch: nothing where an arc is a point, the curves bend alike or come no nearer")
  {
    const Curve line  = Line{.Start = {.X = -3.0, .Y = 1.0}, .End = {.X = 3.0, .Y = 1.0}};
    const Curve point = Arc{
        .Center = {.X = 0.0, .Y = 0.0}, .Radius = 1.0, .Start = Degrees(90), .Sweep = Degrees(0)};
    CHECK_FALSE(Nearest({&point, &line}, {0.5, 0.5}, 1e-9).has_value());
    CHECK_FALSE(Nearest({&line, &point}, {0.5, 0.5}, 1e-9).has_value());
    const Curve other = Line{.Start = {.X = -3.0, .Y = 1.5}, .End = {.X = 3.0, .Y = 1.5}};
    CHECK_FALSE(Nearest({&line, &other}, {0.5, 0.5}, 1.0).has_value());
    const Curve circle =
        Arc{.Center = {.X = 0.0, .Y = 0.0}, .Radius = 0.5, .Start = Degrees(0), .Sweep = Turn};
    CHECK_FALSE(Nearest({&line, &circle}, {0.5, 0.25}, 1e-9).has_value());
    const Curve shard = Line{.Start = {.X = 2.0, .Y = 1.0}, .End = {.X = 3.0, .Y = 1.0}};
    CHECK_FALSE(Nearest({&shard, &circle}, {0.5, 0.25}, 1.0).has_value());
  }
} // namespace Geoxy::Touches
