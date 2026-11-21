/***************************************************************************************************
 *  @file      tests/geoxy/offset.cpp
 *  @brief     Tests of offsets: the gap kinds, the invariants every offset keeps at every scale,
 *             and every curve kind.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/offset.hpp"

#include "belt/carrier.hpp"
#include "draws.hpp"
#include "parallel/bend.hpp"
#include "real.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/area.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/error.hpp"
#include "geoxy/extents.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/length.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/nearest.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <expected>
#include <initializer_list>
#include <limits>
#include <numbers>
#include <span>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  The tolerance most cases are offset within.
     */
    constexpr Tolerance OffsetFine = {.Length = 1e-9, .Radians = 1e-10};

    /**
     *  A vertex as a row: its x, its y and its bulge.
     */
    using Row = std::array<double, 3>;

    /**
     *  A polyline from rows, moved.
     *
     *  @param  rows   The vertices.
     *  @param  closed Whether it is closed.
     *  @param  by     How far it is moved.
     *  @return        The polyline.
     */
    auto Made(std::span<const Row> rows, bool closed, Vector by = {}) -> Polyline
    {
      Polyline polyline = {.Closed = closed};
      for (const Row &row : rows)
      {
        polyline.Vertices.push_back(
            {.Position = Point{.X = row[0], .Y = row[1]} + by, .Bulge = row[2]});
      }
      return polyline;
    }

    /**
     *  A polyline run the other way: its vertices in reverse, each bulge on the segment it bends.
     *
     *  @param  polyline The polyline.
     *  @return          The polyline reversed.
     */
    auto Backward(const Polyline &polyline) -> Polyline
    {
      Polyline          reversed = {.Closed = polyline.Closed};
      const std::size_t count    = polyline.Vertices.size();
      for (std::size_t index = 0; index < count; ++index)
      {
        const std::size_t from  = count - 1 - index;
        double            bulge = 0.0;
        if (polyline.Closed || from > 0)
        {
          bulge = -polyline.Vertices[(from + count - 1) % count].Bulge;
        }
        reversed.Vertices.push_back({.Position = polyline.Vertices[from].Position, .Bulge = bulge});
      }
      return reversed;
    }

    /**
     *  How far a point lies from a curve, by Nearest.
     *
     *  @param  curve The curve.
     *  @param  point The point.
     *  @return       The distance.
     */
    auto Away(const Curve &curve, Point point) -> double
    {
      return Distance(point, Nearest(curve, point).At);
    }

    /**
     *  The farthest a polyline's vertices and the points a quarter, a half and three quarters along
     *  each of its segments lie from a distance from a curve, nearer or farther.
     *
     *  @param  polyline The polyline, an offset.
     *  @param  curve    The curve offset.
     *  @param  distance The distance.
     *  @return          The most, and the least of the distances less the distance.
     */
    auto Spread(const Polyline &polyline, const Curve &curve, double distance)
        -> std::array<double, 2>
    {
      std::array<double, 2> spread = {-std::numeric_limits<double>::infinity(),
                                      std::numeric_limits<double>::infinity()};
      const std::size_t     count  = Count(polyline);
      for (std::size_t index = 0; index < count; ++index)
      {
        for (const double share : {0.0, 0.25, 0.5, 0.75})
        {
          const double off =
              Away(curve, At(polyline, static_cast<double>(index) + share)) - std::abs(distance);
          spread[0] = std::max(spread[0], off);
          spread[1] = std::min(spread[1], off);
        }
      }
      return spread;
    }

    /**
     *  Whether a polyline crosses itself: two of its segments crossing away from the vertices they
     *  share.
     *
     *  @param  polyline  The polyline.
     *  @param  tolerance The tolerance.
     *  @return           Whether it does.
     */
    auto Crossed(const Polyline &polyline, const Tolerance &tolerance) -> bool
    {
      const std::size_t count = Count(polyline);
      for (std::size_t first = 0; first < count; ++first)
      {
        for (std::size_t second = first + 1; second < count; ++second)
        {
          const Curve one =
              std::visit([](const auto &shape) -> Curve { return shape; }, Nth(polyline, first));
          const Curve other =
              std::visit([](const auto &shape) -> Curve { return shape; }, Nth(polyline, second));
          const Intersections met = Intersect(one, other, tolerance);
          for (const Hit &hit : met.Hits)
          {
            if (hit.Kind == Contact::Cross)
            {
              return true;
            }
          }
        }
      }
      return false;
    }

    /**
     *  A polyline's offset, which must be found.
     *
     *  @param  polyline   The polyline.
     *  @param  offsetting The distance and the gap kind.
     *  @param  tolerance  The tolerance.
     *  @return            The offset's polylines.
     */
    auto MustOffset(const Polyline &polyline, const Offsetting &offsetting,
                    const Tolerance &tolerance = OffsetFine) -> std::vector<Polyline>
    {
      const std::expected<std::vector<Polyline>, Error> found =
          Offset(polyline, offsetting, tolerance);
      REQUIRE(found.has_value());
      return *found;
    }

    /**
     *  The total area and length of polylines.
     *
     *  @param  polylines The polylines.
     *  @return           Their areas and their lengths summed.
     */
    auto Totals(std::span<const Polyline> polylines) -> std::array<double, 2>
    {
      std::array<double, 2> totals = {};
      for (const Polyline &polyline : polylines)
      {
        totals[0] += polyline.Closed ? Area(polyline) : 0.0;
        totals[1] += Length(polyline);
      }
      return totals;
    }

  } // namespace

  TEST_CASE("offset: a circle of two bulges offsets to the circle, to a few steps of a double")
  {
    const Polyline circle = Made(std::array<Row, 2>{{{-35, 5, 1}, {45, 5, 1}}}, true);
    for (const double distance : {-21.0, -1.0, 1.0, 21.0, 39.0})
    {
      const std::vector<Polyline> found =
          MustOffset(circle, {.Distance = distance, .Gaps = Gap::Fillet});
      REQUIRE(found.size() == 1);
      const double                radius = 40 + distance;
      const std::array<double, 2> totals = Totals(found);
      CHECK(std::abs(totals[0] - (std::numbers::pi * radius * radius)) < 1e-12 * totals[0]);
      CHECK(std::abs(totals[1] - (2 * std::numbers::pi * radius)) < 1e-13 * totals[1]);
    }
    CHECK(MustOffset(circle, {.Distance = -40}).empty());
    CHECK(MustOffset(circle, {.Distance = -41}).empty());
  }

  TEST_CASE("offset: a bulged segment whose ends are one point is read as a line, as Bent reads it")
  {
    const std::array<Vertex, 2> dot = {Vertex{.Position = {.X = 1, .Y = 2}, .Bulge = 0.5},
                                       Vertex{.Position = {.X = 1, .Y = 2}}};
    CHECK_FALSE(Parallel::Round(Parallel::Bent(dot, 0)));
    CHECK(std::holds_alternative<Line>(Bent(dot[0].Position, dot[1].Position, dot[0].Bulge)));
  }

  TEST_CASE("offset: a nearly straight arc is offset as the arc it is, however slight its bulge")
  {
    /**
     *  A bulge of a ten thousandth, a millionth or a thousand millionth puts its arc's center
     *  thousands to hundreds of millions of chords away, too far to measure from. Such an arc is
     *  measured from its ends, its radius and the way it bends, and its offset is the true one: the
     *  lengths are the truth's, worked out to 120 digits, an open L to one side and a rectangle
     *  with two such sides inward with corners extended and outward with them chamfered.
     */
    struct Bowed
    {
      double                Bulge = 0.0;
      std::array<double, 3> Long  = {};
    };

    for (const Bowed &bowed :
         {Bowed{.Bulge = 9e-5, .Long = {1447.129295203832, 2560.000491897126, 2998.822552290934}},
          Bowed{.Bulge = 1e-6, .Long = {1447.123949804513, 2560.000005400728, 2998.822510323816}},
          Bowed{.Bulge = 1e-9, .Long = {1447.123889863847, 2560.0000000054, 2998.822509939469}}})
    {
      CAPTURE(bowed.Bulge);
      const double   bulge = bowed.Bulge;
      const Polyline open =
          Made(std::array<Row, 3>{{{0, 0, bulge}, {1000, 0, 0}, {1000, 400, 0}}}, false);
      const Polyline shut =
          Made(std::array<Row, 4>{{{0, 0, bulge}, {1000, 0, -bulge}, {1000, 400, 0}, {0, 400, 0}}},
               true);
      const auto total = [](const std::vector<Polyline> &pieces) -> double
      {
        double length = 0.0;
        for (const Polyline &piece : pieces)
        {
          length += Length(piece);
        }
        return length;
      };
      CHECK(std::abs(total(MustOffset(open, {.Distance = 30, .Gaps = Gap::Fillet})) - bowed.Long[0])
            < 1e-9);
      CHECK(
          std::abs(total(MustOffset(shut, {.Distance = -30, .Gaps = Gap::Extend})) - bowed.Long[1])
          < 1e-9);
      CHECK(
          std::abs(total(MustOffset(shut, {.Distance = 30, .Gaps = Gap::Chamfer})) - bowed.Long[2])
          < 1e-9);
    }
    /**
     *  A shape that crosses itself, its long side an arc that rises 0.045 from its chord, offset by
     *  0.01: inward, the edge runs between that arc and its chord, where what the shape fills is
     *  told by the arc and not by the chord.
     */
    const Polyline crossed =
        Made(std::array<Row, 4>{{{600, 400, 0}, {0, 0, 9e-5}, {1000, 0, 0}, {0, 400, 0}}}, true);
    double inward = 0.0;
    for (const Polyline &loop : MustOffset(crossed, {.Distance = -0.01, .Gaps = Gap::Fillet}))
    {
      inward += Length(loop);
    }
    CHECK(std::abs(inward - 3397.782434513633) < 1e-9);
    double outward = 0.0;
    for (const Polyline &loop : MustOffset(crossed, {.Distance = 0.01, .Gaps = Gap::Fillet}))
    {
      outward += Length(loop);
    }
    CHECK(std::abs(outward - 3398.173474226632) < 1e-9);
    /**
     *  The arc and its chord alone, a sliver 0.045 wide: inward by 0.01 a thinner one is left,
     *  every point of it between the arc and the chord.
     */
    const Polyline sliver            = Made(std::array<Row, 2>{{{0, 0, 9e-5}, {1000, 0, 0}}}, true);
    const std::vector<Polyline> thin = MustOffset(sliver, {.Distance = -0.01, .Gaps = Gap::Fillet});
    REQUIRE(thin.size() == 1);
    CHECK(std::abs(Length(thin[0]) - 1490.711987235928) < 1e-9);
    const std::vector<Polyline> wide = MustOffset(sliver, {.Distance = 0.01, .Gaps = Gap::Fillet});
    REQUIRE(wide.size() == 1);
    CHECK(std::abs(Length(wide[0]) - 2000.062837253072) < 1e-9);
  }

  TEST_CASE("offset: how far a point is from an arc past half a turn, by the arc or by its ends")
  {
    /**
     *  Three quarters of the unit circle, from its right round by its top to its bottom: a point to
     *  its left is nearest the arc itself, one in the quarter left open nearest an end.
     */
    const std::array<Vertex, 2> wide = {
        Vertex{.Position = {.X = 1, .Y = 0}, .Bulge = std::tan(3 * Pi / 8)},
        Vertex{.Position = {.X = 0, .Y = -1}}};
    const Parallel::Bend bend = Parallel::Bent(wide, 0);
    CHECK(std::abs(Parallel::Apart(bend, {.X = -3, .Y = 0}) - 2) < 1e-12);
    CHECK(std::abs(Parallel::Apart(bend, {.X = 2, .Y = -2}) - std::sqrt(5.0)) < 1e-12);
  }

  TEST_CASE("offset: a line, a ray and an endless line move to their right")
  {
    const Line line = {.Start = {.X = 1, .Y = 2}, .End = {.X = 11, .Y = 2}};
    CHECK(Offset(line, 3).value()  == Line{.Start = {.X = 1, .Y = -1}, .End = {.X = 11, .Y = -1}});
    CHECK(Offset(line, -3).value() == Line{.Start = {.X = 1, .Y = 5}, .End = {.X = 11, .Y = 5}});
    const Ray ray = {.Start = {.X = 1, .Y = 2}, .Way = {.X = 0, .Y = 4}};
    CHECK(Offset(ray, 2).value() == Ray{.Start = {.X = 3, .Y = 2}, .Way = {.X = 0, .Y = 4}});
    const Xline xline = {.Through = {.X = 1, .Y = 2}, .Way = {.X = -5, .Y = 0}};
    CHECK(Offset(xline, 2).value() == Xline{.Through = {.X = 1, .Y = 4}, .Way = {.X = -5, .Y = 0}});
    const Line slanted = {.Start = {.X = 0, .Y = 0}, .End = {.X = 3, .Y = 4}};
    CHECK(Offset(slanted, 5).value() == Line{.Start = {.X = 4, .Y = -3}, .End = {.X = 7, .Y = 1}});
  }

  TEST_CASE("offset: a line, a ray and an endless line of no way, or not finite, are reported")
  {
    constexpr double Inf  = std::numeric_limits<double>::infinity();
    constexpr double Nan  = std::numeric_limits<double>::quiet_NaN();
    const Point      spot = {.X = 1, .Y = 2};
    CHECK(Offset(Line{.Start = spot, .End = spot}, 1).error() == Error::Degenerate);
    CHECK(Offset(Ray{.Start = spot}, 1).error()     == Error::Degenerate);
    CHECK(Offset(Xline{.Through = spot}, 1).error() == Error::Degenerate);
    CHECK(Offset(Line{.Start = spot, .End = {.X = Nan}}, 1).error()            == Error::NotFinite);
    CHECK(Offset(Line{.Start = {.Y = Inf}, .End = spot}, 1).error()            == Error::NotFinite);
    CHECK(Offset(Ray{.Start = spot, .Way = {.X = 1}}, Inf).error()             == Error::NotFinite);
    CHECK(Offset(Ray{.Start = spot, .Way = {.X = 1, .Y = Nan}}, 1).error()     == Error::NotFinite);
    CHECK(Offset(Line{.Start = {.X = -1e308}, .End = {.X = 1e308}}, 1).error() == Error::NotFinite);
    CHECK(Offset(Line{.Start = {.Y = -1.7e308}, .End = {.X = 1, .Y = -1.7e308}}, 1e308).error()
          == Error::NotFinite);
  }

  TEST_CASE("offset: an arc grows to its right, shrinks to its left, and collapses past its center")
  {
    const Arc wise = {
        .Center = {.X = 1, .Y = 2}, .Radius = 5, .Start = Degrees(30), .Sweep = Degrees(100)};
    CHECK(Offset(wise, 2).value().Radius  == 7);
    CHECK(Offset(wise, -2).value().Radius == 3);
    CHECK(Offset(wise, -2).value().Center == wise.Center);
    CHECK(Offset(wise, -2).value().Start  == wise.Start);
    CHECK(Offset(wise, -2).value().Sweep  == wise.Sweep);
    CHECK(Offset(wise, -5).error()        == Error::TooLarge);
    CHECK(Offset(wise, -6).error()        == Error::TooLarge);
    const Arc back = {
        .Center = {.X = 1, .Y = 2}, .Radius = 5, .Start = Radians(1), .Sweep = Radians(-2)};
    CHECK(Offset(back, 2).value().Radius == 3);
    CHECK(Offset(back, 6).error()        == Error::TooLarge);
    CHECK(Offset(Circle({.X = 3, .Y = 4}, 2), 1).value() == Circle({.X = 3, .Y = 4}, 3));
  }

  TEST_CASE("offset: an arc of no radius or sweep, or not finite, is reported")
  {
    constexpr double Nan  = std::numeric_limits<double>::quiet_NaN();
    const Arc        arc  = {.Radius = 5, .Sweep = Degrees(90)};
    Arc              flat = arc;
    flat.Radius           = 0;
    Arc still             = arc;
    still.Sweep           = Degrees(0);
    CHECK(Offset(flat, 1).error()  == Error::Degenerate);
    CHECK(Offset(still, 1).error() == Error::Degenerate);
    Arc lost      = arc;
    lost.Center.X = Nan;
    Arc turned    = arc;
    turned.Start  = Radians(Nan);
    Arc swept     = arc;
    swept.Sweep   = Radians(Nan);
    Arc wide      = arc;
    wide.Radius   = Nan;
    Arc huge      = arc;
    huge.Radius   = 1.7e308;
    CHECK(Offset(lost, 1).error()       == Error::NotFinite);
    CHECK(Offset(turned, 1).error()     == Error::NotFinite);
    CHECK(Offset(swept, 1).error()      == Error::NotFinite);
    CHECK(Offset(wide, 1).error()       == Error::NotFinite);
    CHECK(Offset(arc, Nan).error()      == Error::NotFinite);
    CHECK(Offset(huge, 1.7e308).error() == Error::NotFinite);
  }

  TEST_CASE("offset: level and upright lines move exactly, and their inner corners lie exactly")
  {
    /**
     *  Level and upright lines move by an exact unit normal, and their crossing is rounded once: a
     *  square's every corner is the double nearest the true one, inward and outward.
     */
    constexpr std::array<Row, 4> Rows = {{{-15, -15, 0}, {25, -15, 0}, {25, 25, 0}, {-15, 25, 0}}};
    const Polyline               square = Made(Rows, true);
    for (const double distance : {0.1, 0.7, 1.9, 3.3, 7.3})
    {
      CAPTURE(distance);
      const std::vector<Polyline> shrunk = MustOffset(square, {.Distance = -distance});
      REQUIRE(shrunk.size() == 1);
      for (const Vertex &vertex : shrunk[0].Vertices)
      {
        const auto inside = [distance](double coordinate) -> bool
        {
          return Exactly(coordinate, -15 + distance) || Exactly(coordinate, 25 - distance);
        };
        CHECK(inside(vertex.Position.X));
        CHECK(inside(vertex.Position.Y));
      }
      const std::vector<Polyline> grown =
          MustOffset(square, {.Distance = distance, .Gaps = Gap::Fillet});
      REQUIRE(grown.size() == 1);
      for (const Vertex &vertex : grown[0].Vertices)
      {
        const std::array<double, 4> sides = {-15 - distance, 25 + distance, -15, 25};
        const auto                  sided = [&sides](double coordinate) -> std::ptrdiff_t
        {
          return std::ranges::count_if(sides, [coordinate](double side) -> bool
                                       { return Exactly(side, coordinate); });
        };
        CHECK(sided(vertex.Position.X) + sided(vertex.Position.Y) == 2);
      }
    }
  }

  TEST_CASE("offset: a loop keeps an edge too short to see, and a loop too small to see is kept")
  {
    /**
     *  A random polyline's offset is a loop of four vertices, two of them 3e-8 apart: the edge
     *  between them is an edge of the true offset, and is kept, whatever tolerance is given.
     */
    constexpr std::array<Row, 3> Rows   = {{
        {0x1.28fea25747ff6p-9, 0x1.a94b96c0ebdbcp-8, -0x1.f315281281668p-5},
        {-0x1.1840abd2c5fadp-8, -0x1.9c01aaf6a68p-14, -0x1.a827999fcef33p-2},
        {0x1.b5e1c8db02d56p-8, -0x1.76cf70ebf9797p-9, 0x1.b75b42d80d0d8p+0},
    }};
    const Tolerance              coarse = {.Length = 1e-6, .Radians = 1e-10};
    const std::vector<Polyline>  loops  = MustOffset(
        Made(Rows, true), {.Distance = -0x1.3ac725835b98p-15, .Gaps = Gap::Chamfer}, coarse);
    REQUIRE(loops.size() == 1);
    CHECK(loops[0].Closed);
    CHECK(loops[0].Vertices.size() == 4);
    CHECK(std::abs(Length(loops[0]) - 0.0456111465454) < 1e-12);
    /**
     *  A circle shrunk to a millionth of a millionth of itself is still a circle: the true offset
     *  has it, and so does this one, two half turns about the center.
     */
    constexpr std::array<Row, 2> Circle = {{{-1, 0, 1}, {1, 0, 1}}};
    const std::vector<Polyline>  dot =
        MustOffset(Made(Circle, true), {.Distance = -(1 - 1e-12), .Gaps = Gap::Fillet});
    REQUIRE(dot.size()             == 1);
    REQUIRE(dot[0].Vertices.size() == 2);
    CHECK(std::abs(Distance(dot[0].Vertices[0].Position, dot[0].Vertices[1].Position) - 2e-12)
          < 1e-15);
    CHECK(std::abs(dot[0].Vertices[0].Bulge) == 1.0);
  }

  TEST_CASE("offset: a doubled segment's offset with fillets runs round both its ends")
  {
    /**
     *  A closed polyline of two vertices runs out along a segment and back: its offset either way
     *  is the stadium about it. Its moved ends at each vertex lie either side of it, in line with
     *  it but for the rounding of each, which once read as an inner corner: the ends were bridged
     *  by a line through the vertex, cut away as too near, and the loop left open or lost.
     */
    const std::array<std::array<Row, 2>, 2> segments = {{
        {{{5, 3, 0}, {2, -3, 0}}},
        {{{0x1.ea13b4272603p+0, -0x1.4ced831e018f6p+2, 0},
          {0x1.32eed2da53668p+3, -0x1.0c93bd42ea12ap+3, 0}}},
    }};
    for (const std::array<Row, 2> &rows : segments)
    {
      for (const Polyline &segment : {Made(rows, true), Backward(Made(rows, true))})
      {
        for (const double distance : {0.25, -0.25, 0x1.580b4f9e0fc85p-2})
        {
          const std::vector<Polyline> loops =
              MustOffset(segment, {.Distance = distance, .Gaps = Gap::Fillet});
          REQUIRE(loops.size() == 1);
          CHECK(loops[0].Closed);
          REQUIRE(loops[0].Vertices.size() == 4);
          const auto half = [](const Vertex &vertex) -> bool
          {
            return std::abs(std::abs(vertex.Bulge) - 1.0) < 1e-12;
          };
          CHECK(std::ranges::count_if(loops[0].Vertices, half) == 2);
        }
      }
    }
  }

  TEST_CASE("offset: the tolerance is fit to the arcs' centers, which a nearly straight arc puts "
            "far off")
  {
    /**
     *  The first segment bends by a bulge of 0.0146 over a chord of some 12,000: its center lies
     *  170,000 away, where a step of a double is 2.9e-11, past the tolerance of 1e-12 asked for.
     *  Where a point lies against that arc was then judged finer than its numbers tell, and two
     *  slices of the offset that meet end to end came out apart, the one's last vertex the other's
     *  first.
     */
    constexpr std::array<Row, 5> Rows = {{
        {-0x1.03d1d6c9c069ap+13, -0x1.0a9d6a80bbaf2p+13, 0x1.dfff34c5e84ap-7},
        {-0x1.d04fcc2758238p+9, 0x1.b871a3b856401p+11, -0x1.137387987b0b4p+1},
        {-0x1.609968e3a740fp+12, 0x1.0928c818d4619p+11, 0},
        {-0x1.1c5899ba58862p+11, 0x1.f61d3e7e7225ep+12, 0x1.a827999fcef33p-2},
        {-0x1.21d98470dc59cp+12, 0x1.0a8bd3faa3c62p+13, -0x1.a827999fcef33p-2},
    }};
    const Tolerance              fine = {.Length = 1e-12, .Radians = 1e-10};
    const std::vector<Polyline>  pieces =
        MustOffset(Made(Rows, true), {.Distance = 0x1.d989efc109dbp+6, .Gaps = Gap::Fillet}, fine);
    for (const Polyline &one : pieces)
    {
      for (const Polyline &other : pieces)
      {
        if (!one.Closed && !other.Closed && &one != &other)
        {
          CHECK_FALSE(Same(one.Vertices.back().Position, other.Vertices.front().Position, fine));
        }
      }
    }
    CHECK(pieces.size() == 5);
  }

  TEST_CASE("offset: an inner corner of two slanted lines is their moved lines' crossing, to a "
            "few steps of a double")
  {
    /**
     *  The corner of a triangle, offset inward, is where its two sides moved by the distance cross,
     *  each end moved square to its side by d / |side| of it: the double Intersect gives for the
     *  two moved lines.
     */
    constexpr std::array<Row, 3> Rows     = {{{0, 0, 0}, {17, 3, 0}, {5, 13, 0}}};
    const Polyline               triangle = Made(Rows, true);
    const auto                   corner   = [&Rows](std::size_t index) -> Point
    {
      return {.X = Rows.at(index % 3)[0], .Y = Rows.at(index % 3)[1]};
    };
    for (const double distance : {0.3, 1.1, 2.7})
    {
      CAPTURE(distance);
      const std::vector<Polyline> found = MustOffset(triangle, {.Distance = -distance});
      REQUIRE(found.size()             == 1);
      REQUIRE(found[0].Vertices.size() == 3);
      const auto moved = [distance](Point start, Point end) -> Curve
      {
        const Vector way   = end - start;
        const Vector right = Vector{.X = way.Y, .Y = -way.X} * (-distance / Length(way));
        return Line{.Start = start + right, .End = end + right};
      };
      for (std::size_t index = 0; index < 3; ++index)
      {
        const Intersections met = Intersect(moved(corner(index + 2), corner(index)),
                                            moved(corner(index), corner(index + 1)), OffsetFine);
        REQUIRE(met.Hits.size() == 1);
        CHECK(std::ranges::any_of(found[0].Vertices, [&met](const Vertex &vertex) -> bool
                                  { return Distance(vertex.Position, met.Hits[0].At) < 1e-14; }));
      }
    }
  }

  TEST_CASE("offset: a rectangle's corners extended, filleted and chamfered, out and in")
  {
    constexpr std::array<Row, 4> Rows = {{{0, 0, 0}, {20, 0, 0}, {20, 10, 0}, {0, 10, 0}}};
    const double                 cut  = 4 * 4 * (3 - (2 * std::numbers::sqrt2));
    for (const Vector by : {Vector{}, Vector{.X = 5e6, .Y = -5e6}, Vector{.X = 1e8, .Y = 1e8}})
    {
      CAPTURE(by.X);
      const Polyline  rectangle = Made(Rows, true, by);
      const Tolerance tolerance = Tolerance::Of(Extents(rectangle));
      const auto      total     = [&rectangle, &tolerance](double distance,
                                                           Gap    gaps) -> std::array<double, 2>
      {
        const std::vector<Polyline> found =
            MustOffset(rectangle, {.Distance = distance, .Gaps = gaps}, tolerance);
        REQUIRE(found.size() == 1);
        return std::array{Totals(found)[0], static_cast<double>(found[0].Vertices.size())};
      };
      const double slack = 1e-15 * 1e8 * 1e3;
      CHECK(std::abs(total(2, Gap::Extend)[0] - (24 * 14)) < slack);
      CHECK(total(2, Gap::Extend)[1] == 4);
      CHECK(std::abs(total(2, Gap::Fillet)[0] - ((20 * 10) + (60 * 2) + (4 * std::numbers::pi)))
            < slack);
      CHECK(total(2, Gap::Fillet)[1] == 8);
      CHECK(std::abs(total(2, Gap::Chamfer)[0] - ((24 * 14) - cut)) < slack);
      CHECK(total(2, Gap::Chamfer)[1] == 8);
      for (const Gap gaps : {Gap::Extend, Gap::Fillet, Gap::Chamfer})
      {
        CHECK(std::abs(total(-2, gaps)[0] - (16 * 6)) < slack);
        CHECK(total(-2, gaps)[1] == 4);
      }
    }
  }

  TEST_CASE("offset: every point of a fillet offset lies at the distance, and the loops do not "
            "cross")
  {
    Draws draws{2026};
    for (int round = 0; round < Rounds(40); ++round)
    {
      const std::size_t count = 3 + (draws.Next() % 9);
      std::vector<Row>  rows;
      for (std::size_t index = 0; index < count; ++index)
      {
        const double angle = (static_cast<double>(index) + draws.Near(0.5, 0.3))
                           * 2.0
                           * std::numbers::pi
                           / static_cast<double>(count);
        const double reach = draws.Near(10, 5);
        rows.push_back({reach * std::cos(angle), reach * std::sin(angle),
                        draws.Next() % 3 == 0 ? draws.Near(0, 0.4) : 0.0});
      }
      const bool     closed   = draws.Next() % 4 != 0;
      const Vector   by       = {.X = draws.Near(0, 1e3), .Y = draws.Near(0, 1e3)};
      const Polyline polyline = Made(rows, closed, by);
      for (const double distance : {draws.Near(2, 1.5), -draws.Near(2, 1.5)})
      {
        CAPTURE(round);
        CAPTURE(distance);
        for (const Gap gaps : {Gap::Fillet, Gap::Extend, Gap::Chamfer})
        {
          for (const Polyline &loop :
               MustOffset(polyline, {.Distance = distance, .Gaps = gaps}, OffsetFine))
          {
            const std::array<double, 2> spread = Spread(loop, polyline, distance);
            CHECK(spread[1] > -1e-8);
            if (gaps == Gap::Fillet)
            {
              CHECK(spread[0] < 1e-8);
            }
            CHECK_FALSE(Crossed(loop, OffsetFine));
          }
        }
      }
    }
  }

  TEST_CASE("offset: a convex loop's fillet offset has area A + d P + pi d^2 and length P + 2 pi d")
  {
    for (const std::size_t sides : {3U, 5U, 12U})
    {
      std::vector<Row> rows;
      for (std::size_t index = 0; index < sides; ++index)
      {
        const double angle =
            static_cast<double>(index) * 2.0 * std::numbers::pi / static_cast<double>(sides);
        rows.push_back({7 * std::cos(angle), 7 * std::sin(angle), 0.0});
      }
      const Polyline polygon = Made(rows, true);
      const double   area    = Area(polygon);
      const double   length  = Length(polygon);
      for (const double distance : {0.5, 3.0, 40.0})
      {
        const std::vector<Polyline> found =
            MustOffset(polygon, {.Distance = distance, .Gaps = Gap::Fillet});
        REQUIRE(found.size() == 1);
        const std::array<double, 2> totals = Totals(found);
        const double                round  = std::numbers::pi * distance * distance;
        CHECK(std::abs(totals[0] - (area + (distance * length) + round)) < 1e-11 * totals[0]);
        CHECK(std::abs(totals[1] - (length + (2 * std::numbers::pi * distance)))
              < 1e-12 * totals[1]);
      }
    }
  }

  TEST_CASE("offset: a hairpin, run out and back over one stretch, has neither side")
  {
    /**
     *  Beside a stretch the polyline runs over twice, one way and back, its nearest points disagree
     *  about which side a point is on: it is on neither, and so is the tip, an end of the stretch
     *  turning right back. Nothing is left, however corners are closed.
     */
    const Polyline hairpin = Made(std::array<Row, 3>{{{0, 0, 0}, {10, 0, 0}, {0, 0, 0}}}, false);
    for (const Gap gaps : {Gap::Extend, Gap::Fillet, Gap::Chamfer})
    {
      for (const double distance : {1.0, -1.0})
      {
        CAPTURE(distance);
        CHECK(MustOffset(hairpin, {.Distance = distance, .Gaps = gaps}).empty());
      }
    }
  }

  TEST_CASE("offset: an arc run back round its own circle turns at its end as its numbers tell")
  {
    /**
     *  A half circle, then a quarter of the same circle back: the two bend neither way apart, and
     *  which side of the end is outer is told by how the moved ends lie, the last steps of the
     *  numbers, so the quarter's bulge is written out: a compiler that folds the tangent and a
     *  library that works it out differ by a step. On the side they tell the offset is one chain of
     *  four vertices, round the half circle and back along the quarter; on the other the quarter's
     *  side alone, a chain of two.
     */
    const Polyline turned =
        Made(std::array<Row, 3>{{{0, 0, 1}, {0, 2, -0x1.a827999fcef33p-2}, {1, 1, 0}}}, false);
    for (const Gap gaps : {Gap::Extend, Gap::Chamfer})
    {
      for (const double distance : {0.5, -0.5})
      {
        CAPTURE(distance);
        const std::vector<Polyline> found =
            MustOffset(turned, {.Distance = distance, .Gaps = gaps});
        REQUIRE(found.size()           == 1);
        CHECK(found[0].Vertices.size() == (distance > 0.0 ? 2 : 4));
      }
    }
    /**
     *  One whose two arcs' centers come out the same double and whose moved ends do not lie in
     *  line: the corner is extended, and two arcs about one center meet nowhere, so it is filleted.
     */
    const Polyline same =
        Made(std::array<Row, 3>{{
                 {-0x1.06b22b1ed086ap+4, 0x1.661f1a591ce16p+3, 0x1.3266fbc333366p-1},
                 {-0x1.8d716e3cb4014p+4, 0x1.d63a1b47457p+4, -0x1.1c0ef078b1819p-2},
                 {-0x1.e54ee129d06e1p+3, 0x1.6c53f48b1696ep+4, 0},
             }},
             false);
    const std::vector<Polyline> found =
        MustOffset(same, {.Distance = -0x1.6c0ee5a83b03ep+1, .Gaps = Gap::Extend});
    REQUIRE(found.size() == 1);
    CHECK_FALSE(found[0].Closed);
  }

  TEST_CASE("offset: corners between lines and arcs extended, filleted and chamfered, both ways")
  {
    const std::vector<std::vector<Row>> shapes = {
        {{0, 0, 0}, {10, 0, 0.6}, {10, 8, 0}, {0, 8, 0}},
        {{0, 0, 0.3}, {10, 0, 0}, {10, 8, 0.2}, {0, 8, -0.1}},
        {{0, 0, 0.4}, {10, 0, 0.5}, {10, 8, 0.3}, {0, 8, 0.2}},
        {{0, 0, -0.3}, {10, 0, 0.7}, {12, 9, -0.2}, {3, 10, 0.45}, {-2, 4, 0}},
        {{0, 0, 0}, {6, 0, 1}, {6, 4, -1}, {6, 8, 0}, {0, 8, 0}},
        {{0, 0, 1.5}, {10, 0, 0}, {10, 8, -1.3}, {0, 8, 0}},
        {{0, 0, 0}, {10, 0, 0.05}, {20, 0, 0}, {20, 8, 0}, {0, 8, 0}},
        {{0, 0, 0}, {10, 0.5, 0}, {10.003, 0.5, 0}, {20, 1, 0}, {20, 8, 0}, {0, 8, 0}},
        {{0, 0, 0.3}, {10, 0, -0.32}, {20, 0, 0}, {20, 8, 0}, {0, 8, 0}},
    };
    for (std::size_t index = 0; index < shapes.size(); ++index)
    {
      CAPTURE(index);
      for (const bool closed : {true, false})
      {
        const Polyline forward = Made(shapes[index], closed);
        for (const Polyline &polyline : {forward, Backward(forward)})
        {
          for (const double distance : {0.7, -0.7, 1.9, -1.9})
          {
            CAPTURE(distance);
            for (const Gap gaps : {Gap::Extend, Gap::Fillet, Gap::Chamfer})
            {
              for (const Polyline &loop :
                   MustOffset(polyline, {.Distance = distance, .Gaps = gaps}))
              {
                const std::array<double, 2> spread = Spread(loop, polyline, distance);
                CHECK(spread[1] > -1e-8);
                CHECK((gaps != Gap::Fillet || spread[0] < 1e-8));
                CHECK_FALSE(Crossed(loop, OffsetFine));
              }
            }
          }
        }
      }
    }
  }

  /**
   *  An inner corner whose two pieces share a stretch, as Intersect finds one where a line or an
   *  arc lies within the tolerance of an arc, is cut where their lines or circles cross on both, or
   *  where they only graze at the start of the piece after; it was bridged by a sliver from one's
   *  end back to the other's start, run along both, and its loop came out cut in two or left open.
   *  In the first case below the arc offset by 0.1 runs 7.5e-8 past the line's start, within 1e-13
   *  of it, and the two cross 3.6e-8 along the line, at (28.499999777045722, 0.90000000000006994)
   *  by MPFR; it came out as three loops, its larger cut along a chord. A grid polyline offset
   *  inward by 1/4, two arcs at its top corner sharing 5.6e-4 at a tolerance of 1e-6, came out as
   *  two open pieces. An open polyline at 1e8 whose arcs' circles cross only past the first's end
   *  came out as two pieces.
   */
  TEST_CASE("offset: an inner corner whose pieces share a stretch is cut where they cross")
  {
    const std::vector<Row>      rows = {{27.804688, 1, 0},
                                        {28.46842055794889, 0.3429054695163245, 0},
                                        {32.34577133994935, 0.9269762697003898, 0},
                                        {32.38116957207762, 1.451312562563487, 0},
                                        {31.5, 1, -0.31783751349740424},
                                        {30.79289310940682, 1.5, 0},
                                        {29.20710689059337, 1.5, -0.31783754777018053},
                                        {28.49999981323106, 1.00000000000007, 0}};
    const std::vector<Polyline> past = MustOffset(Made(rows, true), {.Distance = -0.1});
    REQUIRE(past.size() == 2);
    const auto larger = std::ranges::find_if(past, [](const Polyline &loop) -> bool
                                             { return loop.Vertices.size() == 6; });
    REQUIRE(larger != past.end());
    CHECK(std::ranges::any_of(larger->Vertices,
                              [](const Vertex &vertex) -> bool
                              {
                                return std::abs(vertex.Position.X - 28.499999777045722)  < 1e-9
                                    && std::abs(vertex.Position.Y - 0.90000000000006994) < 1e-15;
                              }));
    const std::vector<Row> grid = {
        {5, 3, 0.5}, {-1, 9, 0.5}, {-7, -2, -1}, {-3, -6, 1}, {5, -5, 0}};
    const std::vector<Polyline> top =
        MustOffset(Made(grid, true), {.Distance = -0.25, .Gaps = Gap::Extend},
                   {.Length = 1e-6, .Radians = 1e-10});
    REQUIRE(top.size() == 1);
    CHECK(top[0].Closed);
    CHECK(top[0].Vertices.size() == 5);
    const std::vector<Row>      far = {{0x1.7d78400010625p+26, -0x1.7d783fffae148p+26, 0x1p-1},
                                       {0x1.7d783fff7ced9p+26, -0x1.7d783fffced91p+26, 0x1p-1},
                                       {0x1.7d783fffced91p+26, -0x1.7d784000624ddp+26, -0x1p+0},
                                       {0x1.7d78400010625p+26, -0x1.7d7840003126fp+26, 0x0p+0}};
    const std::vector<Polyline> graze =
        MustOffset(Made(far, false), {.Distance = -0x1.0624dd2f1a9fcp-12, .Gaps = Gap::Chamfer},
                   {.Length = 1e-9, .Radians = 1e-10});
    REQUIRE(graze.size()             == 1);
    REQUIRE(graze[0].Vertices.size() == 4);
    /**
     *  The double nearest the true vertex, as worked out to 120 digits.
     */
    CHECK(graze[0].Vertices[1].Position.X == 0x1.7d783fff80a9fp+26);
    CHECK(graze[0].Vertices[1].Position.Y == -0x1.7d783fffd0871p+26);
  }

  /**
   *  An inner corner whose pieces Intersect makes one end hit with where they cross is cut where
   *  they cross. A pentagon of radius 0.01 with bulges offset inward by 1.281e-5 at a tolerance of
   *  1e-6 (the offset corpus's 4407): two pieces cross at a sine of 1e-6 and an end of one lies
   *  within the tolerance of the other, so Intersect gives one end hit, at that end; the corner
   *  stood there, 0.14 tolerances off the offset, where every vertex lies within 4e-12 of one of
   *  it.
   */
  TEST_CASE("offset: an inner corner Intersect gives as an end hit is cut where its pieces cross")
  {
    const std::vector<Row> rows = {
        {0x1.47ae147ae147bp-7, 0x0p+0, 0x0p+0},
        {0x1.9508e5b1e9b48p-9, 0x1.37a4672bb17adp-7, 0x1.a827999fcef33p-2},
        {-0x1.091943a9eb10ep-7, 0x1.8136005b4ee8ep-8, 0x1p+0},
        {-0x1.091943a9eb10fp-7, -0x1.8136005b4ee8bp-8, 0x1.4d5b7f33c96b4p+0},
        {0x1.9508e5b1e9b43p-9, -0x1.37a4672bb17aep-7, -0x1.a98384898339cp-2}};
    const Polyline              source   = Made(rows, true);
    const double                distance = -0x1.add5338c4281bp-17;
    const Tolerance             within   = {.Length = 1e-6, .Radians = 1e-10};
    const std::vector<Polyline> found =
        MustOffset(source, {.Distance = distance, .Gaps = Gap::Fillet}, within);
    REQUIRE(found.size() == 1);
    for (const Vertex &vertex : found[0].Vertices)
    {
      const double off = Distance(vertex.Position, Nearest(source, vertex.Position).At);
      CHECK(std::abs(off - std::abs(distance)) < 1e-3 * within.Length);
    }
  }

  TEST_CASE("offset: arcs the distance takes past their centers collapse, their corners filleted")
  {
    const double   bulge   = std::tan(std::numbers::pi / 8.0);
    const Polyline rounded = Made(std::array<Row, 8>{{{-10, -5, 0},
                                                      {20, -5, bulge},
                                                      {25, 0, 0},
                                                      {25, 10, bulge},
                                                      {20, 15, 0},
                                                      {-10, 15, bulge},
                                                      {-15, 10, 0},
                                                      {-15, 0, bulge}}},
                                  true);
    for (const Gap gaps : {Gap::Extend, Gap::Fillet, Gap::Chamfer})
    {
      for (const double distance : {-5.0, -5.5, -7.0})
      {
        const std::vector<Polyline> found =
            MustOffset(rounded, {.Distance = distance, .Gaps = gaps});
        REQUIRE(found.size() == 1);
        const double inset = -distance;
        CHECK(found[0].Vertices.size() == 4);
        CHECK(std::abs(Area(found[0]) - ((40 - (2 * inset)) * (20 - (2 * inset)))) < 1e-9);
      }
    }
    const Polyline turned = Made(std::array<Row, 3>{{{0, 0, 0}, {10, 0, 1}, {10, 2, 0}}}, false);
    for (const Gap gaps : {Gap::Extend, Gap::Fillet, Gap::Chamfer})
    {
      for (const Polyline &loop : MustOffset(turned, {.Distance = -3, .Gaps = gaps}))
      {
        CHECK(Spread(loop, turned, -3)[1] > -1e-8);
      }
    }
  }

  TEST_CASE("offset: a half circle run round and back has neither side")
  {
    /**
     *  Two arcs of one circle, the second back over the first: every point beside it is nearest the
     *  polyline at a stretch run both ways, and is on neither side.
     */
    const Polyline back = Made(std::array<Row, 3>{{{10, 0, 1}, {-10, 0, -1}, {10, 0, 0}}}, false);
    for (const Gap gaps : {Gap::Extend, Gap::Chamfer})
    {
      for (const double distance : {1.0, -1.0})
      {
        CHECK(MustOffset(back, {.Distance = distance, .Gaps = gaps}).empty());
      }
    }
  }

  TEST_CASE("offset: a polyline crossing itself is offset against its other side too")
  {
    const Polyline eight =
        Made(std::array<Row, 4>{{{0, 0, 0}, {10, 10, 0}, {10, 0, 0}, {0, 10, 0}}}, true);
    for (const double distance : {1.0, -1.0, 3.0})
    {
      for (const Gap gaps : {Gap::Extend, Gap::Fillet, Gap::Chamfer})
      {
        for (const Polyline &loop : MustOffset(eight, {.Distance = distance, .Gaps = gaps}))
        {
          CHECK(Spread(loop, eight, distance)[1] > -1e-8);
          CHECK_FALSE(Crossed(loop, OffsetFine));
        }
      }
    }
  }

  TEST_CASE("offset: a polyline's vertices one point within the tolerance are taken as one")
  {
    const Polyline doubled = Made(
        std::array<Row, 6>{
            {{0, 0, 0}, {0, 1e-12, 0}, {10, 0, 0}, {10, 10, 0}, {0, 10, 0}, {0, 0, 0}}},
        true);
    const std::vector<Polyline> same = MustOffset(doubled, {.Distance = 0});
    REQUIRE(same.size()           == 1);
    CHECK(same[0].Vertices.size() == 4);
    const std::vector<Polyline> out = MustOffset(doubled, {.Distance = 1, .Gaps = Gap::Extend});
    REQUIRE(out.size() == 1);
    CHECK(std::abs(Area(out[0]) - 144) < 1e-9);
    const Polyline open = Made(std::array<Row, 3>{{{0, 0, 0}, {5, 0, 0}, {10, 0, 0}}}, false);
    CHECK(MustOffset(open, {.Distance = 0}).front() == open);
    const Polyline repeated =
        Made(std::array<Row, 5>{{{0, 0, 0}, {10, 0, 0}, {10, 10, 0}, {0, 10, 0}, {0, 0, 0}}}, true);
    CHECK(MustOffset(repeated, {.Distance = 0}).front().Vertices.size() == 4);
    const Polyline inside = Made(
        std::array<Row, 5>{{{0, 0, 0}, {0, 1e-12, 0}, {10, 0, 0}, {10, 10, 0}, {0, 10, 0}}}, true);
    CHECK(MustOffset(inside, {.Distance = 0}).front().Vertices.size() == 4);
  }

  TEST_CASE("offset: a short edge between inner corners is cut away where its neighbours cross")
  {
    const Polyline chamfered = Made(
        std::array<Row, 5>{{{0, 0, 0}, {20, 0, 0}, {20, 19, 0}, {19, 20, 0}, {0, 20, 0}}}, true);
    for (const Gap gaps : {Gap::Extend, Gap::Fillet, Gap::Chamfer})
    {
      const std::vector<Polyline> found = MustOffset(chamfered, {.Distance = -5, .Gaps = gaps});
      REQUIRE(found.size()           == 1);
      CHECK(found[0].Vertices.size() == 4);
      CHECK(std::abs(Area(found[0]) - 100) < 1e-9);
    }
  }

  TEST_CASE("offset: two mitres along one line are kept once")
  {
    const Polyline              slot  = Made(std::array<Row, 8>{{{0, 0, 0},
                                                                 {9, 0, 0},
                                                                 {9, 5, 0},
                                                                 {5, 5, 0},
                                                                 {5, 1, 0},
                                                                 {4, 1, 0},
                                                                 {4, 5, 0},
                                                                 {0, 5, 0}}},
                                             true);
    const std::vector<Polyline> found = MustOffset(slot, {.Distance = 0.6, .Gaps = Gap::Extend});
    REQUIRE(found.size() == 1);
    CHECK(std::abs(Area(found[0]) - (10.2 * 6.2)) < 1e-9);
    CHECK_FALSE(Crossed(found[0], OffsetFine));
  }

  TEST_CASE("offset: a stretch run there, back and there again is kept whole, not left open")
  {
    /**
     *  At 5e6 and 1e8 a chamfer at a corner turning by 0.0045 rad comes out 3e-7 long and turned
     *  back by the rounding of its ends, so the two pieces and the chamfer lie along each other
     *  within a step of a double: the offset runs a stretch there, back and there again. Both runs
     *  there are cut at the same points; the second was dropped as a twin of the first, and the
     *  slices stopped at the end of the run back, two open pieces where one loop, or one piece of
     *  an open polyline, belongs.
     */
    constexpr std::array<Row, 9> Loop = {{
        {0x1.312d0006be597p+22, -0x1.312cfffd66422p+22, 0x1.a827999fcef33p-2},
        {0x1.312d0002e2ffap+22, -0x1.312cfffa49d9bp+22, 0x1.d8567473fc218p-2},
        {0x1.312cffffd5633p+22, -0x1.312cfffbf6d73p+22, 0},
        {0x1.312cfffb67b0fp+22, -0x1.312cfffbd4759p+22, 0},
        {0x1.312cfff62d0f4p+22, -0x1.312cffffb7907p+22, 0},
        {0x1.312cfff7f61dp+22, -0x1.312d000427bcdp+22, 1},
        {0x1.312cfffd5cb69p+22, -0x1.312d0009334fdp+22, 0},
        {0x1.312d0005d3f5fp+22, -0x1.312d0005db8fcp+22, 0},
        {0x1.312d00079fae6p+22, -0x1.312d000523a19p+22, -1},
    }};
    const std::vector<Polyline>  loops =
        MustOffset(Made(Loop, true), {.Distance = 0x1.b52beba7a642ap-15, .Gaps = Gap::Chamfer});
    REQUIRE(loops.size() == 1);
    CHECK(loops[0].Closed);
    constexpr std::array<Row, 4> Open = {{
        {0x1.7d7840002137ap+26, -0x1.7d783fffabd8p+26, -0x1.8f619c5481e18p-3},
        {0x1.7d783fff6d914p+26, -0x1.7d783fffe5f78p+26, 0x1.a827999fcef33p-2},
        {0x1.7d783fffe3062p+26, -0x1.7d7840002957ap+26, -0x1.837379922dabp-2},
        {0x1.7d7840006b4b2p+26, -0x1.7d7840006ee89p+26, 0},
    }};
    const std::vector<Polyline>  pieces =
        MustOffset(Made(Open, false), {.Distance = -0x1.7d59aed48a73ep-16, .Gaps = Gap::Chamfer});
    REQUIRE(pieces.size() == 1);
    CHECK_FALSE(pieces[0].Closed);
  }

  TEST_CASE("offset: a piece its slices begin inside is one piece, not two that meet")
  {
    /**
     *  The slices were joined from the first not yet joined, which may lie inside a piece of the
     *  offset: the part of the piece before it, joined later, stopped where it began, two open
     *  pieces where one belongs. A closed polyline's offset inward with chamfers (the corpus's
     *  3059) was a loop and two open pieces, the second ending where the first starts; an open
     *  polyline's with mitres (2657) two pieces, the second ending where the first starts.
     */
    const auto loose = [](const std::vector<Polyline> &pieces, const Tolerance &tolerance) -> bool
    {
      return std::ranges::any_of(
          pieces,
          [&pieces, &tolerance](const Polyline &piece) -> bool
          {
            return !piece.Closed
                && std::ranges::any_of(pieces,
                                       [&piece, &tolerance](const Polyline &other) -> bool
                                       {
                                         return &other != &piece
                                             && !other.Closed
                                             && Same(piece.Vertices.back().Position,
                                                     other.Vertices.front().Position, tolerance);
                                       });
          });
    };
    constexpr std::array<Row, 11> Closed = {{
        {10, 0, -1},
        {8, 5, -0.5},
        {4, 9, 0},
        {-1, 10, 0},
        {-7, 8, -0.5},
        {-10, 3, 0},
        {-10, -3, -0.5},
        {-7, -8, 1},
        {-1, -10, 0},
        {4, -9, 1},
        {8, -5, -0.5},
    }};
    const std::vector<Polyline>   loops =
        MustOffset(Made(Closed, true), {.Distance = -3.25, .Gaps = Gap::Chamfer});
    REQUIRE(loops.size() == 1);
    CHECK(loops[0].Closed);
    CHECK(loops[0].Vertices.size() == 11);
    CHECK_FALSE(loose(loops, OffsetFine));
    constexpr std::array<Row, 8> Open   = {{
        {0, 8, -1},
        {7, 8, -0.5},
        {0, 1, 1},
        {5, -2, -0.5},
        {-5, 6, 0},
        {-9, -3, 0.5},
        {-5, -9, 0.5},
        {-3, 7, 0},
    }};
    constexpr Tolerance          Coarse = {.Length = 1e-6, .Radians = 1e-10};
    const std::vector<Polyline>  pieces =
        MustOffset(Made(Open, false), {.Distance = 4.5, .Gaps = Gap::Extend}, Coarse);
    REQUIRE(pieces.size() == 1);
    CHECK_FALSE(pieces[0].Closed);
    /**
     *  A loop that closes where an open piece begun before it starts (the corpus's 255699, at 1e8)
     *  is kept a loop, not run on through the piece.
     */
    constexpr std::array<Row, 6> Touching = {{
        {0x1.7d7840009e039p+26, -0x1.7d7840002a5ecp+26, -1},
        {0x1.7d783fffedc95p+26, -0x1.7d7840007d02p+26, 0x1.a827999fcef33p-2},
        {0x1.7d783fffdec16p+26, -0x1.7d784000188f5p+26, 0},
        {0x1.7d78400030e85p+26, -0x1.7d783fffebed9p+26, 0x1.2ceb9e22b3c48p-3},
        {0x1.7d78400000a76p+26, -0x1.7d7840004d288p+26, -0x1.70fb3877b7155p-3},
        {0x1.7d78400028bb1p+26, -0x1.7d78400034b49p+26, 0},
    }};
    const std::vector<Polyline>  touching =
        MustOffset(Made(Touching, true), {.Distance = 0x1.781d212c6c4e7p-10, .Gaps = Gap::Chamfer});
    REQUIRE(touching.size() == 1);
    CHECK(touching[0].Closed);
    CHECK(touching[0].Vertices.size() == 4);
  }

  TEST_CASE("offset: a corner whose pieces part within the tolerance turns the way they part")
  {
    /**
     *  A corner turned back past a quarter turn whose two arcs, run back from it side by side,
     *  cross back over each other having strayed less than the tolerance apart is a reversal within
     *  the tolerance, turned the way the arcs part: an outer corner, filleted. Read by its tangents
     *  as an inner corner, it was cut where the moved arcs cross, and the loop round it was left
     *  open. At 1e8 (the corpus's 257472 and 269044) each offset is one loop.
     */
    constexpr std::array<Row, 8> Wide = {{
        {0x1.7d78518p+26, -0x1.7d784p+26, -0x1.edaaef94c080cp-1},
        {0x1.7d784c5fd6a01p+26, -0x1.7d7833a0295ffp+26, 0},
        {0x1.7d784p+26, -0x1.7d782e8p+26, -1},
        {0x1.7d7833a0295ffp+26, -0x1.7d7833a0295ffp+26, -0x1.a827999fcef33p-2},
        {0x1.7d782e8p+26, -0x1.7d784p+26, 1},
        {0x1.7d7833a0295ffp+26, -0x1.7d784c5fd6a01p+26, -0x1.cf641028db872p+0},
        {0x1.7d784p+26, -0x1.7d78518p+26, 0x1.177163d4c3b58p-4},
        {0x1.7d784c5fd6a01p+26, -0x1.7d784c5fd6a01p+26, 1},
    }};
    const std::vector<Polyline>  wide =
        MustOffset(Made(Wide, true), {.Distance = 0x1.7d48533c99b53p+3, .Gaps = Gap::Fillet});
    REQUIRE(wide.size() == 1);
    CHECK(wide[0].Closed);
    constexpr std::array<Row, 13> Narrow = {{
        {0x1.7d7840005a9b5p+26, -0x1.7d783fffee795p+26, 0},
        {0x1.7d784000642ebp+26, -0x1.7d783fffbe10cp+26, 0},
        {0x1.7d78400048251p+26, -0x1.7d783fff86996p+26, 0},
        {0x1.7d783fffe8cbap+26, -0x1.7d783fff7c3bfp+26, 0x1.ccee10c52519cp-2},
        {0x1.7d783fffc75f1p+26, -0x1.7d783fffba679p+26, 0},
        {0x1.7d783fff9a59cp+26, -0x1.7d783fffb993fp+26, 0x1.ab11327f22c6p+0},
        {0x1.7d783fff8adecp+26, -0x1.7d783ffffcc1ep+26, -0x1.03901164e4d4p-1},
        {0x1.7d783fffcfda9p+26, -0x1.7d78400018c5ep+26, -0x1.6fe8c1cce92fp-3},
        {0x1.7d783fffc4834p+26, -0x1.7d78400078257p+26, 0x1.f08342ac1de7p-5},
        {0x1.7d783ffff077bp+26, -0x1.7d78400058aefp+26, -1},
        {0x1.7d78400035104p+26, -0x1.7d78400062e0ap+26, 0x1.b64adae7f1dep-2},
        {0x1.7d784000430c5p+26, -0x1.7d78400042afp+26, 0x1.97283426f79ap-1},
        {0x1.7d7840007dc8fp+26, -0x1.7d784000395edp+26, 0x1.7b5557c36242cp-2},
    }};
    constexpr Tolerance           Finest = {.Length = 1e-12, .Radians = 1e-10};
    const std::vector<Polyline>   narrow = MustOffset(
        Made(Narrow, true), {.Distance = -0x1.8abaa5d74bb0ep-11, .Gaps = Gap::Fillet}, Finest);
    REQUIRE(narrow.size() == 1);
    CHECK(narrow[0].Closed);
    /**
     *  Read by its tangents, such a corner left the offset of a small loop inward empty, though a
     *  point inside lies farther than the distance from it (the corpus's 41645).
     */
    constexpr std::array<Row, 4> Small  = {{
        {0x1.eaa8d9a420b17p-8, 0x1.3951b3a86c117p-9, -1},
        {-0x1.483adf15b0fd3p-9, 0x1.17405363a6485p-7, 0x1.87aa2c1d5da16p-3},
        {-0x1.22755a5362d5p-8, -0x1.9c43a0c7adaf6p-9, 0x1.ddd1049a4a07p-4},
        {-0x1.c1ec80fbc1106p-11, -0x1.02c85325c84d7p-7, 0x1.b272026fcb68p-2},
    }};
    constexpr Tolerance          Coarse = {.Length = 1e-6, .Radians = 1e-10};
    const Polyline               small  = Made(Small, true);
    const std::vector<Polyline>  inside =
        MustOffset(small, {.Distance = 0x1.7b4e619e31112p-10, .Gaps = Gap::Fillet}, Coarse);
    REQUIRE(inside.size() == 1);
    CHECK(inside[0].Closed);
    const std::array<double, 2> spread = Spread(inside[0], small, 0x1.7b4e619e31112p-10);
    CHECK(spread[0] <= Coarse.Length);
    CHECK(spread[1] >= -Coarse.Length);
  }

  TEST_CASE("offset: a spike turned back past a reversal by its parting, or within the tolerance")
  {
    /**
     *  A half disc whose rim, reaching its corner going down, runs back up a line turned t from the
     *  rim's way, then back to the start along the diameter, offset 0.1 to the right. The line
     *  crosses back over the rim at 2t along, having strayed t^2/2 from it.
     *
     *  - t = 1e-6, the line 1 long: the stray is 5e-13, the parting on both pieces, so the corner
     *    is a reversal within the tolerance, turned the way they part: outer, its offset run round
     *    the corner's tip by a half turn from the rim's moved end (0.9, 0).
     *  - t = 1e-11, within the angular tolerance of a reversal: outer by a half turn, the same.
     *  - t = 1e-6, the line 1e-5 long, shorter than where they part: the tangents' reading, an
     *    inner corner, no fillet.
     *
     *  Each time the offset inside the disc is one loop at the distance.
     */
    struct Spike
    {
      double Turn   = 0.0;
      double Length = 0.0;
      bool   Round  = false;
    };

    for (const Spike spike : {Spike{.Turn = 1e-6, .Length = 1, .Round = true},
                              Spike{.Turn = 1e-11, .Length = 1, .Round = true},
                              Spike{.Turn = 1e-6, .Length = 1e-5, .Round = false}})
    {
      CAPTURE(spike.Turn);
      CAPTURE(spike.Length);
      const Polyline disc = {
          .Vertices = {{.Position = {.X = -1}, .Bulge = -1},
                       {.Position = {.X = 1}},
                       {.Position = {.X = 1 - (std::sin(spike.Turn) * spike.Length),
                                     .Y = std::cos(spike.Turn) * spike.Length}}},
          .Closed   = true,
      };
      const std::vector<Polyline> pieces = MustOffset(disc, {.Distance = 0.1, .Gaps = Gap::Fillet});
      const auto                  loop   = std::ranges::find_if(pieces, &Polyline::Closed);
      REQUIRE(loop != pieces.end());
      const std::array<double, 2> spread = Spread(*loop, disc, 0.1);
      CHECK(spread[0] <= OffsetFine.Length);
      CHECK(spread[1] >= -OffsetFine.Length);
      CHECK(pieces.size() == (spike.Round ? 2U : 1U));
    }
  }

  TEST_CASE("offset: two lines turned back within the angular tolerance are filleted")
  {
    /**
     *  Two lines within the angular tolerance of a reversal are parallel within it and meet
     *  nowhere: the gap between their moved ends is closed by a fillet, as for any two parallel
     *  lines. Extended, they met where the rounding of their ends put a crossing, 180 units out
     *  (the corpus's 170862): a segment there and back, offset by 265, was an open piece with a
     *  mitre far off, not the stadium about it.
     */
    constexpr std::array<Row, 2> Stretch = {{{2000, 9000, 0}, {5000, -4000, 0}}};
    constexpr Tolerance          Finest  = {.Length = 1e-12, .Radians = 1e-10};
    const Polyline               there   = Made(Stretch, true);
    const std::vector<Polyline>  stadium =
        MustOffset(there, {.Distance = -265, .Gaps = Gap::Extend}, Finest);
    REQUIRE(stadium.size() == 1);
    CHECK(stadium[0].Closed);
    CHECK(stadium[0].Vertices.size() == 4);
    const std::array<double, 2> spread = Spread(stadium[0], there, -265);
    CHECK(spread[0] <= 1e-9);
    CHECK(spread[1] >= -1e-9);
  }

  TEST_CASE("offset: no piece kept runs on inside the distance from a cut its slice missed")
  {
    /**
     *  A slice was kept by the middles of its pieces alone, so one whose end lay deep inside the
     *  distance was kept with it (the corpus's 72912: an open piece starting 1.03 inside the
     *  distance of 2.75). An end deep inside is now cut back to where the slice leaves the
     *  distance, and a slice with both ends deep inside, a corner another piece buries, dropped:
     *  91779's mitre, poking out across the semicircle 0.3 from its corner, ended 0.36 and 0.48
     *  inside the distance of 0.75; 226583 and 161385, whose truth is empty, each kept a mitre
     *  whose ends lay 74.8 and 2.67 inside. In 15747, at 1e8, a vertex between a slice's ends lay
     *  637 tolerances inside: the slice is dropped. In 78420, crossing itself 174 times, a slice on
     *  the distance ended 180 and 200 tolerances inside at cuts made inexact by a shallow crossing,
     *  not a buried corner: both ends are cut back.
     */
    struct Case
    {
      std::vector<Row> Rows     = {};  /**< The closed polyline.            */
      double           Distance = 0.0; /**< The distance.                   */
      double           Length   = 0.0; /**< The tolerance.                  */
      std::size_t      Pieces   = 0;   /**< How many pieces the offset has. */
    };

    const std::array<Case, 6> cases = {{
        {.Rows     = {{-1, 7, -0.5}, {-6, -4, 0}, {4, -9, -1}},
         .Distance = 0.75,
         .Length   = 1e-6,
         .Pieces   = 2},
        {.Rows     = {{130.45600000000002, -121.456, 1},
                      {123.456, -119.456, 0},
                      {119.456, -121.456, -1},
                      {121.456, -125.456, 0},
                      {123.456, -130.45600000000002, 0},
                      {132.45600000000002, -126.456, 0}},
         .Distance = -2.75,
         .Length   = 1e-9,
         .Pieces   = 1},
        {.Rows     = {{1288.2051641469898, -1078.0136590539983, 0},
                      {1633.0871887451676, 536.52931850340212, 0},
                      {-7780.4392934489033, -2710.1281768022404, -1.8294343361861487}},
         .Distance = 8826.7570140114367,
         .Length   = 1e-9,
         .Pieces   = 0},
        {.Rows     = {{-1.0467550624178834, 0.90888154411338817, 0},
                      {-9.583073847069608, -2.7633770827240456, 0},
                      {0.66309981137412066, 3.3959949793994326, -0.41421356237309509},
                      {-6.1071358392163821, -8.6377406978796429, 0}},
         .Distance = 12.040871539067405,
         .Length   = 1e-12,
         .Pieces   = 0},
        {.Rows     = {{100000000.009, -99999999.997999996, 0},
                      {99999999.991999999, -100000000.006, -1},
                      {99999999.993000001, -100000000.008, -0.5},
                      {99999999.990999997, -100000000.005, 0.5},
                      {99999999.997999996, -99999999.994000003, 0}},
         .Distance = -0.00025000000000000001,
         .Length   = 0x1p-22,
         .Pieces   = 1},
        {.Rows     = {{-7609.7672156796643, 5770.5756131008438, 0},
                      {-9687.3812255237699, 8748.8972794820056, 0},
                      {-1128.1000139833495, -3230.6226710682413, -1},
                      {5478.8494376762028, 4846.089974862979, -0.021654568996478663},
                      {-1497.0852811497793, -5685.986740843141, 0},
                      {2360.6452209936178, -9252.1822680125824, 0},
                      {1072.2111355043307, -5753.0151940552623, 0.049644461895442638},
                      {1168.6317787398277, 9491.9161525592735, -1},
                      {9627.2542419478996, 9560.4643562125566, -1},
                      {6477.9296507953286, -5923.8596098332791, -1.2551229165543665},
                      {7029.3691183823912, -9789.5311266438675, 0},
                      {-2760.3946084687127, 7119.809524566911, 0},
                      {792.7333762012729, 7958.8673365019422, 0.41421356237309509},
                      {-6827.6399974897104, -2650.5979228059555, 1},
                      {2899.6835816363741, 3324.3251352607235, -1.8397851137107475},
                      {3402.8177481688822, -6044.5688770585339, 0},
                      {-3113.3436107756811, -8594.1071600474061, -0.14457371859593504},
                      {5769.6897794051165, 7813.4355177648767, 0.82771315262453182},
                      {-5950.5923737720695, -9599.6420599006797, -1.1802484505857302e-05},
                      {-6117.482657826954, 6608.054581562129, 1},
                      {9870.964079307083, -9618.9513422547607, -0.60448508627725106},
                      {-9861.2948959166297, 5652.9427416186381, 0.41421356237309509},
                      {-6706.5920943109541, 5188.981586783395, 1},
                      {4795.8925600840203, -2381.4560725936394, 2.0978111741584176},
                      {2210.9119578026953, 8245.897184493977, 1},
                      {-20.710551669932187, 4277.2888685199641, 1.3982960464410734},
                      {5751.0382031935578, 6352.5388083648768, 0},
                      {-5537.796564009067, -4498.5390398445015, 1},
                      {6556.9553152752924, -4156.9901829820728, -2.1341299738365427},
                      {-3981.9933937359056, -9777.0838922123348, 0}},
         .Distance = 18.774811457189688,
         .Length   = 0x1p-20,
         .Pieces   = 59},
    }};
    for (const Case &each : cases)
    {
      const Polyline              source = Made(each.Rows, true);
      const std::vector<Polyline> pieces =
          MustOffset(source, {.Distance = each.Distance, .Gaps = Gap::Extend},
                     {.Length = each.Length, .Radians = 1e-10});
      CHECK(pieces.size() == each.Pieces);
      for (const Polyline &piece : pieces)
      {
        /**
         *  At its vertices: between them At reads an arc of a bulge under a millionth as a line.
         */
        for (const Vertex &vertex : piece.Vertices)
        {
          CHECK(Away(source, vertex.Position) - std::abs(each.Distance) >= -2.0 * each.Length);
        }
      }
    }
  }

  TEST_CASE("offset: a segment run over more than once is offset as once")
  {
    /**
     *  A closed polyline running along one segment there, back, there and back again, its offset
     *  run twice, and one running there, back past its start and part of the way there: each is
     *  offset as the segment is, round at each end: a line back over itself turns neither way, so
     *  neither side of its end is outer and no chamfer closes it.
     */
    constexpr std::array<Row, 6> Twice = {{
        {0, 0, 0},
        {2, 0, 0},
        {2, 0, 0},
        {0, 0, -1},
        {0, 0, 0},
        {2, 0, 0},
    }};
    const std::vector<Polyline>  loop =
        MustOffset(Made(Twice, true), {.Distance = -0x1.6666666666667p-1, .Gaps = Gap::Chamfer});
    REQUIRE(loop.size() == 1);
    CHECK(loop[0].Closed);
    CHECK(std::abs(std::abs(Area(loop[0])) - ((2.0 * 1.4) + (Pi * 0.49))) < 1e-9);
    constexpr std::array<Row, 6> Past = {{
        {0, 0, 0},
        {-2, 0, 0},
        {-2, 0, 0},
        {2, 0, 0},
        {2, 0, 0},
        {-1, 0, 0},
    }};
    const std::vector<Polyline>  past =
        MustOffset(Made(Past, true), {.Distance = 0x1.3333333333334p-1, .Gaps = Gap::Chamfer});
    REQUIRE(past.size() == 1);
    CHECK(past[0].Closed);
    CHECK(std::abs(std::abs(Area(past[0])) - ((4.0 * 1.2) + (Pi * 0.36))) < 1e-9);
  }

  TEST_CASE("offset: a polyline of no extent, or not finite, is reported")
  {
    constexpr double Nan = std::numeric_limits<double>::quiet_NaN();
    constexpr double Inf = std::numeric_limits<double>::infinity();
    const Offsetting by  = {.Distance = 1};
    CHECK(Offset(Polyline{}, by, OffsetFine).error() == Error::Degenerate);
    CHECK(Offset(Made(std::array<Row, 1>{{{1, 2, 0}}}, true), by, OffsetFine).error()
          == Error::Degenerate);
    CHECK(Offset(Made(std::array<Row, 2>{{{1, 2, 0}, {1, 2, 0}}}, false), by, OffsetFine).error()
          == Error::Degenerate);
    CHECK(Offset(Made(std::array<Row, 2>{{{1, 2, 0}, {1, 2, 0}}}, false), {}, OffsetFine).error()
          == Error::Degenerate);
    CHECK(Offset(Made(std::array<Row, 2>{{{1, 2, 0}, {1, 2, 0}}}, true), by, OffsetFine).error()
          == Error::Degenerate);
    CHECK(Offset(Made(std::array<Row, 2>{{{1, Nan, 0}, {1, 2, 0}}}, false), by, OffsetFine).error()
          == Error::NotFinite);
    CHECK(Offset(Made(std::array<Row, 2>{{{1, 3, Inf}, {1, 2, 0}}}, false), by, OffsetFine).error()
          == Error::NotFinite);
    CHECK(Offset(Made(std::array<Row, 2>{{{1, 3, 0}, {1, 2, 0}}}, false), {.Distance = Nan},
                 OffsetFine)
              .error()
          == Error::NotFinite);
    CHECK(
        Offset(Made(std::array<Row, 2>{{{1e300, 3, 0}, {1, 2, 0}}}, false), by, OffsetFine).error()
        == Error::NotFinite);
  }

  TEST_CASE("offset: an ellipse's offset lies within the accuracy of the curve at the distance")
  {
    const Ellipse circle = {.Map   = {.X = {.X = 10}, .Y = {.Y = 10}, .Origin = {.X = 1, .Y = 2}},
                            .Start = Radians(0),
                            .Sweep = Radians(2 * Pi)};
    for (const double distance : {2.0, -2.0})
    {
      const std::vector<Polyline> found =
          Offset(circle, {.Distance = distance, .Approximation = {1e-6}}, OffsetFine).value();
      REQUIRE(found.size() == 1);
      const double radius = 10 + distance;
      CHECK(std::abs(Area(found[0]) - (std::numbers::pi * radius * radius)) < 1e-9);
    }
    const Ellipse oval      = {.Map   = {.X = {.X = 6, .Y = 8}, .Y = {.X = -4, .Y = 3}},
                               .Start = Degrees(0),
                               .Sweep = Degrees(360)};
    const double  perimeter = Length(oval);
    for (const double accuracy : {1e-3, 1e-7})
    {
      CAPTURE(accuracy);
      const std::vector<Polyline> out =
          Offset(oval, {.Distance = 2, .Approximation = {accuracy}}, OffsetFine).value();
      REQUIRE(out.size() == 1);
      const double steiner = (std::numbers::pi * 50) + (2 * perimeter) + (4 * std::numbers::pi);
      CHECK(std::abs(Area(out[0]) - steiner) < 100 * accuracy);
      for (const double distance : {2.0, -2.0, -3.5, 1e-4})
      {
        CAPTURE(distance);
        for (const Polyline &loop :
             Offset(oval, {.Distance = distance, .Approximation = {accuracy}}, OffsetFine).value())
        {
          const std::array<double, 2> spread = Spread(loop, oval, distance);
          CHECK(spread[0] < accuracy);
          CHECK(spread[1] > -accuracy);
          CHECK_FALSE(Crossed(loop, OffsetFine));
        }
      }
    }
    const Ellipse half = {
        .Map = {.X = {.X = 10}, .Y = {.Y = 5}}, .Start = Radians(0), .Sweep = Radians(-Pi)};
    for (const double distance : {3.0, -3.0})
    {
      for (const Polyline &loop :
           Offset(half, {.Distance = distance, .Approximation = {1e-6}}, OffsetFine).value())
      {
        CHECK_FALSE(loop.Closed);
        const std::array<double, 2> spread = Spread(loop, half, distance);
        CHECK(spread[0] < 1e-6);
        CHECK(spread[1] > -1e-6);
      }
      /**
       *  Asked no nearer than its own size, the half ellipse is laid as its chord, and the offset
       *  is that chord moved.
       */
      const std::vector<Polyline> coarse =
          Offset(half, {.Distance = distance, .Approximation = {1e3}}, OffsetFine).value();
      REQUIRE(coarse.size()           == 1);
      CHECK(coarse[0].Vertices.size() == 2);
      CHECK(std::abs(Length(coarse[0]) - 20) < 1e-12);
    }
  }

  TEST_CASE("offset: an ellipse of no extent, or not finite, or needing too many arcs, is reported")
  {
    constexpr double Nan  = std::numeric_limits<double>::quiet_NaN();
    const Ellipse    flat = {.Map   = {.X = {.X = 3, .Y = 1}, .Y = {.X = 6, .Y = 2}},
                             .Sweep = Radians(1)};
    CHECK(Offset(flat, {.Distance = 1}, OffsetFine).error() == Error::Degenerate);
    const Ellipse still = {.Map = {.X = {.X = 3}, .Y = {.Y = 2}}, .Sweep = Radians(0)};
    CHECK(Offset(still, {.Distance = 1}, OffsetFine).error() == Error::Degenerate);
    Ellipse lost = still;
    lost.Sweep   = Radians(1);
    lost.Map.Y.X = Nan;
    CHECK(Offset(lost, {.Distance = 1}, OffsetFine).error() == Error::NotFinite);
    Ellipse turned = still;
    turned.Sweep   = Radians(1);
    turned.Start   = Radians(Nan);
    CHECK(Offset(turned, {.Distance = 1}, OffsetFine).error() == Error::NotFinite);
    CHECK(Offset(Ellipse{.Map = {.X = {.X = 3}, .Y = {.Y = 2}}, .Sweep = Radians(1)},
                 {.Distance = Nan}, OffsetFine)
              .error()
          == Error::NotFinite);
    const Ellipse thin = {.Map = {.X = {.X = 1e3}, .Y = {.Y = 1}}, .Sweep = Radians(2 * Pi)};
    CHECK(Offset(thin, {.Distance = 0.5}, {}).error() == Error::Spent);
  }

  TEST_CASE("offset: a spline's offset lies within the accuracy, its kinks closed by the gap kind")
  {
    const Spline wave = {.Degree = 3,
                         .Knots  = {0, 0, 0, 0, 1, 2, 2, 2, 2},
                         .Poles  = {{.X = 0, .Y = 0},
                                    {.X = 3, .Y =  5},
                                    {.X = 6, .Y = -2},
                                    {.X = 9, .Y =  4},
                                    {.X = 12, .Y = 0}}};
    const Spline kink = {.Degree = 2,
                         .Knots  = {0, 0, 0, 1, 1, 2, 2, 2},
                         .Poles  = {{.X = 0, .Y = 0},
                                    {.X = 2, .Y = 3},
                                    {.X = 4, .Y = 0},
                                    {.X = 6, .Y = 3},
                                    {.X = 8, .Y = 0}}};
    const Spline loop = {.Degree  = 2,
                         .Knots   = {0, 0, 0, 1, 1, 2, 2, 2},
                         .Poles   = {{.X = 0, .Y = 0},
                                     {.X =  5, .Y =  5},
                                     {.X = 10, .Y =  0},
                                     {.X =  5, .Y = -5},
                                     {.X = 0, .Y = 0}},
                         .Weights = {1, 0.5, 1, 0.5, 1}};
    for (const Spline *spline : {&wave, &kink, &loop})
    {
      for (const double distance : {0.5, -0.5, 3.0})
      {
        CAPTURE(distance);
        for (const Gap gaps : {Gap::Extend, Gap::Fillet, Gap::Chamfer})
        {
          for (const Polyline &found :
               Offset(*spline, {.Distance = distance, .Gaps = gaps, .Approximation = {1e-6}},
                      OffsetFine)
                   .value())
          {
            const std::array<double, 2> spread = Spread(found, *spline, distance);
            CHECK(spread[1] > -1e-6);
            CHECK((gaps != Gap::Fillet || spread[0] < 1e-6));
            CHECK(found.Closed == (spline == &loop));
          }
        }
      }
    }
  }

  TEST_CASE("offset: a closed spline crossing itself, and one standing still on a span")
  {
    const Spline eight = {.Degree = 2,
                          .Knots  = {0, 0, 0, 1, 2, 3, 4, 5, 5, 5},
                          .Poles  = {{.X = 0, .Y = 0},
                                     {.X =  4, .Y =  4},
                                     {.X =  8, .Y = -4},
                                     {.X = 12, .Y =  0},
                                     {.X =  8, .Y =  4},
                                     {.X =  4, .Y = -4},
                                     {.X = 0, .Y = 0}}};
    for (const double distance : {0.4, -0.4})
    {
      const std::vector<Polyline> found =
          Offset(eight, {.Distance = distance, .Approximation = {1e-6}}, OffsetFine).value();
      REQUIRE_FALSE(found.empty());
      for (const Polyline &loop : found)
      {
        CHECK(Spread(loop, eight, distance)[1] > -1e-6);
      }
    }
    const Spline still = {
        .Degree = 1,
        .Knots  = {0, 0, 1, 2, 3, 3},
        .Poles  = {{.X = 0, .Y = 0}, {.X = 5, .Y = 0}, {.X = 5, .Y = 0}, {.X = 5, .Y = 5}}};
    const std::vector<Polyline> bent =
        Offset(still, {.Distance = 1, .Gaps = Gap::Extend, .Approximation = {1e-6}}, OffsetFine)
            .value();
    REQUIRE(bent.size() == 1);
    CHECK(std::abs(Length(bent[0]) - 12) < 1e-6);
  }

  TEST_CASE("offset: a spline that is not valid, of no extent, or not finite, is reported")
  {
    constexpr double Nan     = std::numeric_limits<double>::quiet_NaN();
    const Spline     invalid = {
        .Degree = 3, .Knots = {0, 1}, .Poles = {{.X = 0, .Y = 0}, {.X = 1, .Y = 1}}};
    CHECK(Offset(invalid, {.Distance = 1}, OffsetFine).error() == Error::Degenerate);
    const Spline point = {
        .Degree = 1, .Knots = {0, 0, 1, 1}, .Poles = {{.X = 2, .Y = 3}, {.X = 2, .Y = 3}}};
    CHECK(Offset(point, {.Distance = 1}, OffsetFine).error() == Error::Degenerate);
    const Spline lost = {
        .Degree = 1, .Knots = {0, 0, 1, 1}, .Poles = {{.X = 2, .Y = Nan}, {.X = 2, .Y = 3}}};
    CHECK(Offset(lost, {.Distance = 1}, OffsetFine).error() == Error::NotFinite);
    const Spline line = {
        .Degree = 1, .Knots = {0, 0, 1, 1}, .Poles = {{.X = 2, .Y = 3}, {.X = 5, .Y = 3}}};
    CHECK(Offset(line, {.Distance = Nan}, OffsetFine).error() == Error::NotFinite);
  }

  TEST_CASE("offset: any curve, each as its own kind gives it")
  {
    const Tolerance  tolerance = OffsetFine;
    const Offsetting by        = {.Distance = 1, .Approximation = {1e-6}};
    const Line       line      = {.End = {.X = 4}};
    CHECK(std::get<Line>(Offset(Curve{line}, by, tolerance).value().front())
          == Offset(line, 1).value());
    const Ray ray = {.Way = {.X = 4}};
    CHECK(std::get<Ray>(Offset(Curve{ray}, by, tolerance).value().front())
          == Offset(ray, 1).value());
    const Xline xline = {.Way = {.X = 4}};
    CHECK(std::get<Xline>(Offset(Curve{xline}, by, tolerance).value().front())
          == Offset(xline, 1).value());
    const Arc arc = Circle({}, 2);
    CHECK(std::get<Arc>(Offset(Curve{arc}, by, tolerance).value().front()) == Circle({}, 3));
    CHECK(Offset(Curve{arc}, {.Distance = -3}, tolerance).value().empty());
    CHECK(Offset(Curve{Arc{}}, by, tolerance).error()  == Error::Degenerate);
    CHECK(Offset(Curve{Line{}}, by, tolerance).error() == Error::Degenerate);
    const Polyline square =
        Made(std::array<Row, 4>{{{0, 0, 0}, {4, 0, 0}, {4, 4, 0}, {0, 4, 0}}}, true);
    CHECK(std::get<Polyline>(Offset(Curve{square}, by, tolerance).value().front()).Vertices.size()
          == 4);
    const Ellipse ellipse = {.Map = {.X = {.X = 3}, .Y = {.Y = 2}}, .Sweep = Radians(2 * Pi)};
    CHECK(Offset(Curve{ellipse}, by, tolerance).value().size() == 1);
    const Spline spline = {
        .Degree = 1, .Knots = {0, 0, 1, 1}, .Poles = {{.X = 2, .Y = 3}, {.X = 5, .Y = 3}}};
    CHECK(Offset(Curve{spline}, by, tolerance).value().size() == 1);
    const Bezier::Quadratic quadratic = {.Start = {}, .Control = {.X = 2, .Y = 3}, .End = {.X = 4}};
    const Bezier::Cubic     cubic     = {
        .Start = {}, .Out = {.X = 1, .Y = 3}, .In = {.X = 3, .Y = 3}, .End = {.X = 4}};
    for (const Curve &curve : {Curve{quadratic}, Curve{cubic}})
    {
      const std::vector<Curve> found = Offset(curve, by, tolerance).value();
      REQUIRE(found.size() == 1);
      CHECK(Spread(std::get<Polyline>(found.front()), curve, 1)[0] < 1e-6);
    }
    const Bezier::Cubic stalled = {
        .Start = {}, .Out = {}, .In = {.X = 4, .Y = 3}, .End = {.X = 4, .Y = 3}};
    const std::vector<Curve> straight = Offset(Curve{stalled}, by, tolerance).value();
    REQUIRE(straight.size() == 1);
    CHECK(std::abs(Length(std::get<Polyline>(straight.front())) - 5) < 1e-6);
    const Bezier::Cubic closed = {
        .Start = {}, .Out = {.X = 4, .Y = 4}, .In = {.X = -4, .Y = 4}, .End = {}};
    CHECK(std::get<Polyline>(Offset(Curve{closed}, by, tolerance).value().front()).Closed);
    const Bezier::Quadratic dot = {.Start = {.X = 1}, .Control = {.X = 1}, .End = {.X = 1}};
    CHECK(Offset(Curve{dot}, by, tolerance).error() == Error::Degenerate);
    const Bezier::Quadratic lost = {.Start = {.X = std::numeric_limits<double>::quiet_NaN()},
                                    .End   = {.X = 1}};
    CHECK(Offset(Curve{lost}, by, tolerance).error() == Error::NotFinite);
  }

  namespace
  {
    /**
     *  How many vertices each loop of a region's offset has, fewest first.
     *
     *  @param  rows      The closed polyline's vertices, each its x, its y and its bulge.
     *  @param  buffering The distance, the side and the filling.
     *  @return           The counts.
     */
    auto Loops(std::initializer_list<Row> rows, const Buffering &buffering)
        -> std::vector<std::size_t>
    {
      Polyline polyline = {.Closed = true};
      for (const Row &row : rows)
      {
        polyline.Vertices.push_back({.Position = {.X = row[0], .Y = row[1]}, .Bulge = row[2]});
      }
      std::vector<std::size_t> counts;
      for (const Polyline &loop : Buffer(polyline, buffering, OffsetFine).value())
      {
        CHECK(loop.Closed);
        counts.push_back(loop.Vertices.size());
      }
      std::ranges::sort(counts);
      return counts;
    }

    /**
     *  A list of counts.
     */
    using Counts = std::vector<std::size_t>;
  } // namespace

  TEST_CASE("buffer: a square shrinks, grows with round corners, and both at once")
  {
    const std::initializer_list<Row> square = {{0, 0, 0}, {10, 0, 0}, {10, 10, 0}, {0, 10, 0}};
    CHECK(Loops(square, {.Distance = 1, .Sides = Side::Inward})  == Counts{4});
    CHECK(Loops(square, {.Distance = 1, .Sides = Side::Outward}) == Counts{8});
    CHECK(Loops(square, {.Distance = 1, .Sides = Side::Both})    == Counts{4, 8});
    CHECK(Loops(square, {.Distance = 6, .Sides = Side::Inward}).empty());
    /**
     *  Run the other way round, its inside is on its right: the same answers.
     */
    const std::initializer_list<Row> clockwise = {{0, 0, 0}, {0, 10, 0}, {10, 10, 0}, {10, 0, 0}};
    CHECK(Loops(clockwise, {.Distance = 1, .Sides = Side::Inward})  == Counts{4});
    CHECK(Loops(clockwise, {.Distance = 1, .Sides = Side::Outward}) == Counts{8});
    CHECK(Loops(clockwise, {.Distance = 1, .Sides = Side::Both})    == Counts{4, 8});
    Polyline polyline = {.Closed = true};
    for (const Row &row : square)
    {
      polyline.Vertices.push_back({.Position = {.X = row[0], .Y = row[1]}});
    }
    const std::vector<Polyline> grown =
        Buffer(polyline, {.Distance = 1, .Sides = Side::Outward}, OffsetFine).value();
    REQUIRE(grown.size()      == 1);
    CHECK(Area(grown.front()) == doctest::Approx(140 + Pi));
    const std::vector<Polyline> band =
        Buffer(polyline, {.Distance = 1, .Sides = Side::Both}, OffsetFine).value();
    REQUIRE(band.size() == 2);
    /**
     *  The region on the left of each loop: its outline counterclockwise, its hole clockwise.
     */
    CHECK(Area(band[0]) + Area(band[1]) == doctest::Approx(140 + Pi - 64));
    polyline.Vertices.insert(polyline.Vertices.begin() + 1, polyline.Vertices[1]);
    CHECK(Buffer(polyline, {.Distance = 1, .Sides = Side::Inward}, OffsetFine).value().size() == 1);
  }

  TEST_CASE("buffer: where two parts lie twice the distance apart the stretch between "
            "them is no edge")
  {
    /**
     *  A U whose bottom is 2 thick, offset inward by 1: the bottom's two sides moved meet along y =
     *  1 from x = 4 to x = 6, the region on neither side, and the answer is two loops.
     */
    const std::initializer_list<Row> u = {{0, 0, 0}, {10, 0, 0}, {10, 10, 0}, {6, 10, 0},
                                          {6, 2, 0}, {4, 2, 0},  {4, 10, 0},  {0, 10, 0}};
    CHECK(Loops(u, {.Distance = 1, .Sides = Side::Inward})  == Counts{5, 5});
    CHECK(Loops(u, {.Distance = 1, .Sides = Side::Outward}) == Counts{11});
    CHECK(Loops(u, {.Distance = 1, .Sides = Side::Both})    == Counts{5, 5, 11});
  }

  TEST_CASE("buffer: a star that crosses itself is offset as it is filled")
  {
    /**
     *  A five-pointed star drawn in one stroke winds twice about its middle: filled there by
     *  nonzero winding, not by even-odd.
     */
    const std::initializer_list<Row> star = {
        {0, 10, 0}, {-5.878, -8.09, 0}, {9.511, 3.09, 0}, {-9.511, 3.09, 0}, {5.878, -8.09, 0}};
    const auto asked = [&star](Side side, Fills fill) -> Counts
    {
      return Loops(star, {.Distance = 0.3, .Sides = side, .Fill = fill});
    };
    CHECK(asked(Side::Inward, Fills::NonZero)  == Counts{3, 3, 3, 3, 3, 5});
    CHECK(asked(Side::Inward, Fills::EvenOdd)  == Counts{3, 3, 3, 3, 3});
    CHECK(asked(Side::Outward, Fills::NonZero) == Counts{15});
    CHECK(asked(Side::Outward, Fills::EvenOdd) == Counts{5, 15});
    CHECK(asked(Side::Both, Fills::NonZero)    == Counts{3, 3, 3, 3, 3, 5, 15});
    CHECK(asked(Side::Both, Fills::EvenOdd)    == Counts{3, 3, 3, 3, 3, 5, 15});
  }

  TEST_CASE("buffer: a stretch reached from both sides of a polyline that runs back over "
            "itself is kept once")
  {
    /**
     *  Up a line, round a half circle, and straight back down over the line: the line's two sides
     *  are the same stretches, and inward nothing is filled but the half disk, too thin to keep.
     */
    const std::initializer_list<Row> back = {{0, 0, 0}, {0, 1, 1}, {0, 2, 0}};
    CHECK(Loops(back, {.Distance = 1, .Sides = Side::Outward}) == Counts{5});
    CHECK(Loops(back, {.Distance = 1, .Sides = Side::Both})    == Counts{5});
    CHECK(Loops(back, {.Distance = 1, .Sides = Side::Inward}).empty());
    /**
     *  Three arcs, two of them meeting at a cusp: the one-sided offset leaves the first hair of the
     *  cusp's arc, which leads to no loop.
     */
    const std::initializer_list<Row> cusp = {{0, 0, 1}, {0, 1, 1}, {1, 1, -0.4142135623730951}};
    CHECK(Loops(cusp, {.Distance = 0.5, .Sides = Side::Outward}) == Counts{6});
    /**
     *  Within a tolerance as fine as the numbers allow the hairs the one-sided offset leaves are
     *  longer than what two ends are one within, and lead to no loop: they are left out.
     */
    const Polyline hairy = {
        .Vertices = {{.Bulge = 1},
                     {.Position = {.Y = 2}, .Bulge = -1},
                     {.Position = {.X = 1, .Y = 1}, .Bulge = -0.4142135623730951}},
        .Closed   = true};
    const std::vector<Polyline> around = Buffer(hairy, {.Distance = 0.5, .Sides = Side::Outward},
                                                {.Length = 0x1p-50, .Radians = 1e-10})
                                             .value();
    REQUIRE(around.size() == 1);
    CHECK(around.front().Vertices.size() == 5);
  }

  TEST_CASE("buffer: a circle, which no cut is made on, is kept whole on the side asked")
  {
    const std::initializer_list<Row> circle = {{0, 0, 1}, {10, 0, 1}};
    CHECK(Loops(circle, {.Distance = 1, .Sides = Side::Inward})  == Counts{2});
    CHECK(Loops(circle, {.Distance = 1, .Sides = Side::Outward}) == Counts{2});
    CHECK(Loops(circle, {.Distance = 1, .Sides = Side::Both})    == Counts{2, 2});
    /**
     *  Drawn round twice, it touches itself all along and is offset both ways: of the two whole
     *  circles each way gives, the one on the wrong side of what it fills is left out.
     */
    const std::initializer_list<Row> twice = {{0, 0, 1}, {10, 0, 1}, {0, 0, 1}, {10, 0, 1}};
    CHECK(Loops(twice, {.Distance = 1, .Sides = Side::Inward})  == Counts{2});
    CHECK(Loops(twice, {.Distance = 1, .Sides = Side::Outward}) == Counts{2});
    CHECK(Loops(twice, {.Distance = 1, .Sides = Side::Both})    == Counts{2, 2});
    /**
     *  Two circles that touch, drawn in one stroke: round the first, then round the second the
     *  other way. The raw offset to one side is the first grown and the second shrunk, touching and
     *  never crossing: yet half of it lies in what the circles fill and half does not, and each
     *  half is told by itself.
     */
    const std::initializer_list<Row> kissing = {{5, 0, 1}, {-5, 0, 1}, {5, 0, -1}, {15, 0, -1}};
    CHECK(Loops(kissing, {.Distance = 1, .Sides = Side::Inward})  == Counts{2, 2});
    CHECK(Loops(kissing, {.Distance = 1, .Sides = Side::Outward}) == Counts{2});
    CHECK(Loops(kissing, {.Distance = 1, .Sides = Side::Both})    == Counts{2, 2, 2});
    /**
     *  A closed line there and back fills nothing, and both its sides give the one loop about it.
     */
    const std::initializer_list<Row> back = {{0, 0, 0}, {10, 0, 0}};
    CHECK(Loops(back, {.Distance = 1, .Sides = Side::Inward}).empty());
    CHECK(Loops(back, {.Distance = 1, .Sides = Side::Outward}) == Counts{4});
    CHECK(Loops(back, {.Distance = 1, .Sides = Side::Both})    == Counts{4});
    /**
     *  A crescent between a half circle and two half its size, each met where the last turns
     *  straight back: what is moved inward from one comes nearer another with no crossing there.
     */
    const std::initializer_list<Row> crescent = {{0, 0, 1}, {0, 1, 1}, {0, 2, -1}};
    CHECK(Loops(crescent, {.Distance = 0.5, .Sides = Side::Inward}).empty());
    CHECK(Loops(crescent, {.Distance = 0.5, .Sides = Side::Outward}) == Counts{4});
    CHECK(Loops(crescent, {.Distance = 0.5, .Sides = Side::Both})    == Counts{4});
    const double                     quarter = std::tan(Pi / 8.0);
    const std::initializer_list<Row> leaf    = {{0, 0, 0}, {0, 2, quarter}, {1, 1, quarter}};
    CHECK(Loops(leaf, {.Distance = 1, .Sides = Side::Inward}).empty());
    /**
     *  Two squares that share a corner, drawn in one stroke. Inward the raw offset is a square in
     *  each and a corner between them that nothing cuts, in what neither fills.
     */
    const std::initializer_list<Row> corner = {{0, 0, 0},   {5, 0, 0},  {5, 5, 0}, {10, 5, 0},
                                               {10, 10, 0}, {5, 10, 0}, {5, 5, 0}, {0, 5, 0}};
    CHECK(Loops(corner, {.Distance = 1, .Sides = Side::Inward})  == Counts{4, 4});
    CHECK(Loops(corner, {.Distance = 1, .Sides = Side::Outward}) == Counts{14});
    CHECK(Loops(corner, {.Distance = 1, .Sides = Side::Both})    == Counts{4, 4, 14});
  }

  TEST_CASE("buffer: where pieces are joined, several may leave one point, and some lead "
            "to no loop")
  {
    const double quarter = 0.4142135623730951;
    /**
     *  A half circle and another back to a point beside the start: where the band's outline touches
     *  itself two edges leave one point, and each loop takes its own, the outline and the hole it
     *  closes there drawn apart.
     */
    const Polyline              touching = {.Vertices = {{.Bulge = 1},
                                                         {.Position = {.X = 1, .Y = 1}, .Bulge = 1},
                                                         {.Position = {.Y = 1}}},
                                            .Closed   = false};
    const std::vector<Polyline> outline  = Buffer(touching, {.Distance = 0.5, .Sides = Side::Both},
                                                  {.Length = 1e-15, .Radians = 1e-10})
                                               .value();
    CHECK(outline.size() == 2);
    CHECK(std::ranges::all_of(outline, &Polyline::Closed));
    /**
     *  Within a tolerance of a twentieth of the shape, ends far apart are one point and an edge can
     *  lead where every edge onward is taken: it is left out, and what comes back is closed loops
     *  still.
     */
    const Polyline coarse = {.Vertices = {{.Bulge = 1},
                                          {.Position = {.Y = 1}, .Bulge = 1},
                                          {.Position = {.X = 1, .Y = 1}, .Bulge = -quarter}},
                             .Closed   = true};
    const std::vector<Polyline> loose =
        Buffer(coarse, {.Distance = 1, .Sides = Side::Both}, {.Length = 0.05, .Radians = 1e-10})
            .value();
    CHECK(!loose.empty());
    CHECK(std::ranges::all_of(loose, &Polyline::Closed));
  }

  TEST_CASE("buffer: an open polyline's band has round ends")
  {
    const auto band = [](std::initializer_list<Row> rows, double distance) -> Counts
    {
      Polyline polyline = {.Closed = false};
      for (const Row &row : rows)
      {
        polyline.Vertices.push_back({.Position = {.X = row[0], .Y = row[1]}, .Bulge = row[2]});
      }
      Counts counts;
      for (const Polyline &loop :
           Buffer(polyline, {.Distance = distance, .Sides = Side::Both}, OffsetFine).value())
      {
        CHECK(loop.Closed);
        counts.push_back(loop.Vertices.size());
      }
      std::ranges::sort(counts);
      return counts;
    };
    /**
     *  A line: its two sides and a half circle round each end. An L: a corner's arc besides. A
     *  hairpin whose legs are nearer than twice the distance: one outline round both.
     */
    CHECK(band({{0, 0, 0}, {10, 0, 0}}, 1) == Counts{4});
    CHECK(band({{0, 0, 0}, {1, 0, 0}}, 1)  == Counts{4});
    CHECK(band({{0, 0, 0}, {10, 0, 0}, {10, 10, 0}}, 1) == Counts{7});
    CHECK(band({{0, 0, 0}, {10, 0, 0}, {10, 1.5, 0}, {0, 1.5, 0}}, 1) == Counts{7});
    /**
     *  Two quarter arcs out and back to the start: both ends are one point, and the one arc of
     *  their circle the band's edge takes is reached from both.
     */
    const double quarter = 0.4142135623730951;
    CHECK(band({{0, 0, quarter}, {0, 1, quarter}, {0, 0, 0}}, 0.5) == Counts{4});
    /**
     *  Two half circles whose offsets end on the end circles' half circles, touching them: the half
     *  circles are cut there all the same.
     */
    CHECK(band({{0, 0, -1}, {0, 1, 1}, {1, 1, 0}}, 0.5) == Counts{5});
    /**
     *  A half circle and a quarter arc back: the one-sided offset leaves a hair at the start of a
     *  piece as well as at an end, which nothing arrives at.
     */
    CHECK(band({{0, 0, 1}, {0, 2, -quarter}, {1, 1, 0}}, 0.5) == Counts{5});
    /**
     *  A line out from a point and an arc of radius 2 three quarters round it, offset by 1: the
     *  arc's inner side runs along the start's circle, the band on both sides of it there.
     */
    CHECK(band({{0, 0, 0}, {2, 0, 2.414213562373095}, {0, -2, 0}}, 1) == Counts{4});
    const Polyline line = {.Vertices = {{}, {.Position = {.X = 4}}}};
    CHECK(Area(Buffer(line, {.Distance = 1, .Sides = Side::Both}, OffsetFine).value().front())
          == doctest::Approx(8 + Pi));
    CHECK(Buffer(line, {.Distance = 1, .Sides = Side::Inward}, OffsetFine).error()
          == Error::Degenerate);
    CHECK(Buffer(line, {.Distance = 1, .Sides = Side::Outward}, OffsetFine).error()
          == Error::Degenerate);
  }

  TEST_CASE("buffer: outer corners are extended or chamfered where they point the way "
            "asked")
  {
    const auto with = [](Side sides, Gap gaps, double distance) -> Buffering
    {
      return {.Distance = distance, .Sides = sides, .Gaps = gaps};
    };
    const std::initializer_list<Row> square = {{0, 0, 0}, {10, 0, 0}, {10, 10, 0}, {0, 10, 0}};
    CHECK(Loops(square, with(Side::Outward, Gap::Extend, 1))  == Counts{4});
    CHECK(Loops(square, with(Side::Both, Gap::Extend, 1))     == Counts{4, 4});
    CHECK(Loops(square, with(Side::Outward, Gap::Chamfer, 1)) == Counts{8});
    CHECK(Loops(square, with(Side::Inward, Gap::Chamfer, 1))  == Counts{4});
    /**
     *  An L whose arms are twice the distance wide: inward its inner corner's mitre takes the last
     *  of it, and a chamfer leaves the corner the line across cuts off.
     */
    const std::initializer_list<Row> bent = {{0, 0, 0}, {10, 0, 0}, {10, 3, 0},
                                             {3, 3, 0}, {3, 10, 0}, {0, 10, 0}};
    CHECK(Loops(bent, with(Side::Inward, Gap::Extend, 1.5)).empty());
    CHECK(Loops(bent, with(Side::Outward, Gap::Extend, 1.5)) == Counts{6});
    CHECK(Loops(bent, with(Side::Inward, Gap::Chamfer, 1.5)) == Counts{3});
    CHECK(Loops(bent, with(Side::Both, Gap::Chamfer, 1.5))   == Counts{3, 11});
    /**
     *  Two squares that share a corner: the corners of each that point into the other stay round,
     *  outward, and inward none points into what they fill.
     */
    const std::initializer_list<Row> corner = {{0, 0, 0},   {5, 0, 0},  {5, 5, 0}, {10, 5, 0},
                                               {10, 10, 0}, {5, 10, 0}, {5, 5, 0}, {0, 5, 0}};
    CHECK(Loops(corner, with(Side::Outward, Gap::Extend, 1))  == Counts{8});
    CHECK(Loops(corner, with(Side::Inward, Gap::Extend, 1))   == Counts{4, 4});
    CHECK(Loops(corner, with(Side::Both, Gap::Extend, 1))     == Counts{4, 4, 8});
    CHECK(Loops(corner, with(Side::Outward, Gap::Chamfer, 1)) == Counts{14});
    CHECK(Loops(corner, {.Distance = 1, .Fill = Fills::EvenOdd, .Gaps = Gap::Extend}) == Counts{8});
    /**
     *  A corner that turns right back, a line into a quarter circle and that into a line: it is
     *  outer on the side its pieces bend away from alone, and a line there and back on neither.
     */
    const double                     quarter = 0x1.a827999fcef33p-2;
    const std::initializer_list<Row> cusped  = {{0, 0, 0}, {0, 1, 0}, {1, 1, quarter}};
    CHECK(Loops(cusped, with(Side::Outward, Gap::Extend, 0.5))  == Counts{7});
    CHECK(Loops(cusped, with(Side::Outward, Gap::Chamfer, 0.5)) == Counts{10});
    CHECK(Loops(cusped, with(Side::Both, Gap::Chamfer, 0.5))    == Counts{10});
    CHECK(Loops(cusped, with(Side::Inward, Gap::Extend, 0.5)).empty());
    /**
     *  A closed polyline has no ends, however they are asked; and an open one's ends are round
     *  beside corners extended: two strips less the square they share, the corner's square, and a
     *  half disk past each end.
     */
    CHECK(
        Loops(square, {.Distance = 1, .Sides = Side::Both, .Gaps = Gap::Extend, .Ends = Cap::Flat})
        == Counts{4, 4});
    const Polyline              elbow = {.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                                      {.Position = {.X = 10, .Y = 0}},
                                                      {.Position = {.X = 10, .Y = 10}}},
                                         .Closed   = false};
    const std::vector<Polyline> banded =
        Buffer(elbow, {.Distance = 1, .Sides = Side::Both, .Gaps = Gap::Extend}, OffsetFine)
            .value();
    REQUIRE(banded.size() == 1);
    CHECK(std::abs(Area(banded[0]) - (40.0 + Pi)) < 1e-9);
    const std::initializer_list<Row> back = {{0, 0, 0}, {10, 0, 0}};
    CHECK(Loops(back, with(Side::Outward, Gap::Chamfer, 1)) == Counts{4});
    CHECK(Loops(back, with(Side::Both, Gap::Extend, 1))     == Counts{4});
    /**
     *  A small shape far from the origin, simple, whose outward offset loses a stretch to a seam
     *  and is left in open pieces: they are joined into its one loop like any others.
     */
    const Polyline far = {
        .Vertices =
            {{.Position = {.X = 0x1.7d7840009374cp+26, .Y = -0x1.7d783fffdf3b6p+26}},
             {.Position = {.X = 0x1.7d783fff7ced9p+26, .Y = -0x1.7d784000624ddp+26}, .Bulge = -1.0},
             {.Position = {.X = 0x1.7d783fff8d4fep+26, .Y = -0x1.7d78400083127p+26}, .Bulge = -0.5},
             {.Position = {.X = 0x1.7d783fff6c8b4p+26, .Y = -0x1.7d78400051eb8p+26}, .Bulge = 0.5},
             {.Position = {.X = 0x1.7d783fffdf3b6p+26, .Y = -0x1.7d783fff9db23p+26}}},
        .Closed = true};
    const std::vector<Polyline> grown = Buffer(far, with(Side::Outward, Gap::Chamfer, 0.00025),
                                               {.Length = 0x1p-22, .Radians = 1e-10})
                                            .value();
    REQUIRE(grown.size() == 1);
    CHECK(grown[0].Closed);
  }

  TEST_CASE("buffer: an open polyline's band ends flat, or in a square past each end")
  {
    /**
     *  The area of an open polyline's band, one loop.
     */
    const auto area = [](std::initializer_list<Row> rows, double distance, Cap ends) -> double
    {
      Polyline polyline;
      for (const Row &row : rows)
      {
        polyline.Vertices.push_back({.Position = {.X = row[0], .Y = row[1]}, .Bulge = row[2]});
      }
      const std::vector<Polyline> loops =
          Buffer(polyline, {.Distance = distance, .Sides = Side::Both, .Ends = ends}, OffsetFine)
              .value();
      REQUIRE(loops.size() == 1);
      CHECK(loops[0].Closed);
      return Area(loops[0]);
    };
    const std::initializer_list<Row> line = {{0, 0, 0}, {10, 0, 0}};
    CHECK(std::abs(area(line, 1, Cap::Flat) - 20.0) < 1e-9);
    CHECK(std::abs(area(line, 1, Cap::Square) - 24.0) < 1e-9);
    /**
     *  Two strips, less the square they share, with the quarter disk about the corner.
     */
    const std::initializer_list<Row> bent  = {{0, 0, 0}, {10, 0, 0}, {10, 10, 0}};
    const double                     elbow = 39.0 + (Pi / 4.0);
    CHECK(std::abs(area(bent, 1, Cap::Flat) - elbow) < 1e-9);
    CHECK(std::abs(area(bent, 1, Cap::Square) - (elbow + 4.0)) < 1e-9);
    /**
     *  A vertex on the way changes nothing; and a line run back to its middle ends there in a
     *  square the band already holds, round at the far end it turns back at.
     */
    const std::initializer_list<Row> halved = {{0, 0, 0}, {0, 1, 0}, {0, 2, 0}};
    CHECK(std::abs(area(halved, 1, Cap::Flat) - 4.0) < 1e-9);
    CHECK(std::abs(area(halved, 1, Cap::Square) - 8.0) < 1e-9);
    const std::initializer_list<Row> back = {{0, 0, 0}, {0, 2, 0}, {0, 1, 0}};
    CHECK(std::abs(area(back, 1, Cap::Square) - (6.0 + (Pi / 2.0))) < 1e-9);
    /**
     *  A quarter circle tighter than the distance at the end: the line across it stops at its
     *  center, and the band's edge takes the circle about the vertex before it, past where the
     *  corner's fillet ends.
     */
    const double                     quarter = std::tan(Pi / 8.0);
    const std::initializer_list<Row> tight   = {{0, 0, 0}, {0, 1, -quarter}, {0, 2, 0}};
    CHECK(std::abs(area(tight, 1, Cap::Flat) - 4.418774416982109) < 1e-9);
    CHECK(std::abs(area(tight, 1, Cap::Square) - 8.194532219380358) < 1e-9);
    /**
     *  A half circle of the distance's radius into a line: the square past its start lies against
     *  the line's side of the band, and the side they share is no edge.
     */
    const std::initializer_list<Row> mirror = {{0, 0, 0}, {0, 1, quarter}, {0, 2, 0}};
    CHECK(std::abs(area(mirror, 1, Cap::Flat) - 4.418774416982109) < 1e-9);
    /**
     *  A half circle wider than the distance, alone: half a ring, flat across both ends.
     */
    const std::initializer_list<Row> bow = {{0, 0, 1}, {10, 0, 0}};
    CHECK(std::abs(area(bow, 1, Cap::Flat) - (10.0 * Pi)) < 1e-9);
    CHECK(std::abs(area(bow, 1, Cap::Square) - ((10.0 * Pi) + 4.0)) < 1e-9);
    /**
     *  A closed polyline has no ends, however they are asked.
     */
    const std::initializer_list<Row> square = {{0, 0, 0}, {10, 0, 0}, {10, 10, 0}, {0, 10, 0}};
    CHECK(Loops(square, {.Distance = 1, .Sides = Side::Both, .Ends = Cap::Flat}) == Counts{4, 8});
    const std::initializer_list<Row> hooked = {{0, 0, -1}, {0, 1, 0}, {1, 1, 0}};
    CHECK(std::abs(area(hooked, 0.5, Cap::Square) - (2.0 + (Pi / 2.0))) < 1e-9);
  }

  TEST_CASE("buffer: shapes whose curves touch, lie along each other or come the distance "
            "from another part give the area the definition gives")
  {
    /**
     *  A shape, what is asked of it, and the answer worked out to 120 digits: the area its loops
     *  enclose, a hole's taken off. The loops are not counted against it: where two touch they may
     *  be drawn as one or as two, and a hole too small to see may be left out.
     */
    struct Asked
    {
      bool             Closed   = true;
      Side             Sides    = Side::Both;
      Fills            Fill     = Fills::NonZero;
      Gap              Gaps     = Gap::Fillet;
      Cap              Ends     = Cap::Round;
      double           Distance = 0.0;
      double           Area     = 0.0;
      std::vector<Row> Rows     = {};
    };

    const std::vector<Asked> cases = {
        {.Closed   = false,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Square,
         .Distance = 0x1.592a31b5dd997p-17,
         .Area     = 4.905594506253391e-07,
         .Rows =
             {{-0x1.265207891da16p-7, -0x1.b02b4ea68cb76p-8, 0x1.a827999fcef33p-2},
              {0x1.1efc236084f0ap-7, 0x1.65eec7e01382dp-8,
               0x1.6b19671b06100p-10}}}, // corpus/100900/1.028671327295005e-05/Both/NonZero/Chamfer/Square
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Fillet,
         .Ends     = Cap::Round,
         .Distance = 0.25,
         .Area     = 14.160490540377934,
         .Rows     = {{4.0, 6.0, -1.0},
                      {-4.0, 4.0, -0.5},
                      {4.0, -2.0, 0.0}}}, // corpus/101909/0.25/Both/NonZero/Round/Round
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Extend,
         .Ends     = Cap::Round,
         .Distance = 0.75,
         .Area     = 77.0955810546875,
         .Rows     = {{100000010.0, -0x1.7d78400000000p+26, 1.0},
                      {99999990.0, -0x1.7d78400000000p+26,
                       0.0}}}, // corpus/108972/0.75/Both/NonZero/Mitre/Round
        {.Closed   = true,
         .Sides    = Side::Inward,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0.75,
         .Area     = 120.54147839546204,
         .Rows     = {{100000010.0, -0x1.7d78400000000p+26, 1.0},
                      {99999990.0, -0x1.7d78400000000p+26,
                       0.0}}}, // corpus/108972/0.75/Inward/NonZero/Chamfer/Round
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Extend,
         .Ends     = Cap::Round,
         .Distance = 0x1.139612f6389a0p+1,
         .Area     = 154.03250317433003,
         .Rows =
             {{0x1.a1cb798816080p-1, 0x1.1b53d2c6ce7c4p+3, 1.0},
              {-0x1.98c624c04387ep+2, 0x1.d4af4c3a96bd0p+2, 1.0},
              {0x1.28e2b93e23d80p-1, -0x1.3417ae37349f0p+0,
               -0x1.dba32251ab376p+0}}}, // corpus/122089/2.1530173971415962/Both/NonZero/Mitre/Round
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0x1.6a86c2bd4eeecp-7,
         .Area     = 0.0013619002260799253,
         .Rows =
             {{0.01, 0.0, -0x1.a827999fcef33p-2},
              {0x1.f608ee4e066afp-8, 0x1.a541dd9e682b3p-8, 0.0},
              {0x1.c73551c379e24p-10, 0x1.42b3a975cdec7p-7, 0.0},
              {-0x1.47ae147ae1479p-8, 0x1.1bc779ef0d4e9p-7, 0.0},
              {-0x1.33eb215f7271cp-7, 0x1.c04aea9a675b8p-9, 0.0},
              {-0x1.33eb215f7271dp-7, -0x1.c04aea9a675b3p-9, -0x1.15fc79ce98b92p+1},
              {-0x1.47ae147ae1480p-8, -0x1.1bc779ef0d4e7p-7, -0x1.a3941d7aa0ff7p-4},
              {0x1.c73551c379e10p-10, -0x1.42b3a975cdec8p-7, 0.0},
              {0x1.f608ee4e066adp-8, -0x1.a541dd9e682b7p-8,
               -0x1.02680c321653cp-1}}}, // corpus/132179/0.01106342800868939/Both/NonZero/Chamfer/Round
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Extend,
         .Ends     = Cap::Round,
         .Distance = 2.0,
         .Area     = 612.2564045124368,
         .Rows     = {{100000007.0, -99999965.0, 0.5},
                      {99999972.0, -0x1.7d78454000000p+26,
                       0.5}}}, // corpus/204827/2.0/Both/NonZero/Mitre/Round
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Fillet,
         .Ends     = Cap::Round,
         .Distance = 0x1.2c65a0b2d3221p+4,
         .Area     = 18350338.83386236,
         .Rows     = {{-0x1.db9c4683f2d23p+12, 0x1.68a935b615346p+12, 0.0},
                      {-0x1.2ebb0cbff7a7fp+13, 0x1.11672da0dd74bp+13, 0.0},
                      {-0x1.1a0666a10ce84p+10, -0x1.93d3ecebe0486p+11, -1.0},
                      {0x1.566d974bf5f48p+12, 0x1.2ee170897b5f5p+12, -0x1.62c9d86986c80p-6},
                      {-0x1.7645753f11512p+10, -0x1.635fc9b0c42ebp+12, 0.0},
                      {0x1.2714a5a67f490p+11, -0x1.21217548ee893p+13, 0.0},
                      {0x1.0c0d833e7d880p+10, -0x1.67903e3c1f272p+12, 0x1.96affb8897188p-5},
                      {0x1.24286f1018779p+10, 0x1.289f5447cb01dp+13, -1.0},
                      {0x1.2cda08b0009c0p+13, 0x1.2ac3b70063d50p+13, -1.0},
                      {0x1.94dedfd9832a3p+12, -0x1.723dc0f63d941p+12, -0x1.414fbc470fa6fp+0},
                      {0x1.b755e7e8ad4b9p+12, -0x1.31ec3fbf536b9p+13, 0.0},
                      {-0x1.590ca0a1f07b1p+11, 0x1.bcfcf3d008431p+12, 0.0},
                      {0x1.8c5ddf457810ep+9, 0x1.f16de09c3d678p+12, 0x1.a827999fcef33p-2},
                      {-0x1.aaba3d6e01fd4p+12, -0x1.4b53222f02239p+11, 1.0},
                      {0x1.6a75dfe6988bfp+11, 0x1.9f8a67820ff2bp+11, -0x1.d6fc283f0e284p+0},
                      {0x1.a95a2afe3536dp+11, -0x1.79c91a1ed49d9p+12, 0.0},
                      {-0x1.852afedc06835p+11, -0x1.0c90db76ba186p+13, -0x1.28164409d85e0p-3},
                      {0x1.689b09562126ep+12, 0x1.e856f7e179cf9p+12, 0x1.a7ca04b1fb9b8p-1},
                      {-0x1.73e97a5ceba0cp+12, -0x1.2bfd22f04d1bfp+13, -0x1.8c0691e6d4000p-17},
                      {-0x1.7e57b8f769dedp+12, 0x1.9d00df90ea84fp+12, 1.0},
                      {0x1.3477b66f36356p+13, -0x1.2c979c5953fc0p+13, -0x1.357f11b8f6084p-1},
                      {-0x1.342a5bf263ed3p+13, 0x1.614f15783c4a1p+12, 0x1.a827999fcef33p-2},
                      {-0x1.a3297937e25b2p+12, 0x1.444fb49457cdep+12, 1.0},
                      {0x1.2bbe47ed15295p+12, -0x1.29ae98258d491p+11, 0x1.0c85139918cd4p+1},
                      {0x1.145d2ec2213d1p+11, 0x1.01af2d6f1060ep+13, 1.0},
                      {-0x1.4b5e6b6d87a18p+4, 0x1.0b549f3498e17p+12, 0x1.65f6bacd972b8p+0},
                      {0x1.67709c7af3aefp+12, 0x1.8d089ef5851f5p+12, 0.0},
                      {-0x1.5a1cbeb9e701dp+12, -0x1.19289fe83e760p+12, 1.0},
                      {0x1.99cf48f8ab8c0p+12, -0x1.03cfd7ca1c50fp+12, -0x1.112b2bc58548ep+1},
                      {-0x1.f1bfc9e1a8f88p+11, -0x1.3188abcfae22fp+13,
                       0.0}}}, // corpus/78420/18.77481145718969/Both/NonZero/Round/Round
        {.Closed   = false,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Fillet,
         .Ends     = Cap::Flat,
         .Distance = 1.0,
         .Area     = 59.57079632679489,
         .Rows     = {{0.0, 0.0, 0.0},
                      {10.0, 0.0, 0.0},
                      {10.0, 10.0, 0.0},
                      {20.0, 10.0, 0.0}}}, // known/zigzag/1/Both/NonZero/Round/Butt
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Fillet,
         .Ends     = Cap::Round,
         .Distance = 0x1.8000000000001p+0,
         .Area     = 117.10287520586556,
         .Rows     = {{0.0, 0.0, 0.0},
                      {10.0, 0.0, 0.0},
                      {10.0, 3.0, 0.0},
                      {3.0, 3.0, 0.0},
                      {3.0, 10.0, 0.0},
                      {0.0, 10.0, 0.0}}}, // near/L/1.5000000000000002/Both/NonZero/Round/Round
        {.Closed   = false,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Fillet,
         .Ends     = Cap::Round,
         .Distance = 0x1.4000000000001p+2,
         .Area     = 235.61944901923448,
         .Rows = {{0.0, 0.0, 1.0},
                  {10.0, 0.0, 0.0}}}, // near/half-circle/5.000000000000001/Both/NonZero/Round/Round
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0x1.fffffffffffffp+0,
         .Area     = 130.2654824574367,
         .Rows = {{0.0, 0.0, 0.0},
                  {10.0, 0.0, 1.0},
                  {10.0, 4.0, 0.0},
                  {0.0, 4.0, 1.0}}}, // near/stadium/1.9999999999999998/Both/NonZero/Chamfer/Round
        {.Closed   = false,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Fillet,
         .Ends     = Cap::Square,
         .Distance = 1.0,
         .Area     = 7.308662255316469,
         .Rows     = {{0.0, 0.0, 1.0},
                      {0.0, 1.0, 1.0},
                      {0.0, 0.0, 0.0}}}, // small/00+_01+_000o/1.0/Both/NonZero/Round/Square
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 3.9450754879692242,
         .Rows     = {{0.0, 0.0, 1.0},
                      {0.0, 1.0, 1.0},
                      {1.0, 0.0, -1.0}}}, // small/00+_01+_10-c/0.5/Both/NonZero/Chamfer/Round
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Extend,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 3.2534819882443373,
         .Rows     = {{0.0, 0.0, 1.0},
                      {0.0, 1.0, 1.0},
                      {1.0, 1.0,
                       -0x1.a827999fcef33p-2}}}, // small/00+_01+_11Qc/0.5/Both/NonZero/Mitre/Round
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 3.945075487969224,
         .Rows     = {{0.0, 0.0, 1.0},
                      {0.0, 1.0, -1.0},
                      {1.0, 1.0, 1.0}}}, // small/00+_01-_11+c/0.5/Both/NonZero/Chamfer/Round
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Fillet,
         .Ends     = Cap::Round,
         .Distance = 1.0,
         .Area     = 18.14490970987465,
         .Rows     = {{0.0, 0.0, 1.0},
                      {0.0, 1.0, -1.0},
                      {2.0, 2.0, -1.0}}}, // small/00+_01-_22-c/1.0/Both/NonZero/Round/Round
        {.Closed   = false,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Fillet,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 6.093182583916945,
         .Rows     = {{0.0, 0.0, 1.0},
                      {0.0, 2.0, 1.0},
                      {1.0, 0.0, 0.0}}}, // small/00+_02+_100o/0.5/Both/NonZero/Round/Round
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 7.405799484818728,
         .Rows = {{0.0, 0.0, 1.0},
                  {0.0, 2.0, 1.0},
                  {2.0, 0.0,
                   -0x1.a827999fcef33p-2}}}, // small/00+_02+_20Qc/0.5/Both/NonZero/Chamfer/Round
        {.Closed   = false,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Fillet,
         .Ends     = Cap::Round,
         .Distance = 1.0,
         .Area     = 10.566090624703588,
         .Rows     = {{0.0, 0.0, 1.0},
                      {0.0, 2.0, -1.0},
                      {1.0, 1.0, 0.0}}}, // small/00+_02-_110o/1.0/Both/NonZero/Round/Round
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Fillet,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 4.690426395626385,
         .Rows     = {{0.0, 0.0, -1.0},
                      {0.0, 1.0, -1.0},
                      {1.0, 1.0,
                       -0x1.a827999fcef33p-2}}}, // small/00-_01-_11Qc/0.5/Both/NonZero/Round/Round
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Fillet,
         .Ends     = Cap::Round,
         .Distance = 1.0,
         .Area     = 7.433830449463873,
         .Rows     = {{0.0, 0.0, 0.0},
                      {0.0, 1.0, 1.0},
                      {1.0, 0.0, 0.0}}}, // small/000_01+_100c/1.0/Both/NonZero/Round/Round
        {.Closed   = false,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Fillet,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 3.446953020859589,
         .Rows     = {{0.0, 0.0, -0x1.a827999fcef33p-2},
                      {0.0, 1.0, -0x1.a827999fcef33p-2},
                      {1.0, 0.0, 0.0}}}, // small/00Q_01Q_100o/0.5/Both/NonZero/Round/Round
        {.Closed   = false,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Fillet,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 3.9269908169872414,
         .Rows     = {{0.0, 1.0, 1.0},
                      {0.0, 0.0, 0x1.a827999fcef33p-2},
                      {1.0, 1.0, 0.0}}}, // small/01+_00q_110o/0.5/Both/NonZero/Round/Round
        {.Closed   = false,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Fillet,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 3.4469530208595893,
         .Rows     = {{0.0, 1.0, 0x1.a827999fcef33p-2},
                      {0.0, 0.0, 0x1.a827999fcef33p-2},
                      {1.0, 1.0, 0.0}}}, // small/01q_00q_110o/0.5/Both/NonZero/Round/Round
        {.Closed   = false,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Square,
         .Distance = 0x1.fffffffffe000p-1,
         .Area     = 60.46189071954072,
         .Rows = {{0.0, 0.0, 1.0},
                  {4.0, 0.0, 1.0},
                  {-2.0, 0.0, 1.0},
                  {6.0, 0.0, 0.0}}}, // near/spiral/0.9999999999990905/Both/NonZero/Chamfer/Square
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Extend,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 5.497787143782138,
         .Rows     = {{0.0, 0.0, 1.0},
                      {0.0, 1.0, -1.0},
                      {0.0, 2.0, 1.0}}}, // small/00+_01-_02+c/0.5/Both/NonZero/Mitre/Round
        {.Closed   = false,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Flat,
         .Distance = 0x1.b9ba29b0b0f32p+5,
         .Area     = 33316.902008370605,
         .Rows =
             {{70.0, 0.0, 0x1.a827999fcef33p-2},
              {0x1.34db9c7de44a8p-48, 70.0, 0.0},
              {-70.0, 0x1.34db9c7de44a8p-47, -0x1.581ddf4191e10p-3},
              {-0x1.cf496abcd66fap-47, -70.0,
               0x1.90f9a7c089fb4p-1}}}, // corpus/206845/55.215899830246045/Both/NonZero/Chamfer/Butt
        {.Closed   = true,
         .Sides    = Side::Inward,
         .Fill     = Fills::EvenOdd,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 0,
         .Rows     = {{0.0, 0.0, 1.0},
                      {0.0, 1.0, 1.0},
                      {1.0, 0.0, 1.0}}}, // small/00+_01+_10+c/0.5/Inward/EvenOdd/Chamfer/Round
        {.Closed   = true,
         .Sides    = Side::Inward,
         .Fill     = Fills::EvenOdd,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 0,
         .Rows = {{0.0, 0.0, 1.0},
                  {0.0, 1.0, 1.0},
                  {1.0, 2.0,
                   0x1.a827999fcef33p-2}}}, // small/00+_01+_12qc/0.5/Inward/EvenOdd/Chamfer/Round
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Extend,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 4.108870492746559,
         .Rows     = {{0.0, 0.0, 1.0},
                      {0.0, 1.0, -1.0},
                      {1.0, 0.0,
                       0x1.a827999fcef33p-2}}}, // small/00+_01-_10qc/0.5/Both/NonZero/Mitre/Round
        {.Closed   = true,
         .Sides    = Side::Inward,
         .Fill     = Fills::EvenOdd,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 0,
         .Rows     = {{0.0, 0.0, 1.0},
                      {1.0, 1.0, 1.0},
                      {2.0, 2.0, 0.0}}}, // small/00+_11+_220c/0.5/Inward/EvenOdd/Chamfer/Round
        {.Closed   = false,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Extend,
         .Ends     = Cap::Flat,
         .Distance = 0.5,
         .Area     = 2.8869927500604446,
         .Rows     = {{0.0, 0.0, 1.0},
                      {1.0, 1.0, -0x1.a827999fcef33p-2},
                      {0.0, 0.0, 0.0}}}, // small/00+_11Q_000o/0.5/Both/NonZero/Mitre/Butt
        {.Closed   = true,
         .Sides    = Side::Outward,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Extend,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 43.52817459346924,
         .Rows     = {{0.0, 0.0, 0.0},
                      {10.0, 0.0, 0.0},
                      {10.0, 5.0, 0.0},
                      {5.0, 1e-10, 0.0},
                      {0.0, 5.0, 0.0}}}, // a notch whose tip comes 1e-10 from the far side
        {.Closed   = true,
         .Sides    = Side::Inward,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Extend,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 10.843145750836905,
         .Rows     = {{0.0, 0.0, 0.0},
                      {10.0, 0.0, 0.0},
                      {10.0, 5.0, 0.0},
                      {5.0, 1e-10, 0.0},
                      {0.0, 5.0, 0.0}}},
        {.Closed   = true,
         .Sides    = Side::Outward,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Extend,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 42.00295482139265,
         .Rows     = {{0.0, 0.0, 0.0},
                      {5.0, 0.0, -0.2},
                      {10.0, 0.0, 0.0},
                      {10.0, 5.0, 0.0},
                      {0x1.3fffffffe4832p+2, 1e-10, 0.0},
                      {0.0, 5.0, 0.0}}}, // the same beside where a line runs into an arc
        {.Closed   = true,
         .Sides    = Side::Inward,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Extend,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 9.527618512385924,
         .Rows     = {{0.0, 0.0, 0.0},
                      {5.0, 0.0, -0.2},
                      {10.0, 0.0, 0.0},
                      {10.0, 5.0, 0.0},
                      {0x1.3fffffffe4832p+2, 1e-10, 0.0},
                      {0.0, 5.0, 0.0}}},
        {.Closed   = true,
         .Sides    = Side::Inward,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 2.0,
         .Area     = 4884.725474388321,
         .Rows     = {{0.0, 28.0, 0.0},
                      {-42.0, -42.0, 1.0},
                      {42.0, -14.0, 0.0}}}, // corpus/172539/2.0/Inward/NonZero/Chamfer/Round
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Extend,
         .Ends     = Cap::Round,
         .Distance = 0x1.13d60e5d5cadep+13,
         .Area     = 648359656.7835279,
         .Rows =
             {{0x1.420d2168cd688p+10, -0x1.0d80dfca398dep+10, 0.0},
              {0x1.984594801a44cp+10, 0x1.0c43c0b56ea3dp+9, 0.0},
              {-0x1.e647075891463p+12, -0x1.52c41a063cb75p+11,
               -0x1.d455cf041956ap+0}}}, // corpus/226583/8826.757014011437/Both/NonZero/Mitre/Round
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Extend,
         .Ends     = Cap::Round,
         .Distance = 0x1.13d60e5d5cadep+13,
         .Area     = 648359656.7835279,
         .Rows     = {{0x1.420d2168cd688p+10, -0x1.0d80dfca398dep+10, 0x1.d455cf041956ap+0},
                      {-0x1.e647075891463p+12, -0x1.52c41a063cb75p+11, -0.0},
                      {0x1.984594801a44cp+10, 0x1.0c43c0b56ea3dp+9,
                       -0.0}}}, // the same run the other way
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Fillet,
         .Ends     = Cap::Round,
         .Distance = 0.25,
         .Area     = 14.160490540377934,
         .Rows     = {{4.0, 6.0, -1.0},
                      {-4.0, 4.0, -0.5},
                      {4.0, -2.0, 0.0}}}, // corpus/101909/0.25/Both/NonZero/Round/Round
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0x1.2c65a0b2d3221p+4,
         .Area     = 18351490.428638887,
         .Rows     = {{-0x1.db9c4683f2d23p+12, 0x1.68a935b615346p+12, 0.0},
                      {-0x1.2ebb0cbff7a7fp+13, 0x1.11672da0dd74bp+13, 0.0},
                      {-0x1.1a0666a10ce84p+10, -0x1.93d3ecebe0486p+11, -1.0},
                      {0x1.566d974bf5f48p+12, 0x1.2ee170897b5f5p+12, -0x1.62c9d86986c80p-6},
                      {-0x1.7645753f11512p+10, -0x1.635fc9b0c42ebp+12, 0.0},
                      {0x1.2714a5a67f490p+11, -0x1.21217548ee893p+13, 0.0},
                      {0x1.0c0d833e7d880p+10, -0x1.67903e3c1f272p+12, 0x1.96affb8897188p-5},
                      {0x1.24286f1018779p+10, 0x1.289f5447cb01dp+13, -1.0},
                      {0x1.2cda08b0009c0p+13, 0x1.2ac3b70063d50p+13, -1.0},
                      {0x1.94dedfd9832a3p+12, -0x1.723dc0f63d941p+12, -0x1.414fbc470fa6fp+0},
                      {0x1.b755e7e8ad4b9p+12, -0x1.31ec3fbf536b9p+13, 0.0},
                      {-0x1.590ca0a1f07b1p+11, 0x1.bcfcf3d008431p+12, 0.0},
                      {0x1.8c5ddf457810ep+9, 0x1.f16de09c3d678p+12, 0x1.a827999fcef33p-2},
                      {-0x1.aaba3d6e01fd4p+12, -0x1.4b53222f02239p+11, 1.0},
                      {0x1.6a75dfe6988bfp+11, 0x1.9f8a67820ff2bp+11, -0x1.d6fc283f0e284p+0},
                      {0x1.a95a2afe3536dp+11, -0x1.79c91a1ed49d9p+12, 0.0},
                      {-0x1.852afedc06835p+11, -0x1.0c90db76ba186p+13, -0x1.28164409d85e0p-3},
                      {0x1.689b09562126ep+12, 0x1.e856f7e179cf9p+12, 0x1.a7ca04b1fb9b8p-1},
                      {-0x1.73e97a5ceba0cp+12, -0x1.2bfd22f04d1bfp+13, -0x1.8c0691e6d4000p-17},
                      {-0x1.7e57b8f769dedp+12, 0x1.9d00df90ea84fp+12, 1.0},
                      {0x1.3477b66f36356p+13, -0x1.2c979c5953fc0p+13, -0x1.357f11b8f6084p-1},
                      {-0x1.342a5bf263ed3p+13, 0x1.614f15783c4a1p+12, 0x1.a827999fcef33p-2},
                      {-0x1.a3297937e25b2p+12, 0x1.444fb49457cdep+12, 1.0},
                      {0x1.2bbe47ed15295p+12, -0x1.29ae98258d491p+11, 0x1.0c85139918cd4p+1},
                      {0x1.145d2ec2213d1p+11, 0x1.01af2d6f1060ep+13, 1.0},
                      {-0x1.4b5e6b6d87a18p+4, 0x1.0b549f3498e17p+12, 0x1.65f6bacd972b8p+0},
                      {0x1.67709c7af3aefp+12, 0x1.8d089ef5851f5p+12, 0.0},
                      {-0x1.5a1cbeb9e701dp+12, -0x1.19289fe83e760p+12, 1.0},
                      {0x1.99cf48f8ab8c0p+12, -0x1.03cfd7ca1c50fp+12, -0x1.112b2bc58548ep+1},
                      {-0x1.f1bfc9e1a8f88p+11, -0x1.3188abcfae22fp+13,
                       0.0}}}, // corpus/78420/18.77481145718969/Both/NonZero/Chamfer/Round
        {.Closed   = false,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Square,
         .Distance = 0x1.b652357add754p+3,
         .Area     = 1171497.6370326905,
         .Rows =
             {{10000.0, 0.0, -0x1.596e6525c6fa4p-3},
              {0x1.58b51a15a5a57p-41, 10000.0, 0x1.1bf542792f734p-4},
              {-10000.0, 0x1.58b51a15a5a57p-40, 0.0},
              {-0x1.0287d3903c3c0p-39, -10000.0,
               0x1.042ab70482efep+1}}}, // corpus/81729/13.697535266848412/Both/NonZero/Chamfer/Square
        {.Closed   = false,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Square,
         .Distance = 0x1.eef97b61d4b5fp+3,
         .Area     = 2106.124354659055,
         .Rows =
             {{10.0, 0.0, 0.0},
              {0x1.c48c600000000p+2, 0x1.c48c600000000p+2, 0x1.2bf461ef83bb0p-2},
              {0.0, 10.0, 0x1.a827999fcef33p-2},
              {-0x1.c48c600000000p+2, 0x1.c48c600000000p+2, 0.0},
              {-10.0, 0.0, 0.0},
              {-0x1.c48c600000000p+2, -0x1.c48c600000000p+2, -0x1.2f7a043849f4fp-3},
              {0.0, -10.0, -1.0},
              {0x1.c48c600000000p+2, -0x1.c48c600000000p+2,
               1.0}}}, // corpus/90810/15.467954341007212/Both/NonZero/Chamfer/Square, moved to 0
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0x1.3ffffffffec00p+2,
         .Area     = 314.1592653586936,
         .Rows = {{-5.0, 0.0, 1.0},
                  {5.0, 0.0, 1.0}}}, // near/circle/4.9999999999954525/Both/NonZero/Chamfer/Round
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0x1.fffffffffe000p-1,
         .Area     = 57.34124422294579,
         .Rows     = {{0.0, 0.0, 0.0},
                      {10.0, 0.0, 0.0},
                      {10.0, 4.0, 0.0},
                      {6.0, 4.0, -1.0},
                      {4.0, 4.0, 0.0},
                      {0.0, 4.0, 0.0}}}, // near/dent/0.9999999999990905/Both/NonZero/Chamfer/Round
        {.Closed   = false,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0x1.4000000000001p+2,
         .Area     = 235.61944901923448,
         .Rows     = {{0.0, 0.0, 1.0},
                      {10.0, 0.0,
                       0.0}}}, // near/half-circle/5.000000000000001/Both/NonZero/Chamfer/Round
        {.Closed   = true,
         .Sides    = Side::Inward,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 0,
         .Rows     = {{0.0, 0.0, 1.0},
                      {0.0, 1.0, -1.0},
                      {0.0, 2.0, 1.0}}}, // small/00+_01-_02+c/0.5/Inward/NonZero/Chamfer/Round
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 5.161614153387452,
         .Rows = {{0.0, 0.0, 1.0},
                  {0.0, 2.0, 1.0},
                  {1.0, 1.0,
                   -0x1.a827999fcef33p-2}}}, // small/00+_02+_11Qc/0.5/Both/NonZero/Chamfer/Round
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 4.940633624641554,
         .Rows = {{0.0, 0.0, 1.0},
                  {0.0, 2.0, -1.0},
                  {1.0, 1.0,
                   -0x1.a827999fcef33p-2}}}, // small/00+_02-_11Qc/0.5/Both/NonZero/Chamfer/Round
        {.Closed   = false,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Flat,
         .Distance = 0.5,
         .Area     = 3.6811130264695913,
         .Rows     = {{0.0, 0.0, 1.0},
                      {0.0, 2.0, -0x1.a827999fcef33p-2},
                      {1.0, 1.0, 0.0}}}, // small/00+_02Q_110o/0.5/Both/NonZero/Chamfer/Butt
        {.Closed   = true,
         .Sides    = Side::Inward,
         .Fill     = Fills::EvenOdd,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 0,
         .Rows = {{0.0, 0.0, 1.0},
                  {0.0, 2.0, -0x1.a827999fcef33p-2},
                  {1.0, 1.0,
                   -0x1.a827999fcef33p-2}}}, // small/00+_02Q_11Qc/0.5/Inward/EvenOdd/Chamfer/Round
        {.Closed   = true,
         .Sides    = Side::Inward,
         .Fill     = Fills::EvenOdd,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 0,
         .Rows     = {{0.0, 0.0, 0.0},
                      {0.0, 1.0, 1.0},
                      {0.0, 2.0, 0.0}}}, // small/000_01+_020c/0.5/Inward/EvenOdd/Chamfer/Round
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 3.122494661309621,
         .Rows     = {{0.0, 0.0, 0.0},
                      {0.0, 1.0, -0x1.a827999fcef33p-2},
                      {1.0, 1.0,
                       0x1.a827999fcef33p-2}}}, // small/000_01Q_11qc/0.5/Both/NonZero/Chamfer/Round
        {.Closed   = true,
         .Sides    = Side::Outward,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Extend,
         .Ends     = Cap::Round,
         .Distance = 0x1.3d921c3334320p+2,
         .Area     = 157.30142025744615,
         .Rows = {{0x1.15cbba77d34c6p+1, 0x1.15cf93b94c277p+2, 0x1.a827999fcef33p-2},
                  {0x1.15ce321549a71p+2, 0x1.15d0d2b5ce260p+1, 0x1.a827999fcef33p-2},
                  {0x1.a0b6ba7d8ac1cp+2, 0x1.15cc250c99f36p+2, 0.0},
                  {0x1.15d02a4a5b702p+2, 0x1.a0b8c232686e1p+2,
                   0x1.a827999fcef33p-2}}}, // wfc/10/4.962042856220052/Outward/NonZero/Mitre/Round
        {.Closed   = true,
         .Sides    = Side::Inward,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Extend,
         .Ends     = Cap::Round,
         .Distance = 0x1.37b85a406e373p+1,
         .Area     = 1.901750865108422,
         .Rows     = {{0x1.bcd6d0ff35206p-1, 0x1.a109640816508p+3, 0.0},
                      {0x1.1606429d698acp+2, 0x1.a10963f5b5c99p+3, 0.0},
                      {0x1.160642c013601p+2, 0x1.d8a43e305cc68p+3, 0.0},
                      {0x1.f471ab1760bdep+2, 0x1.d8a43e252a596p+3, 0.0},
                      {0x1.f471ab58fbba4p+2, 0x1.696e89cb00ff9p+3, 0.0},
                      {0x1.31d3afc72f030p+3, 0x1.696e89d7d9513p+3, 0.0},
                      {0x1.31d3afb2f2f00p+3, 0x1.31d3afc855e87p+3, 0.0},
                      {0x1.696e89eacb51bp+3, 0x1.31d3afc77bf33p+3, 0.0},
                      {0x1.696e89ca99b12p+3, 0x1.081f8c2348435p+4, 0.0},
                      {0x1.31d3afb6ea43ap+3, 0x1.081f8c1721495p+4, 0.0},
                      {0x1.31d3afbd92500p+3, 0x1.d8a43e0cb2692p+3, 0.0},
                      {0x1.d8a43e1af161dp+3, 0x1.d8a43e2c94cd3p+3, 0.0},
                      {0x1.d8a43e1519f3ep+3, 0x1.a10963f4c6deap+3, 0.0},
                      {0x1.a10963fc38d62p+3, 0x1.a10964048fa99p+3, 0.0},
                      {0x1.a109640b34d46p+3, 0x1.081f8c1c6c935p+4, 0.0},
                      {0x1.081f8c29ebfd9p+4, 0x1.081f8c1f6add0p+4, 0.0},
                      {0x1.081f8c286e14cp+4, 0x1.f471ab0aa05c7p+2, 0.0},
                      {0x1.d8a43e22f2619p+3, 0x1.f471ab0606921p+2, 0.0},
                      {0x1.d8a43e1c413c0p+3, 0x1.31d3afb83e5cdp+3, 0.0},
                      {0x1.a10963fb11e49p+3, 0x1.31d3afc76af85p+3, 0.0},
                      {0x1.a10963f36adbbp+3, 0x1.f471ab0eae8a0p+2, 0.0},
                      {0x1.31d3afaa13416p+3, 0x1.f471ab1464e6fp+2, 0.0},
                      {0x1.31d3afb31d5a6p+3, 0x1.853bf6e682046p+2, 0.0},
                      {0x1.f471ab19e5513p+2, 0x1.853bf6f2c3a7cp+2, 0.0},
                      {0x1.f471ab1c532c1p+2, 0x1.160642adcb368p+2, 0.0},
                      {0x1.4da11d121229bp+1, 0x1.160642b7d7943p+2, 0.0},
                      {0x1.4da11cb676693p+1, 0x1.d8a43e336e45ep+3, 0.0},
                      {0x1.bcd6d14eec262p-1, 0x1.d8a43e27e5301p+3,
                       0.0}}}, // wfc/47/2.435313493213863/Inward/NonZero/Mitre/Round
        {.Closed   = false,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Square,
         .Distance = 0.125,
         .Area     = 2.342534795552452,
         .Rows     = {{0x1.0000248b209bdp-1, 0x1.c00000019bb54p+1, 0.0},
                      {0x1.0000055b34119p-1, 0x1.1ffff9e87f34bp+2, 0.0},
                      {0x1.80000956f0465p+1, 0x1.1ffff8d51cc3ep+2, -0x1.a827999fcef33p-2},
                      {0x1.bffffb5012bc1p+1, 0x1.0000076c469dap+2, -0x1.a827999fcef33p-2},
                      {0x1.8000065135c6ep+1, 0x1.c00007b1a0dbfp+1, 0.0},
                      {0x1.4000044c7f790p+1, 0x1.bfffff2bf0e0fp+1, 0.0},
                      {0x1.3ffff2a13020bp+1, 0x1.0000060070aa7p+1, -0x1.a827999fcef33p-2},
                      {0x1.0000035b73b65p+1, 0x1.8000183b25addp+0, -0x1.a827999fcef33p-2},
                      {0x1.7ffff7d43a614p+0, 0x1.00000ed7137b0p+1, 0.0},
                      {0x1.7fffefb781387p+0, 0x1.40000cabfa56dp+1,
                       0.0}}}, // wfc/63/0.125/Both/NonZero/Chamfer/Square
        {.Closed   = false,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Extend,
         .Ends     = Cap::Square,
         .Distance = 0x1.21732b272bd58p-14,
         .Area     = 6.208771749788403e-08,
         .Rows     = {{-0x1.323c74bbdd10cp-12, 0x1.268a6178fe212p-13, 0.0},
                      {-0x1.a64b6916cfb40p-12, 0x1.111628cdbb404p-13, 0.0},
                      {-0x1.b10587b688e75p-12, 0x1.f93400090d29cp-13, 0.0},
                      {-0x1.3cf6952bcc809p-12, 0x1.0754180df4b73p-12,
                       0.0}}}, // wfc/68/6.901022844617766e-05/Both/NonZero/Mitre/Square
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0x1.87e56722ed2c8p-2,
         .Area     = 29.34587865594737,
         .Rows     = {{0x1.fffffc0f161f8p-2, 0x1.ffffff7188908p+2, 0x1.a827999fcef33p-2},
                      {0x1.0000019a7df63p+0, 0x1.dffffffc5b831p+2, 0x1.a827999fcef33p-2},
                      {0x1.7ffffc6dbeba3p+0, 0x1.ffffff1c921a2p+2, -0x1.a827999fcef33p-2},
                      {0x1.000001c03e760p+1, 0x1.0fffffea55b51p+3, -0x1.a827999fcef33p-2},
                      {0x1.3fffffa28734ep+1, 0x1.ffffff97b25b1p+2, 0.0},
                      {0x1.3fffff3d36498p+1, 0x1.bfffffa00c73ep+2, 0x1.a827999fcef33p-2},
                      {0x1.8000000ffc737p+1, 0x1.9fffff7c9ebf7p+2, 0x1.a827999fcef33p-2},
                      {0x1.bfffff9f10de3p+1, 0x1.c00000b43c444p+2, -0x1.a827999fcef33p-2},
                      {0x1.000000652cb24p+2, 0x1.dfffff6fbfb3dp+2, 0.0},
                      {0x1.7fffff7960266p+2, 0x1.dfffff7a24ce9p+2, -0x1.a827999fcef33p-2},
                      {0x1.9fffff9b08653p+2, 0x1.bfffff8dd9ac2p+2, 0.0},
                      {0x1.a00000082ab2ap+2, 0x1.7fffff3997734p+2, 0x1.a827999fcef33p-2},
                      {0x1.bfffffdf2652dp+2, 0x1.5ffffffafed13p+2, 0.0},
                      {0x1.ffffffdb572b0p+2, 0x1.5fffffd24a59dp+2, 0x1.a827999fcef33p-2},
                      {0x1.0fffff9c67d2ap+3, 0x1.7fffff2bd1d55p+2, 0.0},
                      {0x1.0fffffa016faep+3, 0x1.a000008b32fbcp+2, 0.0},
                      {0x1.dfffff71283e5p+2, 0x1.a00000754d887p+2, 0.0},
                      {0x1.e00000b3a35fcp+2, 0x1.000000079cc32p+3, -0x1.a827999fcef33p-2},
                      {0x1.ffffff5193991p+2, 0x1.0ffffff6a19bbp+3, -0x1.a827999fcef33p-2},
                      {0x1.1000000d0e12dp+3, 0x1.0000006e13909p+3, 0x1.a827999fcef33p-2},
                      {0x1.2000001512557p+3, 0x1.e000008e6c761p+2, -0x1.a827999fcef33p-2},
                      {0x1.300000709056fp+3, 0x1.c00000dc0e835p+2, 0x1.a827999fcef33p-2},
                      {0x1.3ffffff613f35p+3, 0x1.a00000469e93fp+2, 0.0},
                      {0x1.5000005209d6dp+3, 0x1.a0000046d064bp+2, 0.0},
                      {0x1.5000006619449p+3, 0x1.00000054659e4p+3, 0x1.a827999fcef33p-2},
                      {0x1.3ffffff2b29e8p+3, 0x1.0fffffa80425fp+3, 0.0},
                      {0x1.2fffff8e04cf4p+3, 0x1.0fffff9f66ed2p+3, 0.0},
                      {0x1.2fffffda58ec4p+3, 0x1.3000000a18eabp+3, 0.0},
                      {0x1.0000000649be2p+3, 0x1.3000000c352dfp+3, -0x1.a827999fcef33p-2},
                      {0x1.dfffff55038d7p+2, 0x1.3fffffde5e4e3p+3, 0x1.a827999fcef33p-2},
                      {0x1.c00000e5b0f6ap+2, 0x1.5000006342991p+3, 0x1.a827999fcef33p-2},
                      {0x1.9fffffb438833p+2, 0x1.3ffffff54e7c0p+3, 0.0},
                      {0x1.a00000c17a2b0p+2, 0x1.2fffffb0cd0cap+3, 0.0},
                      {0x1.5fffff918cf08p+2, 0x1.300000565a52cp+3, 0.0},
                      {0x1.5fffffcbaab70p+2, 0x1.500000674e7ebp+3, 0.0},
                      {0x1.3fffff8efcad7p+2, 0x1.500000384abf6p+3, 0x1.a827999fcef33p-2},
                      {0x1.2000005816420p+2, 0x1.4000005fdd0b3p+3, -0x1.a827999fcef33p-2},
                      {0x1.0000000fb5314p+2, 0x1.2fffffd42eec2p+3, 0.0},
                      {0x1.8000005982f3dp+1, 0x1.3000000f836f5p+3, -0x1.a827999fcef33p-2},
                      {0x1.3fffffcd9bb56p+1, 0x1.4000000b999a6p+3, 0x1.a827999fcef33p-2},
                      {0x1.000000b577971p+1, 0x1.5000004d9b0c8p+3, 0x1.a827999fcef33p-2},
                      {0x1.800000918aa28p+0, 0x1.3fffffffb205fp+3, 0.0},
                      {0x1.7fffff2fe4f3ep+0, 0x1.30000051c44bfp+3, 0.0},
                      {0x1.fffffddceb7c4p-1, 0x1.300000425b46fp+3, 0x1.a827999fcef33p-2},
                      {0x1.00000078efd89p-1, 0x1.1fffffefba748p+3,
                       0.0}}}, // wfc/62/0.38271104002641865/Both/NonZero/Chamfer/Round
        {.Closed   = false,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Extend,
         .Ends     = Cap::Square,
         .Distance = 0x1.9b43f6a19c3eep+15,
         .Area     = 28070838233.975655,
         .Rows     = {{0x1.c7ab0d9bff51dp+14, 0x1.c7ab49304e6ffp+14, 0.0},
                      {0x1.c7ab3639b9dffp+14, 0x1.c7ab28015b349p+16, -0x1.a827999fcef33p-2},
                      {0x1.c7ab1ffdda973p+15, 0x1.1ccafe579d801p+17, 0.0},
                      {0x1.55c0622d7ded4p+16, 0x1.1ccafd6974314p+17,
                       0.0}}}, // wfc/101/52641.981701739/Both/NonZero/Mitre/Square
        {.Closed   = false,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0x1.c7ab2ec466f2fp+14,
         .Area     = 12147399730.96443,
         .Rows     = {{0x1.c7ab0d9bff51dp+14, 0x1.c7ab49304e6ffp+14, 0.0},
                      {0x1.c7ab3639b9dffp+14, 0x1.c7ab28015b349p+16, -0x1.a827999fcef33p-2},
                      {0x1.c7ab1ffdda973p+15, 0x1.1ccafe579d801p+17, 0.0},
                      {0x1.55c0622d7ded4p+16, 0x1.1ccafd6974314p+17,
                       0.0}}}, // wfc/101/29162.795671089723/Both/NonZero/Chamfer/Round
        {.Closed   = false,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Flat,
         .Distance = 0x1.26b800f0cf163p-10,
         .Area     = 5.064129662919661e-06,
         .Rows     = {{0x1.8b174a41f96a8p-4, -0x1.ffc2df3f1d821p-5, 0.0},
                      {0x1.888da16e95d6dp-4, -0x1.f796c8cd8652dp-5, 0.0},
                      {0x1.84779635ca3f3p-4, -0x1.fcaa1a744d7a3p-5,
                       0.0}}}, // wfc/186/0.00112426286450703/Both/NonZero/Chamfer/Butt
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Extend,
         .Ends     = Cap::Round,
         .Distance = 1.0,
         .Area     = 15.0,
         .Rows     = {{0.5, 0.5, 0.0},
                      {1.5, 0.5, 0.0},
                      {1.5, 1.5, 0.0},
                      {2.5, 1.5, 0.0},
                      {2.5, 2.5, 0.0},
                      {0.5, 2.5, 0.0}}}, // wfc/193/1.0/Both/NonZero/Mitre/Round
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0.5,
         .Area     = 5.313301260054891,
         .Rows     = {{-0x1.1f45bcd040011p+10, 0x1.8f072c9bdfc0cp+8, 0x1.a827999fcef33p-2},
                      {-0x1.1f25bcd040bfdp+10, 0x1.8e872c9bce7d0p+8, 0.0},
                      {-0x1.1ee5bcd042619p+10, 0x1.8e872c9bc4093p+8, 0x1.a827999fcef33p-2},
                      {-0x1.1ec5bcd046b02p+10, 0x1.8f072c9bd8b48p+8, 0.0},
                      {-0x1.1ec5bcd03f856p+10, 0x1.8f872c9bcf40cp+8, 0.0},
                      {-0x1.1f25bcd047b63p+10, 0x1.8f872c9bb39bcp+8,
                       0x1.a827999fcef33p-2}}}, // wfc/209/0.5/Both/NonZero/Chamfer/Round
        {.Closed   = true,
         .Sides    = Side::Outward,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Fillet,
         .Ends     = Cap::Round,
         .Distance = 0x1.aa1092374f1fep-1,
         .Area     = 5.575205798438948,
         .Rows = {{0x1.ffffcd742634dp-2, 0x1.ffffc5da70779p-1, 0x1.a827999fcef33p-2},
                  {0x1.ffffd03501867p-1, 0x1.00001c1781db6p-1, 0x1.a827999fcef33p-2},
                  {0x1.7fffec37ca2d4p+0, 0x1.0000085dfda62p+0, 0x1.a827999fcef33p-2},
                  {0x1.0000001150e78p+0, 0x1.7ffffd1f17e80p+0,
                   0x1.a827999fcef33p-2}}}, // wfc/215/0.832157677890393/Outward/NonZero/Round/Round
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Extend,
         .Ends     = Cap::Round,
         .Distance = 0x1.8c399e912cf2ap-5,
         .Area     = 0.30390060897140425,
         .Rows = {{0x1.ffffcd742634dp-2, 0x1.ffffc5da70779p-1, 0x1.a827999fcef33p-2},
                  {0x1.ffffd03501867p-1, 0x1.00001c1781db6p-1, 0x1.a827999fcef33p-2},
                  {0x1.7fffec37ca2d4p+0, 0x1.0000085dfda62p+0, 0x1.a827999fcef33p-2},
                  {0x1.0000001150e78p+0, 0x1.7ffffd1f17e80p+0,
                   0x1.a827999fcef33p-2}}}, // wfc/215/0.0483673188222425/Both/NonZero/Mitre/Round
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Extend,
         .Ends     = Cap::Round,
         .Distance = 0x1.fe5329a3e43dep+17,
         .Area     = 1405722440050.4888,
         .Rows     = {{-0x1.107553e281213p+18, -0x1.e1a53352ee8ecp+14, 0x1.a827999fcef33p-2},
                      {-0x1.d80e838bfc2c3p+17, 0x1.16c912ec9ce65p+17, 0.0},
                      {-0x1.f628d6d698957p+18, 0x1.139f33614078bp+19, 0x1.a827999fcef33p-2},
                      {-0x1.4fd3d9b5c2da4p+19, 0x1.25d63c720e958p+19, 0.0},
                      {-0x1.0ee0641eb9129p+20, 0x1.418ae410a04c5p+18, 0x1.a827999fcef33p-2},
                      {-0x1.17fbe8a36e4c9p+20, 0x1.30180ec3592ddp+17, 0.0},
                      {-0x1.ed6f6c088e8acp+19, -0x1.af0739e452348p+15, 0x1.a827999fcef33p-2},
                      {-0x1.98affdcdc425ap+19, -0x1.693be63fbc37ap+16, 0.0},
                      {-0x1.31b98695d9dbap+19, 0x1.560e87c7b8fbbp+15, 0x1.a827999fcef33p-2},
                      {-0x1.1f827d62b10b6p+19, 0x1.a8815b85df598p+17, 0x1.a827999fcef33p-2},
                      {-0x1.7441ebad9ec19p+19, 0x1.f15d7f6bc1cd5p+17, 0x1.a827999fcef33p-2},
                      {-0x1.8678f4d69e7e6p+19, 0x1.3cbf8c718773fp+16, 0.0},
                      {-0x1.43f08f796b08fp+19, -0x1.faf42e603b81ep+16, 0x1.a827999fcef33p-2},
                      {-0x1.de6242691763cp+18, -0x1.46563b1da823cp+17,
                       0.0}}}, // wfc/212/261286.32531407377/Both/NonZero/Mitre/Round
        {.Closed   = true,
         .Sides    = Side::Outward,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Fillet,
         .Ends     = Cap::Round,
         .Distance = 0x1.ea5aa1eacd860p+14,
         .Area     = 576321796621.3948,
         .Rows     = {{-0x1.107553e281213p+18, -0x1.e1a53352ee8ecp+14, 0x1.a827999fcef33p-2},
                      {-0x1.d80e838bfc2c3p+17, 0x1.16c912ec9ce65p+17, 0.0},
                      {-0x1.f628d6d698957p+18, 0x1.139f33614078bp+19, 0x1.a827999fcef33p-2},
                      {-0x1.4fd3d9b5c2da4p+19, 0x1.25d63c720e958p+19, 0.0},
                      {-0x1.0ee0641eb9129p+20, 0x1.418ae410a04c5p+18, 0x1.a827999fcef33p-2},
                      {-0x1.17fbe8a36e4c9p+20, 0x1.30180ec3592ddp+17, 0.0},
                      {-0x1.ed6f6c088e8acp+19, -0x1.af0739e452348p+15, 0x1.a827999fcef33p-2},
                      {-0x1.98affdcdc425ap+19, -0x1.693be63fbc37ap+16, 0.0},
                      {-0x1.31b98695d9dbap+19, 0x1.560e87c7b8fbbp+15, 0x1.a827999fcef33p-2},
                      {-0x1.1f827d62b10b6p+19, 0x1.a8815b85df598p+17, 0x1.a827999fcef33p-2},
                      {-0x1.7441ebad9ec19p+19, 0x1.f15d7f6bc1cd5p+17, 0x1.a827999fcef33p-2},
                      {-0x1.8678f4d69e7e6p+19, 0x1.3cbf8c718773fp+16, 0.0},
                      {-0x1.43f08f796b08fp+19, -0x1.faf42e603b81ep+16, 0x1.a827999fcef33p-2},
                      {-0x1.de6242691763cp+18, -0x1.46563b1da823cp+17,
                       0.0}}}, // wfc/212/31382.658122264198/Outward/NonZero/Round/Round
        {.Closed   = true,
         .Sides    = Side::Inward,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Extend,
         .Ends     = Cap::Round,
         .Distance = 0x1.e571516c75dd0p-8,
         .Area     = 0.05413488353558422,
         .Rows     = {{0x1.a3b60addef73dp-2, -0x1.4a535d3aaef36p-3, 0x1.a827999fcef33p-2},
                      {0x1.846b55c3e5218p-2, -0x1.850a7e0d5e51fp-3, -0x1.a827999fcef33p-2},
                      {0x1.6520a0a9dacf1p-2, -0x1.bfc19ee00db08p-3, 0x1.a827999fcef33p-2},
                      {0x1.45d5eb8fd07ccp-2, -0x1.fa78bfb2bd0f1p-3, 0.0},
                      {0x1.43e6c6df1dd9ap-2, -0x1.39e2a55cc0893p-2, -0x1.a827999fcef33p-2},
                      {0x1.249c11c513875p-2, -0x1.573e35c618388p-2, -0x1.a827999fcef33p-2},
                      {0x1.0740815bbbd7fp-2, -0x1.37f380ac0de62p-2, 0x1.a827999fcef33p-2},
                      {0x1.d3c9e1e4c8517p-3, -0x1.18a8cb920393bp-2, 0x1.a827999fcef33p-2},
                      {0x1.953477b0b3accp-3, -0x1.36045bfb5b431p-2, -0x1.a827999fcef33p-2},
                      {0x1.569f0d7c9f07fp-3, -0x1.535fec64b2f25p-2, 0x1.a827999fcef33p-2},
                      {0x1.1809a3488a634p-3, -0x1.70bb7cce0aa19p-2, 0x1.a827999fcef33p-2},
                      {0x1.52c0c41b39c1cp-3, -0x1.900631e814f3fp-2, -0x1.a827999fcef33p-2},
                      {0x1.8d77e4ede9205p-3, -0x1.af50e7021f465p-2, 0.0},
                      {0x1.89999b8c83da3p-3, -0x1.ebf72c8581480p-2, 0x1.a827999fcef33p-2},
                      {0x1.c450bc5f3338cp-3, -0x1.05a0f0cfc5cd2p-1, 0x1.a827999fcef33p-2},
                      {0x1.01731349a3eecp-2, -0x1.ede6513633eb0p-2, -0x1.a827999fcef33p-2},
                      {0x1.20bdc863ae412p-2, -0x1.d08ac0ccdc3bdp-2, 0x1.a827999fcef33p-2},
                      {0x1.40087d7db8938p-2, -0x1.b32f3063848c7p-2, -0x1.a827999fcef33p-2},
                      {0x1.5f533297c2e5dp-2, -0x1.95d39ffa2cdd3p-2, -0x1.a827999fcef33p-2},
                      {0x1.7caec3011a952p-2, -0x1.b51e5514372f9p-2, 0.0},
                      {0x1.7abf9e5067f21p-2, -0x1.f1c49a9799313p-2, 0x1.a827999fcef33p-2},
                      {0x1.981b2eb9bfa15p-2, -0x1.0887a7d8d1c1dp-1, 0x1.05d2bbe221515p-2},
                      {0x1.d4c1743d21a30p-2, -0x1.097f3a312b136p-1, 0x1.a827999fcef33p-2},
                      {0x1.f40c29572bf56p-2, -0x1.f5a2e3f8fe777p-2, 0x1.a827999fcef33p-2},
                      {0x1.d6b098edd4461p-2, -0x1.d6582edef4251p-2, -0x1.a827999fcef33p-2},
                      {0x1.b95508847c96dp-2, -0x1.b70d79c4e9d2ap-2, 0.0},
                      {0x1.bb442d352f39ep-2, -0x1.7a67344187d10p-2, -0x1.a827999fcef33p-2},
                      {0x1.da8ee24f398c4p-2, -0x1.5d0ba3d83021bp-2, 0x1.a827999fcef33p-2},
                      {0x1.f9d9976943deap-2, -0x1.3fb0136ed8727p-2, 0x1.a827999fcef33p-2},
                      {0x1.dc7e06ffec2f5p-2, -0x1.20655e54ce201p-2, -0x1.a827999fcef33p-2},
                      {0x1.bf22769694801p-2, -0x1.011aa93ac3cdcp-2, -0x1.a827999fcef33p-2},
                      {0x1.de6d2bb09ed26p-2, -0x1.c77e31a2d83cdp-3, 0x1.a827999fcef33p-2},
                      {0x1.fdb7e0caa924dp-2, -0x1.8cc710d028de4p-3, 0x1.a827999fcef33p-2},
                      {0x1.e05c506151757p-2, -0x1.4e31a69c14399p-3,
                       0.0}}}, // wfc/204/0.007407266978098417/Inward/NonZero/Mitre/Round
        {.Closed   = true,
         .Sides    = Side::Inward,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0x1.438c0c9743923p-11,
         .Area     = 1.3212314896470838e-06,
         .Rows =
             {{0x1.438c0d573ef69p-11, 0x1.438c0c0852573p-10, 0x1.a827999fcef33p-2},
              {0x1.438c0bd377a4dp-10, 0x1.438c0df83bf69p-11, 0.0},
              {0x1.438c0c8b2759ep-9, 0x1.438c0bbe20eb9p-11, 0x1.a827999fcef33p-2},
              {0x1.946f0fe43af2ap-9, 0x1.438c0bf5691f8p-10, 0.0},
              {0x1.946f100dac699p-9, 0x1.438c0c4295655p-9, 0x1.1225ebc2dcea4p-3},
              {0x1.438c0cd26b92bp-9, 0x1.946f0fc4eacdcp-9, -0x1.fe6405f508b9cp-3},
              {0x1.438c0bf88b5b3p-10, 0x1.946f0ffdef720p-9, -0x1.e0459c6c1f35ap-4},
              {0x1.438c0bde14c89p-11, 0x1.438c0c9772020p-9,
               -0x1.63c389b6a87aep-5}}}, // wfc/227/0.0006171170561246924/Inward/NonZero/Chamfer/Round
        {.Closed   = true,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0x1.08af109e75304p-2,
         .Area     = 0.442581335815995,
         .Rows = {{0x1.66f7d5514e475p-2, 0x1.dea705f6259f8p-2, 0x1.a827999fcef33p-2},
                  {0x1.de9b732e527fcp-2, 0x1.67002e73635a3p-2, 0x1.a827999fcef33p-2},
                  {0x1.2b298156e5c37p-1, 0x1.de9e82ded19eap-2, 0x1.a827999fcef33p-2},
                  {0x1.de9bc8c343fecp-2, 0x1.2b25295d8cd81p-1,
                   0x1.a827999fcef33p-2}}}, // wfc/336/0.2584803196635195/Both/NonZero/Chamfer/Round
        {.Closed   = false,
         .Sides    = Side::Both,
         .Fill     = Fills::NonZero,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Flat,
         .Distance = 0x1.ca48d8df3b554p+2,
         .Area     = 607.1649489488462,
         .Rows =
             {{0x1.c071096185862p+1, 0x1.2772689edf4a2p+0, 0x1.a827999fcef33p-2},
              {0x1.053ab49675137p+2, 0x1.392efba81ab3ep+2, -0x1.35db448efd76fp+1},
              {0x1.76a3e51519a73p-1, 0x1.5f9428f462aedp+2, -0x1.b3cb39c929670p-1},
              {-0x1.043191d5957b9p+0, 0x1.d158066797564p+1, 0.0},
              {-0x1.67b33d4f047f1p+2, 0x1.1c690ab974ec4p+2, 0.0},
              {-0x1.d20b1909aa100p+2, -0x1.4a7f49d4b483bp-2, 0.0},
              {-0x1.416d6691c16c0p+2, -0x1.4b6f61c1f58bap+1, -0x1.a827999fcef33p-2},
              {-0x1.abb1443653305p+1, -0x1.9496125ed9fc4p+2, 0x1.76891023f46d4p-1},
              {0x1.9f174be8ef7c5p-1, -0x1.e5b4ade831667p+2, 0.0},
              {0x1.80eb9143b5f05p+2, -0x1.32c9e15be9703p+2, -0x1.2116151514ecep-4},
              {0x1.04577fec5bde2p+2, -0x1.0743057bcea5fp+1,
               -0x1.36871decf77a0p-1}}}, // corpus/227025/7.160696237563872/Both/NonZero/Chamfer/Butt
        {.Closed   = true,
         .Sides    = Side::Inward,
         .Fill     = Fills::EvenOdd,
         .Gaps     = Gap::Chamfer,
         .Ends     = Cap::Round,
         .Distance = 0x1.139612f6389a0p+1,
         .Area     = 0,
         .Rows =
             {{0x1.a1cb798816080p-1, 0x1.1b53d2c6ce7c4p+3, 1.0},
              {-0x1.98c624c04387ep+2, 0x1.d4af4c3a96bd0p+2, 1.0},
              {0x1.28e2b93e23d80p-1, -0x1.3417ae37349f0p+0,
               -0x1.dba32251ab376p+0}}}, // corpus/122089/2.1530173971415962/Inward/EvenOdd/Chamfer/Round
    };
    for (const Asked &asked : cases)
    {
      CAPTURE(asked.Area);
      Polyline polyline = {.Closed = asked.Closed};
      for (const Row &row : asked.Rows)
      {
        polyline.Vertices.push_back({.Position = {.X = row[0], .Y = row[1]}, .Bulge = row[2]});
      }
      const std::vector<Polyline> loops = Buffer(polyline,
                                                 {.Distance = asked.Distance,
                                                  .Sides    = asked.Sides,
                                                  .Fill     = asked.Fill,
                                                  .Gaps     = asked.Gaps,
                                                  .Ends     = asked.Ends},
                                                 {.Length = 1e-12, .Radians = 1e-10})
                                              .value();
      double                      area  = 0.0;
      for (const Polyline &loop : loops)
      {
        CHECK(loop.Closed);
        area += Area(loop);
      }
      CHECK(std::abs(area - asked.Area) < 1e-8 * std::max(1.0, std::abs(asked.Area)));
    }
  }

  TEST_CASE("buffer: a shape of many segments is cut and measured through grids")
  {
    /**
     *  A gear of 64 vertices, every fourth segment an arc: more curves than are each tried against
     *  every other, more segments than are each measured, more points than are searched one by one.
     *  The areas are the definition's, worked out to 120 digits.
     */
    Polyline gear = {.Closed = true};
    for (int index = 0; index < 64; ++index)
    {
      const double angle  = index * 2 * Pi / 64;
      const double radius = index % 2 != 0 ? 10 : 8;
      gear.Vertices.push_back(
          {.Position = {.X = radius * std::cos(angle), .Y = radius * std::sin(angle)},
           .Bulge    = index % 4 == 0 ? 0.2 : 0.0});
    }
    const auto areas = [&gear](Side sides) -> std::vector<double>
    {
      std::vector<double> found;
      for (const Polyline &loop :
           Buffer(gear, {.Distance = 0.7, .Sides = sides}, OffsetFine).value())
      {
        CHECK(loop.Closed);
        found.push_back(Area(loop));
      }
      std::ranges::sort(found);
      return found;
    };
    const std::vector<double> outward = areas(Side::Outward);
    const std::vector<double> inward  = areas(Side::Inward);
    const std::vector<double> both    = areas(Side::Both);
    REQUIRE(outward.size() == 1);
    REQUIRE(inward.size()  == 1);
    REQUIRE(both.size()    == 2);
    /**
     *  A gear of more pieces than are looked at for whether a shape crosses itself: outward it is
     *  offset on both sides like any other, and gives what both sides give outside it.
     */
    Polyline fine = {.Closed = true};
    for (int index = 0; index < 200; ++index)
    {
      const double angle  = index * 2 * Pi / 200;
      const double radius = index % 2 != 0 ? 10 : 9.5;
      fine.Vertices.push_back(
          {.Position = {.X = radius * std::cos(angle), .Y = radius * std::sin(angle)}});
    }
    const std::vector<Polyline> grown =
        Buffer(fine, {.Distance = 0.3, .Sides = Side::Outward}, OffsetFine).value();
    const std::vector<Polyline> banded =
        Buffer(fine, {.Distance = 0.3, .Sides = Side::Both}, OffsetFine).value();
    REQUIRE(grown.size()  == 1);
    REQUIRE(banded.size() == 2);
    CHECK(std::abs(Area(grown[0]) - std::max(Area(banded[0]), Area(banded[1]))) < 1e-9);
    CHECK(std::abs(outward[0] - 338.0408350803273) < 1e-8);
    CHECK(std::abs(inward[0] - 176.8077853206636) < 1e-8);
    CHECK(std::abs(both[0] + 176.8077853206636) < 1e-8);
    CHECK(std::abs(both[1] - 338.0408350803273) < 1e-8);
  }

  TEST_CASE("buffer: what cannot be offset is reported")
  {
    const Buffering by     = {.Distance = 1};
    const Polyline  square = {.Vertices = {{},
                                           {.Position = {.X = 4}},
                                           {.Position = {.X = 4, .Y = 4}},
                                           {.Position = {.Y = 4}}},
                              .Closed   = true};
    Polyline        open   = square;
    open.Closed            = false;
    CHECK(Buffer(open, by, OffsetFine).error() == Error::Degenerate);
    CHECK(Buffer(Polyline{.Closed = true}, by, OffsetFine).error() == Error::Degenerate);
    CHECK(Buffer(square, {.Distance = 0}, OffsetFine).error()      == Error::Degenerate);
    CHECK(Buffer(square, {.Distance = -1}, OffsetFine).error()     == Error::Degenerate);
    CHECK(Buffer(square, {.Distance = std::numeric_limits<double>::quiet_NaN()}, OffsetFine).error()
          == Error::NotFinite);
  }

  TEST_CASE("buffer: a gap a million millionth wide between two edges is kept open")
  {
    /**
     *  Three half circles in a spiral, banded by a hair under 1: the first's outer edge and the
     *  third's inner one are circles about one center 1.8e-12 apart, and the gap between them opens
     *  through a notch 2.7e-6 tall beside the square end. Points within 1e-12 of the size were one,
     *  which closed the gap into a loop of no area and lost the notch; and the square's edge beside
     *  it was looked at 1e-7 to its right, past the gap. Each vertex is the truth's, to the bit.
     */
    const Polyline spiral = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 1},
                                          {.Position = {.X = 4, .Y = 0}, .Bulge = 1},
                                          {.Position = {.X = -2, .Y = 0}, .Bulge = 1},
                                          {.Position = {.X = 6, .Y = 0}}}};
    const double   near   = 0x1.fffffffffep-1;
    const auto     band   = Buffer(
        spiral, {.Distance = near, .Sides = Side::Both, .Gaps = Gap::Chamfer, .Ends = Cap::Square},
        OffsetFine);
    REQUIRE(band.has_value());
    REQUIRE(band->size() == 2);
    const double                under = 1 - near;
    const std::array<Vertex, 9> outer = {
        Vertex{.Position = {.X = -near, .Y = 0}, .Bulge = 1},
        Vertex{.Position = {.X = 4 + near, .Y = 0}, .Bulge = 1},
        Vertex{.Position = {.X = -2 - near, .Y = 0}, .Bulge = 1},
        Vertex{.Position = {.X = 6 + near, .Y = 0}, .Bulge = 0},
        Vertex{.Position = {.X = 6 + near, .Y = near}, .Bulge = 0},
        Vertex{.Position = {.X = 5 + under, .Y = near}, .Bulge = 0},
        Vertex{.Position = {.X = 5 + under, .Y = 0}, .Bulge = -1},
        Vertex{.Position = {.X = -1 - under, .Y = 0}, .Bulge = -3.3717478808711395e-07},
        Vertex{.Position = {.X = -near, .Y = 2.6973983046972182e-06}, .Bulge = 0}};
    const Polyline &loop = band->front();
    REQUIRE(loop.Vertices.size() == outer.size());
    for (std::size_t index = 0; index < outer.size(); ++index)
    {
      CHECK(loop.Vertices[index].Position == outer.at(index).Position);
      CHECK(loop.Vertices[index].Bulge    == doctest::Approx(outer.at(index).Bulge).epsilon(1e-12));
    }
    CHECK(band->back().Vertices.size() == 4);
  }

  TEST_CASE("buffer: a hole thinner at its tips than a million millionth keeps its shape")
  {
    /**
     *  A rounded shape near (-1148, 399) banded by half: its band leaves a hole 0.07 long and 4e-10
     *  wide, a sliver between two lines and an arc that meet at tips thinner than 1e-12 of the
     *  size. The tips' points were made one and the hole a dot; it has its four vertices, each the
     *  truth's to the bit.
     */
    const Polyline shape = {
        .Vertices = {{.Position = {.X = -0x1.1f45bcd040011p+10, .Y = 0x1.8f072c9bdfc0cp+8},
                      .Bulge    = 0x1.a827999fcef33p-2},
                     {.Position = {.X = -0x1.1f25bcd040bfdp+10, .Y = 0x1.8e872c9bce7dp+8}},
                     {.Position = {.X = -0x1.1ee5bcd042619p+10, .Y = 0x1.8e872c9bc4093p+8},
                      .Bulge    = 0x1.a827999fcef33p-2},
                     {.Position = {.X = -0x1.1ec5bcd046b02p+10, .Y = 0x1.8f072c9bd8b48p+8}},
                     {.Position = {.X = -0x1.1ec5bcd03f856p+10, .Y = 0x1.8f872c9bcf40cp+8}},
                     {.Position = {.X = -0x1.1f25bcd047b63p+10, .Y = 0x1.8f872c9bb39bcp+8},
                      .Bulge    = 0x1.a827999fcef33p-2}},
        .Closed   = true};
    const auto band =
        Buffer(shape, {.Distance = 0.5, .Sides = Side::Both, .Gaps = Gap::Chamfer}, OffsetFine);
    REQUIRE(band.has_value());
    REQUIRE(band->size() == 2);
    CHECK(band->front().Vertices.size() == 9);
    const Polyline &hole = band->back();
    REQUIRE(hole.Vertices.size()    == 4);
    CHECK(hole.Vertices[0].Position == Point{.X = -1147.5896492662728, .Y = 399.02802442148266});
    CHECK(hole.Vertices[1].Position == Point{.X = -1147.5896492666732, .Y = 399.02802442108185});
    CHECK(hole.Vertices[2].Position == Point{.X = -1147.6589724286418, .Y = 399.02802442125056});
    CHECK(hole.Vertices[3].Position == Point{.X = -1147.5896492662728, .Y = 399.02802442154803});
  }

  TEST_CASE("buffer: a sharp corner cut across keeps its two run-on edges apart")
  {
    /**
     *  A corner of 150 degrees, slanted, banded with its corners cut across: the two edges run on
     *  from the moved pieces to the line across are edges of one mitre and meet nowhere, though
     *  their boxes lie over each other. One loop, with the cut's two vertices on it.
     */
    const Polyline bent = {.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                        {.Position = {.X = 4, .Y = 3}},
                                        {.Position = {.X = 0.25, .Y = 0.75}}}};
    for (const Gap gaps : {Gap::Chamfer, Gap::Extend})
    {
      const auto band =
          Buffer(bent, {.Distance = 0.5, .Sides = Side::Both, .Gaps = gaps}, OffsetFine);
      REQUIRE(band.has_value());
      REQUIRE(band->size() == 1);
      CHECK(band->front().Closed);
    }
  }

  TEST_CASE("buffer: a line cut across a hairpin crosses the circle it all but touches")
  {
    /**
     *  A half circle and a quarter arc back along it, the arc's bulge tan(pi / 8) as a double holds
     *  it: the arc's circle lies 5.6e-17 off the half circle's, and the shape turns back on itself
     *  by a hair under half a turn. The line cut across the corner leans by half that, and crosses
     *  the inner circle 5.9e-9 above where a line square to the heading in touches it; the way out
     *  of the corner was taken for that heading, the sum of the two moved ends' ways coming to
     *  nothing there, and the two outer arcs, their circles one as doubles hold them, were made
     *  one. Every vertex is the truth's, to the bit.
     */
    const double   quarter = 0x1.a827999fcef33p-2;
    const Polyline open    = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 1},
                                           {.Position = {.X = 0, .Y = 2}, .Bulge = -quarter},
                                           {.Position = {.X = 1, .Y = 1}}}};
    const auto     band    = Buffer(
        open, {.Distance = 0.5, .Sides = Side::Both, .Gaps = Gap::Chamfer, .Ends = Cap::Flat},
        OffsetFine);
    REQUIRE(band.has_value());
    REQUIRE(band->size() == 1);
    const std::array<Point, 5> truth = {
        Point{.X =    0, .Y =              -0.5},
        Point{.X =  1.5, .Y =                 1},
        Point{.X = -0.5, .Y = 2.414213562373095},
        Point{.X = -0.49999999999999994, .Y = 1.0000000059273753}, Point{.X = 0, .Y = 0.5}};
    REQUIRE(band->front().Vertices.size() == truth.size());
    for (std::size_t index = 0; index < truth.size(); ++index)
    {
      CHECK(band->front().Vertices[index].Position == truth.at(index));
    }

    /**
     *  The same closed by a half circle back: seven vertices, every one the truth's to the bit.
     *  Each is made over in two doubles a number by how it was made, where the doubles that found
     *  it lay a step or two off.
     */
    const Polyline shut = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 1},
                                        {.Position = {.X = 0, .Y = 2}, .Bulge = -1},
                                        {.Position = {.X = 1, .Y = 1}, .Bulge = -quarter}},
                           .Closed   = true};
    const auto     ring =
        Buffer(shut, {.Distance = 0.5, .Sides = Side::Both, .Gaps = Gap::Chamfer}, OffsetFine);
    REQUIRE(ring.has_value());
    REQUIRE(ring->size() == 1);
    const std::array<Point, 7> round = {Point{.X = -0.5, .Y = -0.41421356237309503},
                                        Point{.X =   1.4873430246604409, .Y = 0.8055501941530117},
                                        Point{.X =  -0.5675030147051942, .Y =  2.063510509912521},
                                        Point{.X = -0.29438659539145867, .Y = 1.4041491462985238},
                                        Point{.X =   0.4930977365747153, .Y = 1.0827926457174348},
                                        Point{.X =                  0.5, .Y =                  1},
                                        Point{.X = -0.49999999999999994, .Y = 0.9999999940726246}};
    REQUIRE(ring->front().Vertices.size() == round.size());
    for (std::size_t index = 0; index < round.size(); ++index)
    {
      CHECK(ring->front().Vertices[index].Position == round.at(index));
    }
  }

  TEST_CASE("buffer: two square ends whose far sides are all but one line meet at a vertex")
  {
    /**
     *  Two quarter arcs and their ends squared off: the far side of the one square and a side of
     *  the other run within a step of a double of one line, as the arcs' bulge, tan(pi / 8) as a
     *  double holds it, makes them. They are two lines and cross at one point, a vertex of the
     *  edge; taken for one line by their rounded ends, the edge turned at the corner of a square
     *  instead, a fifth of the distance along.
     */
    const double   quarter = 0x1.a827999fcef33p-2;
    const Polyline open    = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = -quarter},
                                           {.Position = {.X = 0, .Y = 1}, .Bulge = -quarter},
                                           {.Position = {.X = 1, .Y = 1}}}};
    const auto     band =
        Buffer(open, {.Distance = 1, .Sides = Side::Both, .Gaps = Gap::Fillet, .Ends = Cap::Square},
               OffsetFine);
    REQUIRE(band.has_value());
    REQUIRE(band->size() == 1);
    const std::array<Point, 6> truth = {
        Point{.X =   -0.7071067811865476, .Y = 1.7071067811865475  },
        Point{.X =   -0.7071067811865475, .Y = -0.7071067811865476 },
        Point{.X = 9.937333212232101e-17, .Y = -std::numbers::sqrt2},
        Point{.X =    1.2071067811865475, .Y = -0.2071067811865475 },
        Point{.X =     2.414213562373095, .Y = 0.9999999999999999  },
        Point{.X = 1.7071067811865477, .Y = 1.7071067811865475}};
    REQUIRE(band->front().Vertices.size() == truth.size());
    for (std::size_t index = 0; index < truth.size(); ++index)
    {
      CHECK(Distance(band->front().Vertices[index].Position, truth.at(index)) <= 0x1p-50);
    }
  }

  TEST_CASE("buffer: a slight arc run out and back is banded once, its two runs one")
  {
    /**
     *  An arc that bows a hundred thousandth of its chord, run out and then back along itself: each
     *  run's moved arcs lie just the distance from the other run, whose band is on their near side
     *  and not past them, so they are the edge. The slight arc's center is too far to measure from,
     *  and which side of it a point lies on is read from its ends and its bulge.
     */
    const Polyline twice = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 1e-6},
                                         {.Position = {.X = 10, .Y = 0}, .Bulge = -1e-6}},
                            .Closed   = true};
    const auto     band  = Buffer(twice, {.Distance = 1, .Sides = Side::Both}, OffsetFine);
    REQUIRE(band.has_value());
    REQUIRE(band->size() == 1);
    CHECK(band->front().Vertices.size() == 4);
    CHECK(std::abs(std::abs(Area(band->front())) - (20 + std::numbers::pi)) < 1e-6);
  }

  TEST_CASE("buffer: a flat end's line that just cuts a circle leaves the edge to the circle")
  {
    /**
     *  A quarter arc out and the same arc back, the ends flat: the line across the ends leans by a
     *  step of a double, as the arcs' bulge makes it, and starts 5.5e-17 within the circle about
     *  the far vertex, which it leaves 1.2e-8 on. Up to there the circle is the edge. How far the
     *  line lies from that vertex was read in doubles, which tell 5.5e-17 from nothing no better
     *  than the line from the circle, and the edge ran along the line from its start.
     */
    const double   quarter = 0x1.a827999fcef33p-2;
    const Polyline open    = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = quarter},
                                           {.Position = {.X = 1, .Y = 1}, .Bulge = -quarter},
                                           {.Position = {.X = 0, .Y = 0}}}};
    const auto     band =
        Buffer(open, {.Distance = 1, .Sides = Side::Both, .Gaps = Gap::Fillet, .Ends = Cap::Flat},
               OffsetFine);
    REQUIRE(band.has_value());
    REQUIRE(band->size() == 1);
    REQUIRE(band->front().Vertices.size()    == 3);
    CHECK(band->front().Vertices[1].Position == Point{.X = 2, .Y = 1});
    CHECK(band->front().Vertices[2].Position
          == Point{.X = 7.0267556179791788e-17, .Y = 0.99999998814524937});
    CHECK(std::abs(band->front().Vertices[1].Bulge - 1.0000000059273753) <= 0x1p-50);
  }

  TEST_CASE(
      "offset: a point two carriers share is stepped onto both, but not where they all but touch")
  {
    /**
     *  Two lines square to each other: the point doubles found a little off is put on both.
     */
    const Belt::Carrier level   = {.WayX = {.High = 1}};
    const Belt::Carrier upright = {.X = {.High = 2}, .WayY = {.High = 1}};
    CHECK(Belt::Sharpened(level, upright, {.X = 2 + 0x1p-40, .Y = 0x1p-40})
          == Point{.X = 2, .Y = 0});
    /**
     *  A circle and a line that touches it run one way where they touch: no step leads onto both,
     *  and the point is left as found. Nor is a step taken that leaps away from the point, as the
     *  two give a little along from where they touch.
     */
    const Belt::Carrier round = {.Round = true, .Radius = {.High = 1}};
    const Belt::Carrier touch = {.Y = {.High = 1}, .WayX = {.High = 1}};
    CHECK(Belt::Sharpened(round, touch, {.X = 0, .Y = 1})    == Point{.X = 0, .Y = 1});
    CHECK(Belt::Sharpened(round, touch, {.X = 1e-3, .Y = 1}) == Point{.X = 1e-3, .Y = 1});
  }
} // namespace Geoxy
