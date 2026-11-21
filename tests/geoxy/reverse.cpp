/***************************************************************************************************
 *  @file      tests/geoxy/reverse.cpp
 *  @brief     Tests of each curve run the other way: its ends swapped to the bit, and twice the
 *             same curve.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/reverse.hpp"

#include "curves.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <array>
#include <bit>
#include <cstdint>
#include <optional>
#include <variant>

namespace Geoxy
{
  namespace
  {
    /**
     *  Whether two points are within a distance.
     *
     *  @param  lhs      One point.
     *  @param  rhs      The other.
     *  @param  distance The distance.
     *  @return          Whether they are.
     */
    auto PointsNear(Point lhs, Point rhs, double distance) -> bool
    {
      return Distance(lhs, rhs) <= distance;
    }

    /**
     *  Whether a curve reversed runs back over it: its ends swapped to the bit and each point of it
     *  at the mirrored parameter within a distance, and twice reversed it is itself.
     *
     *  @tparam ShapeT   The curve's kind.
     *  @tparam MirrorT  What gives the reversed curve's parameter.
     *  @param  shape    The curve.
     *  @param  mirror   The parameter of the reversed curve's point at a parameter of the curve.
     *  @param  distance How far apart the points within may be.
     *  @return          Whether it does.
     */
    template <typename ShapeT, typename MirrorT>
    auto RunsBack(const ShapeT &shape, const MirrorT &mirror, double distance) -> bool
    {
      const ShapeT reversed = Reversed(shape);
      const Range  domain   = Domain(shape);
      bool         near     = At(reversed, mirror(domain.End)) == At(shape, domain.End);
      for (const double along : {0.1, 0.3, 0.5, 0.7, 0.9})
      {
        const double place = domain.Start + ((domain.End - domain.Start) * along);
        near = near && PointsNear(At(reversed, mirror(place)), At(shape, place), distance);
      }
      return near && Reversed(reversed) == shape;
    }

    /**
     *  Whether two angles are one to the bit: held in one unit, as one double.
     *
     *  @param  lhs One angle.
     *  @param  rhs The other.
     *  @return     Whether they are.
     */
    auto Twin(Angle lhs, Angle rhs) -> bool
    {
      const auto held = [](Angle angle) -> std::uint64_t
      {
        return std::bit_cast<std::uint64_t>(angle.HeldInDegrees() ? angle.InDegrees()
                                                                  : angle.InRadians());
      };
      return lhs.HeldInDegrees() == rhs.HeldInDegrees() && held(lhs) == held(rhs);
    }

    /**
     *  The parameter at the other end of a unit domain.
     *
     *  @param  along The parameter.
     *  @return       1 less it.
     */
    auto Opposite(double along) -> double
    {
      return 1.0 - along;
    }
  } // namespace

  TEST_CASE("reverse: a line, an xline and the Bézier curves, twice the same to the bit")
  {
    const Line line = {.Start = {.X = 5e6, .Y = 1}, .End = {.X = 5e6 + 3, .Y = -2}};
    CHECK(Reversed(line) == Line{.Start = line.End, .End = line.Start});
    CHECK(RunsBack(line, Opposite, 1e-9));
    const Xline xline = {.Through = {.X = 1, .Y = 2}, .Way = {.X = 3, .Y = -4}};
    CHECK(Reversed(xline)           == Xline{.Through = xline.Through, .Way = {.X = -3, .Y = 4}});
    CHECK(At(Reversed(xline), -2)   == At(xline, 2));
    CHECK(Reversed(Reversed(xline)) == xline);
    const Bezier::Quadratic quadratic = {
        .Start = {.X = 0, .Y = 0}, .Control = {.X = 1, .Y = 3}, .End = {.X = 4, .Y = 0.5}};
    CHECK(RunsBack(quadratic, Opposite, 1e-15));
    CHECK(Reversed(quadratic).Control == quadratic.Control);
    CHECK(RunsBack(Tests::Wave, Opposite, 1e-14));
    CHECK(Reversed(Tests::Wave).Out == Tests::Wave.In);
  }

  TEST_CASE("reverse: an arc starts at its end angle, twice the same where the angles add exactly")
  {
    const Arc arc      = {.Center = {.X = 5e6, .Y = -3e6},
                          .Radius = 2.5,
                          .Start  = Degrees(30.25),
                          .Sweep  = Degrees(-200.5)};
    const Arc reversed = Reversed(arc);
    CHECK(reversed.Start == Degrees(-170.25));
    CHECK(reversed.Sweep == Degrees(200.5));
    CHECK(RunsBack(arc, Opposite, 1e-9));

    const Arc radians = {
        .Center = {.X = 1, .Y = 2}, .Radius = -0.75, .Start = Radians(0.3), .Sweep = Radians(0.7)};
    const Arc back = Reversed(radians);
    CHECK(At(back, 0) == At(radians, 1));
    CHECK(PointsNear(At(back, 1), At(radians, 0), 1e-15));
    CHECK(PointsNear(At(Reversed(back), 0.25), At(radians, 0.25), 1e-15));

    const Ellipse ellipse = {.Map   = Rotation(Degrees(20), {.X = 4, .Y = 1}) * Scaling(3, 1),
                             .Start = Degrees(-45),
                             .Sweep = Degrees(270)};
    CHECK(RunsBack(ellipse, Opposite, 1e-14));
    CHECK(Reversed(ellipse).Map == ellipse.Map);
  }

  TEST_CASE("reverse: an arc starts at its start and sweep added, to the bit, as a region's whole "
            "circle is turned about")
  {
    const std::array arcs = {
        Arc{.Center = {.X = 5e6, .Y = -3e6},
            .Radius = 2.5,
            .Start  = Degrees(30.25),
            .Sweep  = Degrees(-200.5)},
        Arc{.Center = {.X = 1, .Y = 2},
            .Radius = 0.75,
            .Start  = Radians(0.3),
            .Sweep  = Radians(0.7)},
        Arc{.Center = {}, .Radius = 1, .Start = Degrees(10.1), .Sweep = Radians(-2.9)},
        Arc{.Center = {}, .Radius = 3, .Start = Radians(1e-300), .Sweep = Degrees(359.123456789)},
        Arc{.Center = {}, .Radius = 1e-9, .Start = Degrees(-359.999999), .Sweep = Degrees(360)},
    };
    for (const Arc &arc : arcs)
    {
      const Arc reversed = Reversed(arc);
      CHECK(Twin(reversed.Start, arc.Start + arc.Sweep));
      CHECK(Twin(reversed.Sweep, -arc.Sweep));
    }
  }

  TEST_CASE("reverse: a polyline, open or closed, its bulges negated and moved")
  {
    const Polyline open     = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0.5},
                                            {.Position = {.X = 4, .Y = 0}, .Bulge = 0},
                                            {.Position = {.X = 4, .Y = 3}, .Bulge = -1},
                                            {.Position = {.X = 0, .Y = 3}, .Bulge = 0.25}},
                               .Closed   = false};
    const Polyline reversed = Reversed(open);
    REQUIRE(reversed.Vertices.size() == 4);
    CHECK(reversed.Vertices[0]       == Vertex{.Position = {.X = 0, .Y = 3}, .Bulge = 1});
    CHECK(reversed.Vertices[1]       == Vertex{.Position = {.X = 4, .Y = 3}, .Bulge = 0});
    CHECK(reversed.Vertices[2]       == Vertex{.Position = {.X = 4, .Y = 0}, .Bulge = -0.5});
    CHECK(reversed.Vertices[3]       == Vertex{.Position = {.X = 0, .Y = 0}, .Bulge = -0.25});
    CHECK(RunsBack(open, [](double along) -> double { return 3.0 - along; }, 1e-14));

    Polyline closed = open;
    closed.Closed   = true;
    CHECK(Reversed(closed).Vertices[0] == Vertex{.Position = {.X = 0, .Y = 0}, .Bulge = -0.25});
    CHECK(Reversed(closed).Vertices[1] == Vertex{.Position = {.X = 0, .Y = 3}, .Bulge = 1});
    CHECK(RunsBack(closed, [](double along) -> double { return 4.0 - along; }, 1e-14));

    CHECK(Reversed(Polyline{}) == Polyline{});
    const Polyline point = {.Vertices = {{.Position = {.X = 1, .Y = 1}, .Bulge = 2}},
                            .Closed   = false};
    CHECK(Reversed(Reversed(point)) == point);
  }

  TEST_CASE("reverse: a spline, its knots negated, twice the same to the bit")
  {
    const auto negated = [](double along) -> double
    {
      return -along;
    };
    CHECK(RunsBack(Tests::Clamped(), negated, 1e-14));
    CHECK(RunsBack(Tests::Quarter(), negated, 1e-15));
    const Spline uneven = Tests::Uneven();
    CHECK(RunsBack(uneven, negated, 1e-14));
    CHECK(Domain(Reversed(uneven))         == Range{.Start = -2, .End = 0});
    CHECK(Reversed(uneven).Weights.front() == doctest::Approx(1));
    CHECK(Reversed(uneven).Weights[2]      == doctest::Approx(3));
  }

  TEST_CASE("reverse: any curve as its own kind, but a ray")
  {
    const std::array<Curve, 3> curves = {
        Curve{Line{.Start = {.X = 1, .Y = 0}, .End = {.X = 2, .Y = 0}}},
        Curve{Circle({.X = 1, .Y = 1}, 2)},
        Curve{Tests::Clamped()},
    };
    for (const Curve &curve : curves)
    {
      const std::optional<Curve> reversed = Reversed(curve);
      REQUIRE(reversed);
      CHECK(reversed->index()   == curve.index());
      CHECK(Reversed(*reversed) == curve);
    }
    CHECK_FALSE(Reversed(Curve{Ray{.Start = {}, .Way = {.X = 1, .Y = 0}}}));
  }
} // namespace Geoxy
