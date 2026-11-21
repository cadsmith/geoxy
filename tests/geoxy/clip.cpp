/***************************************************************************************************
 *  @file      tests/geoxy/clip.cpp
 *  @brief     Tests of curves cut to a box, against dense samples of their points.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/clip.hpp"

#include "curves.hpp"
#include "stiff.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/box.hpp"
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

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  The box from 0,0 to 2,1.
     */
    constexpr Box Wide = {.Min = {.X = 0, .Y = 0}, .Max = {.X = 2, .Y = 1}};

    /**
     *  Whether a parameter lies in one of a run of stretches.
     *
     *  @param  ranges The stretches.
     *  @param  along  The parameter.
     *  @return        Whether it does.
     */
    auto Within(const std::vector<Range> &ranges, double along) -> bool
    {
      return std::ranges::any_of(ranges, [along](const Range &range) -> bool
                                 { return range.Start <= along && along <= range.End; });
    }

    /**
     *  Whether a curve's stretches in a box are right: in order and apart, every dense sample
     *  clearly within the box among them, and every one clearly outside it not.
     *
     *  @param  curve The curve.
     *  @param  box   The box.
     *  @return       Whether they are.
     */
    auto Agrees(const Curve &curve, const Box &box) -> bool
    {
      const std::vector<Range> ranges  = Clip(curve, box);
      const Range              domain  = Domain(curve);
      const Box                inner   = Padded(box, -1e-9);
      const Box                outer   = Padded(box, 1e-9);
      constexpr int            Samples = 20'000;
      for (std::size_t index = 0; index < ranges.size(); ++index)
      {
        const bool ordered = ranges[index].Start < ranges[index].End
                          && (index == 0 || ranges[index - 1].End < ranges[index].Start);
        if (!ordered)
        {
          return false;
        }
      }
      for (int index = 0; index <= Samples; ++index)
      {
        const double along = domain.Start + ((domain.End - domain.Start) * index / Samples);
        const Point  point = At(curve, along);
        if ((Contains(inner, point) && !Within(ranges, along))
            || (!Contains(outer, point) && Within(ranges, along)))
        {
          return false;
        }
      }
      return true;
    }

    /**
     *  Whether two stretches are equal to within a few bits.
     *
     *  @param  range    The stretch found.
     *  @param  expected The stretch expected.
     *  @return          Whether they are.
     */
    auto Near(const Range &range, const Range &expected) -> bool
    {
      return std::abs(range.Start - expected.Start) < 1e-12
          && std::abs(range.End - expected.End)     < 1e-12;
    }
  } // namespace

  TEST_CASE("clip: a line, a ray and an xline by their slabs")
  {
    const Line across = {.Start = {.X = -1, .Y = 0.5}, .End = {.X = 3, .Y = 0.5}};
    CHECK(Clip(across, Wide) == std::vector<Range>{{.Start = 0.25, .End = 0.75}});

    const Line inside = {.Start = {.X = 0.5, .Y = 0.5}, .End = {.X = 1, .Y = 1}};
    CHECK(Clip(inside, Wide) == std::vector<Range>{{.Start = 0, .End = 1}});

    const Line edge = {.Start = {.X = -1, .Y = 1}, .End = {.X = 3, .Y = 1}};
    CHECK(Clip(edge, Wide) == std::vector<Range>{{.Start = 0.25, .End = 0.75}});

    const Line above = {.Start = {.X = -1, .Y = 2}, .End = {.X = 3, .Y = 2}};
    CHECK(Clip(above, Wide).empty());
    CHECK(Clip(Line{.Start = {.X = -1, .Y = -1}, .End = {.X = 3, .Y = -1}}, Wide).empty());

    const Line point = {.Start = {.X = 1, .Y = 1}, .End = {.X = 1, .Y = 1}};
    CHECK(Clip(point, Wide) == std::vector<Range>{{.Start = 0, .End = 1}});

    const Line corner = {.Start = {.X = 2, .Y = 1}, .End = {.X = 3, .Y = 2}};
    CHECK(Clip(corner, Wide).empty());
    CHECK(Clip(inside, Box{}).empty());

    const Ray ray = {.Start = {.X = 1, .Y = 0.5}, .Way = {.X = 2, .Y = 0}};
    CHECK(Clip(ray, Wide) == std::vector<Range>{{.Start = 0, .End = 0.5}});

    const Ray away = {.Start = {.X = 2, .Y = 1}, .Way = {.X = 1, .Y = 1}};
    CHECK(Clip(away, Wide).empty());

    const Xline xline = {.Through = {.X = 1, .Y = 0.5}, .Way = {.X = 1, .Y = 1}};
    CHECK(Clip(xline, Wide) == std::vector<Range>{{.Start = -0.5, .End = 0.5}});

    const double infinity = std::numeric_limits<double>::infinity();
    const Box    endless  = {.Min = {.X = -infinity, .Y = 0}, .Max = {.X = infinity, .Y = 1}};
    const Xline  level    = {.Through = {.X = 1, .Y = 0.5}, .Way = {.X = 1, .Y = 0}};
    CHECK(Clip(level, endless) == std::vector<Range>{{.Start = -infinity, .End = infinity}});
  }

  TEST_CASE("clip: an arc where it crosses the box's edges, every turn of it")
  {
    const Box                right  = {.Min = {.X = 0, .Y = -2}, .Max = {.X = 2, .Y = 2}};
    const std::vector<Range> halves = Clip(Circle({}, 1), right);
    REQUIRE(halves.size() == 2);
    CHECK(Near(halves[0], {.Start = 0, .End = 0.25}));
    CHECK(Near(halves[1], {.Start = 0.75, .End = 1}));

    CHECK(Clip(Circle({}, 1), Padded(right, 5)) == std::vector<Range>{{.Start = 0, .End = 1}});
    CHECK(Clip(Circle({.X = 9, .Y = 9}, 1), right).empty());
    CHECK(Clip(Circle({}, 1), Box{}).empty());

    const Arc point = {.Center = {}, .Radius = 1, .Start = Degrees(10), .Sweep = {}};
    CHECK(Clip(point, right) == std::vector<Range>{{.Start = 0, .End = 1}});

    const Box below              = {.Min = {.X = -2, .Y = -2}, .Max = {.X = 2, .Y = -0.5}};
    const Arc under              = {.Center = {}, .Radius = 1, .Start = {}, .Sweep = Degrees(-180)};
    const std::vector<Range> low = Clip(under, below);
    REQUIRE(low.size() == 1);
    CHECK(Near(low[0], {.Start = 1.0 / 6.0, .End = 5.0 / 6.0}));

    const Arc thrice = {.Center = {}, .Radius = 1, .Start = {}, .Sweep = Turn * 3.0};
    CHECK(Clip(thrice, right).size() == 4);
    CHECK(Agrees(thrice, right));
    CHECK(Agrees(Arc{.Center = {.X = 1, .Y = 0.5},
                     .Radius = 0.8,
                     .Start  = Degrees(30),
                     .Sweep  = Degrees(-300)},
                 Wide));

    const Arc endless             = {.Center = {}, .Radius = 1, .Start = {}, .Sweep = Turn * 1e5};
    const std::vector<Range> many = Clip(endless, right);
    CHECK(many.size() > 65'536);
    CHECK(std::ranges::is_sorted(many, {}, &Range::Start));
  }

  TEST_CASE("clip: an ellipse, turned, folded and flattened")
  {
    const Ellipse turned = {
        .Map   = Translation({.X = 1, .Y = 0.5}) * Rotation(Degrees(30)) * Scaling(1.5, 0.4),
        .Start = Degrees(20),
        .Sweep = Degrees(300)};
    CHECK(Agrees(turned, Wide));
    CHECK(!Clip(turned, Wide).empty());

    const Ellipse folded = {.Map = Scaling(-1.5, 0.8), .Start = {}, .Sweep = Turn};
    CHECK(Agrees(folded, Wide));

    const Ellipse flat = {.Map = Scaling(3, 0), .Start = {}, .Sweep = Turn};
    CHECK(Agrees(flat, Wide));
  }

  TEST_CASE("clip: a polyline, its stretches running on across vertices")
  {
    const Polyline square = {.Vertices = {{.Position = {.X = 0.5, .Y = 0.2}, .Bulge = 0},
                                          {.Position = {.X = 1.5, .Y = 0.2}, .Bulge = 0},
                                          {.Position = {.X = 1.5, .Y = 0.8}, .Bulge = 0},
                                          {.Position = {.X = 0.5, .Y = 0.8}, .Bulge = 0}},
                             .Closed   = true};
    CHECK(Clip(square, Wide) == std::vector<Range>{{.Start = 0, .End = 4}});

    const Polyline bulged = {.Vertices = {{.Position = {.X = -0.5, .Y = 0.2}, .Bulge = 0.8},
                                          {.Position = {.X = 1.5, .Y = 0.2}, .Bulge = -1.5},
                                          {.Position = {.X = 2.5, .Y = 0.9}, .Bulge = 0},
                                          {.Position = {.X = 0.5, .Y = 1.5}, .Bulge = 0}},
                             .Closed   = true};
    CHECK(Agrees(bulged, Wide));
    CHECK(Clip(bulged, Wide).size() > 1);

    const Polyline lone = {.Vertices = {{.Position = {.X = 1, .Y = 1}, .Bulge = 0}},
                           .Closed   = false};
    CHECK(Clip(lone, Wide) == std::vector<Range>{{.Start = 0, .End = 0}});
    CHECK(Clip(lone, Padded(Wide, -0.5)).empty());
    CHECK(Clip(Polyline{}, Wide).empty());
  }

  TEST_CASE("clip: splines and Bézier curves by the roots of their pieces")
  {
    const Box middle = {.Min = {.X = 0.5, .Y = 0.2}, .Max = {.X = 6.5, .Y = 1.5}};
    CHECK(Agrees(Tests::Uneven(), middle));
    CHECK(Agrees(Tests::Quarter(), {.Min = {.X = 0, .Y = 0}, .Max = {.X = 0.6, .Y = 2}}));
    CHECK(Agrees(Tests::Wave, {.Min = {.X = 0, .Y = -0.2}, .Max = {.X = 5, .Y = 0.5}}));

    const Bezier::Quadratic arch = {
        .Start = {.X = 0, .Y = 0}, .Control = {.X = 1, .Y = 2}, .End = {.X = 2, .Y = 0}};
    const std::vector<Range> feet =
        Clip(arch, {.Min = {.X = -1, .Y = -1}, .Max = {.X = 3, .Y = 0.5}});
    const double foot = (1 - std::sqrt(0.5)) / 2;
    REQUIRE(feet.size() == 2);
    CHECK(Near(feet[0], {.Start = 0, .End = foot}));
    CHECK(Near(feet[1], {.Start = 1 - foot, .End = 1}));

    const Spline bent = {.Degree  = 1,
                         .Knots   = {0, 0, 1, 2, 2},
                         .Poles   = {{.X = -1, .Y = -1}, {.X = 0, .Y = 0}, {.X = 1, .Y = 1}},
                         .Weights = {}};
    const Box    half = {.Min = {.X = 0, .Y = -5}, .Max = {.X = 5, .Y = 5}};
    CHECK(Clip(bent, half) == std::vector<Range>{{.Start = 1, .End = 2}});
  }

  /**
   *  With control points near the largest double a piece's slope overflows, and Newton's step from
   *  it is nothing: the root finder then went on a double at a time, for as long as anyone waited.
   */
  TEST_CASE("clip: a spline of any weights, its stretches as exact arithmetic gives them")
  {
    /**
     *  Before, stretches beside a weight of 2^700 or more were missed.
     */
    for (const Tests::Boxed &boxed : Tests::Boxes(Tests::Clips))
    {
      CAPTURE(boxed.Curve.Weights[1]);
      CAPTURE(boxed.Within.Min);
      CAPTURE(boxed.Within.Max);
      const std::vector<Range> ranges = Clip(boxed.Curve, boxed.Within);
      REQUIRE(2 * ranges.size() == boxed.Ends.size());
      for (std::size_t index = 0; index < ranges.size(); ++index)
      {
        CHECK(Tests::Tells(ranges[index].Start, boxed.Ends[2 * index]));
        CHECK(Tests::Tells(ranges[index].End, boxed.Ends[(2 * index) + 1]));
      }
    }

    /**
     *  Beside a weight of 2^100 the legs of the curve lie within a sliver of its parameter at each
     *  end. A box about the far leg holds a stretch no double can tell from the end: none is given.
     *  One about the near leg holds a stretch from the start, told to its last digits.
     */
    const Spline stiff = Tests::Stiff("conic", 100);
    CHECK(Clip(stiff, {.Min = {.X = 1.25, .Y = -1}, .Max = {.X = 3, .Y = 0.75}}).empty());
    const std::vector<Range> near =
        Clip(stiff, {.Min = {.X = -1, .Y = -1}, .Max = {.X = 0.5, .Y = 2}});
    REQUIRE(near.size() == 1);
    CHECK(near[0].Start == 0);
    CHECK(near[0].End   == doctest::Approx(0x1p-101).epsilon(1e-14));
    CHECK(Clip(stiff, Box{}).empty());

    /**
     *  A box's edge through the very place two pieces of the cut spline meet: the crossing, found
     *  at the end of one and the start of the next, is one crossing.
     */
    const std::vector<Range> half =
        Clip(Tests::Stiff("conic", 3), {.Min = {.X = -1, .Y = -1}, .Max = {.X = 1, .Y = 2}});
    REQUIRE(half.size() == 1);
    CHECK(half[0].Start == 0);
    CHECK(half[0].End   == 0.5);

    /**
     *  And one through the spline's very end: the stretch runs to the end, the end on the edge.
     */
    const std::vector<Range> last =
        Clip(Tests::Stiff("conic", 3), {.Min = {.X = 1.5, .Y = -1}, .Max = {.X = 2, .Y = 1}});
    REQUIRE(last.size() == 1);
    CHECK(last[0].Start > 0.5);
    CHECK(last[0].End == 1);
  }

  TEST_CASE("clip: a curve with control points near the largest double ends")
  {
    const Box near = {.Min = {.X = -1, .Y = -1}, .Max = {.X = 0, .Y = 1}};
    for (const double far : {1e308, -1e308})
    {
      CAPTURE(far);
      const Bezier::Cubic     cubic = {.Start = {.X = 0, .Y = 0},
                                       .Out   = {.X = far, .Y = 1},
                                       .In    = {.X = 1, .Y = far},
                                       .End   = {.X = 2, .Y = 0}};
      const Bezier::Quadratic arch  = {
          .Start = {.X = 0, .Y = 0}, .Control = {.X = far, .Y = 1}, .End = {.X = 2, .Y = 0}};
      for (const Range &kept : Clip(cubic, near))
      {
        CHECK(kept.Start <= kept.End);
      }
      for (const Range &kept : Clip(arch, near))
      {
        CHECK(kept.Start <= kept.End);
      }
    }
  }

  TEST_CASE("clip: any curve as its own kind")
  {
    const Curve line = Line{.Start = {.X = -1, .Y = 0.5}, .End = {.X = 3, .Y = 0.5}};
    CHECK(Clip(line, Wide) == std::vector<Range>{{.Start = 0.25, .End = 0.75}});
    CHECK(Clip(Curve{Tests::Clamped()}, Wide) == Clip(Tests::Clamped(), Wide));
  }
} // namespace Geoxy
