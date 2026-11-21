/***************************************************************************************************
 *  @file      tests/geoxy/area.cpp
 *  @brief     Tests of the area a curve encloses and its centroid, against closed forms and an
 *             independent quadrature of Green's theorem.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/area.hpp"

#include "curves.hpp"
#include "real.hpp"
#include "splines.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  What a curve encloses with its chord, by Simpson's rule on Green's theorem over its
     *  stretches: an evaluator that shares nothing with the library's.
     */
    struct GreenOracle
    {
      double Area     = 0.0; /**< The signed area. */
      Point  Centroid = {};  /**< The centroid.    */
    };

    /**
     *  What a curve encloses with its chord, by Simpson's rule.
     *
     *  @param  curve  The curve.
     *  @param  ranges The stretches it is smooth over, from its start to its end.
     *  @return        The area and the centroid.
     */
    auto GreenSimpson(const Curve &curve, const std::vector<Range> &ranges) -> GreenOracle
    {
      constexpr int Panels = 20'000;
      double        area   = 0.0;
      Vector        moment = {};
      for (const Range range : ranges)
      {
        const double step = (range.End - range.Start) / Panels;
        for (int index = 0; index <= Panels; ++index)
        {
          double weight = index % 2 == 1 ? 4.0 : 2.0;
          if (index == 0 || index == Panels)
          {
            weight = 1.0;
          }
          const double along  = range.Start + (step * index);
          const Vector place  = At(curve, along) - Point{};
          const double cross  = Cross(place, Derivative(curve, along)) * weight * step / 3.0;
          area               += cross / 2.0;
          moment              = moment + (place * (cross / 3.0));
        }
      }
      const Range  domain  = Domain(curve);
      const Vector end     = At(curve, domain.End) - Point{};
      const Vector start   = At(curve, domain.Start) - Point{};
      const double chord   = Cross(end, start);
      area                += chord / 2.0;
      moment               = moment + ((end + start) * (chord / 6.0));
      return {.Area = area, .Centroid = Point{} + (moment / area)};
    }

    /**
     *  Whether a centroid is a point to within a tolerance.
     *
     *  @param  centroid  The centroid, if any.
     *  @param  point     The point expected.
     *  @param  tolerance How far it may be.
     *  @return           Whether it is.
     */
    auto Near(const std::optional<Point> &centroid, Point point, double tolerance) -> bool
    {
      return centroid.has_value() && Distance(*centroid, point) <= tolerance;
    }

    /**
     *  A spline's reverse: the same curve run the other way, its poles and weights turned about and
     *  its knots turned about and negated.
     *
     *  @param  spline The spline.
     *  @return        Its reverse.
     */
    auto AreaReversed(Spline spline) -> Spline
    {
      std::ranges::reverse(spline.Poles);
      std::ranges::reverse(spline.Weights);
      std::ranges::reverse(spline.Knots);
      for (double &knot : spline.Knots)
      {
        knot = -knot;
      }
      return spline;
    }

    /**
     *  How many steps of a double a value is from a reference.
     *
     *  @param  value     The value.
     *  @param  reference The reference.
     *  @return           The distance in steps of the reference.
     */
    auto Steps(double value, double reference) -> double
    {
      const double magnitude = std::abs(reference);
      const double step =
          std::nextafter(magnitude, std::numeric_limits<double>::infinity()) - magnitude;
      return std::abs(value - reference) / step;
    }

    /**
     *  An arc's sweep, where its center sits, and what it encloses with its chord to 30 digits.
     */
    struct ArcCase
    {
      double Sweep  = 0.0; /**< The sweep, in radians. */
      double Offset = 0.0; /**< The center's x and y.  */
      double Area   = 0.0; /**< The area, rounded.     */
    };
  } // namespace

  TEST_CASE("area: what is straight or endless encloses nothing")
  {
    const Line line = {.Start = {}, .End = {.X = 1, .Y = 1}};
    CHECK(Area(line) == 0);
    CHECK(!Centroid(line).has_value());
    CHECK(Area(Ray{.Start = {}, .Way = {.X = 1, .Y = 0}}) == 0);
    CHECK(!Centroid(Ray{.Start = {}, .Way = {.X = 1, .Y = 0}}).has_value());
    CHECK(Area(Xline{.Through = {}, .Way = {.X = 1, .Y = 0}}) == 0);
    CHECK(!Centroid(Xline{.Through = {}, .Way = {.X = 1, .Y = 0}}).has_value());
  }

  TEST_CASE("area: a polyline by the shoelace, open ones closed, its bulges as arcs")
  {
    const Polyline square = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0},
                                          {.Position = {.X = 2, .Y = 0}, .Bulge = 0},
                                          {.Position = {.X = 2, .Y = 2}, .Bulge = 0},
                                          {.Position = {.X = 0, .Y = 2}, .Bulge = 0}},
                             .Closed   = true};
    CHECK(Area(square)     == 4);
    CHECK(Centroid(square) == Point{.X = 1, .Y = 1});

    Polyline far = square;
    for (Vertex &vertex : far.Vertices)
    {
      vertex.Position = vertex.Position + Vector{.X = 5e6, .Y = -3e6};
    }
    CHECK(Area(far)     == 4);
    CHECK(Centroid(far) == Point{.X = 5e6 + 1, .Y = -3e6 + 1});

    Polyline backward = square;
    backward.Vertices = {square.Vertices[3], square.Vertices[2], square.Vertices[1],
                         square.Vertices[0]};
    CHECK(Area(backward)     == -4);
    CHECK(Centroid(backward) == Point{.X = 1, .Y = 1});

    Polyline open = square;
    open.Closed   = false;
    CHECK(Area(open) == 4);

    Polyline bulged          = square;
    bulged.Vertices[1].Bulge = 1;
    CHECK(Area(bulged) == doctest::Approx(4 + (std::numbers::pi / 2)).epsilon(1e-15));
    const double lobe = (4 + (std::numbers::pi / 2));
    CHECK(
        Near(Centroid(bulged),
             {.X = ((4 * 1) + ((std::numbers::pi / 2) * (2 + (4 / (3 * std::numbers::pi))))) / lobe,
              .Y = 1},
             1e-14));

    CHECK(Area(Polyline{}) == 0);
    CHECK(!Centroid(Polyline{}).has_value());
    CHECK(!Centroid(
               Polyline{.Vertices = {{.Position = {.X = 3, .Y = 4}, .Bulge = 0}}, .Closed = false})
               .has_value());
  }

  /**
   *  A polyline's area was gathered about its first vertex from rounded differences, its bulges'
   *  segments of their circles in doubles: four steps of a double off the exact area of this one.
   *  The differences are now exact, the cross products and the bulges' segments in double-double,
   *  and the sum rounded once.
   */
  TEST_CASE("area: a polyline the double nearest its exact area")
  {
    const Polyline polyline = {
        .Vertices = {{.Position = {.X = 2.419068749223258, .Y = 0.5313921566584238},
                      .Bulge    = -0.37345047957760136},
                     {.Position = {.X = -1.3001715751056093, .Y = 2.226754315116843},
                      .Bulge    = -0.7909890461865743},
                     {.Position = {.X = -2.304949548347449, .Y = -3.1886925450701806},
                      .Bulge    = 0.2499666638798037}},
        .Closed   = true};
    CHECK(Area(polyline) == 0x1.7497948930d4ap+1);
    const Polyline shallow = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0.01},
                                           {.Position = {.X = 1, .Y = 0}, .Bulge = 0}},
                              .Closed   = true};
    CHECK(Area(shallow) == doctest::Approx(0.01 * ((1.0 / 3.0) + (1e-4 / 15.0))).epsilon(1e-12));
    const Polyline huge = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0},
                                        {.Position = {.X = 0x1p460, .Y = 0}, .Bulge = 0},
                                        {.Position = {.X = 0, .Y = 0x1p460}, .Bulge = 0}},
                           .Closed   = true};
    CHECK(Area(huge) == 0x1p919);
    const Polyline bowed = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 1},
                                         {.Position = {.X = 0x1p460, .Y = 0}, .Bulge = 0},
                                         {.Position = {.X = 0, .Y = 0x1p460}, .Bulge = 0}},
                            .Closed   = true};
    CHECK(Area(bowed) == doctest::Approx(0x1p919 * (1.0 + (std::numbers::pi / 4.0))));
    const Polyline tiny = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0},
                                        {.Position = {.X = 0x1p-460, .Y = 0}, .Bulge = 0},
                                        {.Position = {.X = 0, .Y = 0x1p-460}, .Bulge = 0}},
                           .Closed   = true};
    CHECK(Area(tiny) == 0x1p-921);
    const Polyline broken = {
        .Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0},
                     {.Position = {.X = std::numeric_limits<double>::quiet_NaN(), .Y = 0},
                      .Bulge    = 0},
                     {.Position = {.X = 0, .Y = 1}, .Bulge = 0}},
        .Closed   = true};
    CHECK(std::isnan(Area(broken)));
  }

  TEST_CASE("area: arcs and ellipses in closed form, with their chords")
  {
    CHECK(Area(Circle({.X = 3, .Y = 4}, 2))
          == doctest::Approx(4 * std::numbers::pi).epsilon(1e-15));
    CHECK(Near(Centroid(Circle({.X = 3, .Y = 4}, 2)), {.X = 3, .Y = 4}, 1e-14));

    const Arc half = {.Center = {}, .Radius = 1, .Start = {}, .Sweep = Degrees(180)};
    CHECK(Area(half) == doctest::Approx(std::numbers::pi / 2).epsilon(1e-15));
    CHECK(Near(Centroid(half), {.X = 0, .Y = 4 / (3 * std::numbers::pi)}, 1e-15));

    const Arc back = {.Center = {}, .Radius = 1, .Start = {}, .Sweep = Turn * -2};
    CHECK(Area(back) == doctest::Approx(-2 * std::numbers::pi).epsilon(1e-15));

    const Matrix  map   = Translation({.X = 1e6, .Y = 2}) * Rotation(Degrees(30)) * Scaling(3, 1);
    const Ellipse whole = {.Map = map, .Start = Degrees(40), .Sweep = Turn};
    CHECK(Area(whole) == doctest::Approx(3 * std::numbers::pi).epsilon(1e-12));
    CHECK(Near(Centroid(whole), {.X = 1e6, .Y = 2}, 1e-8));

    const Ellipse folded = {.Map = Scaling(-3, 1), .Start = {}, .Sweep = Turn};
    CHECK(Area(folded) == doctest::Approx(-3 * std::numbers::pi).epsilon(1e-15));

    for (const Ellipse &part :
         {Ellipse{.Map = Scaling(3, 1), .Start = Degrees(20), .Sweep = Degrees(130)},
          Ellipse{.Map = Translation({.X = 1, .Y = 2}) * Rotation(Degrees(30)) * Scaling(1.5, 0.4),
                  .Start = Degrees(-50),
                  .Sweep = Degrees(-250)}})
    {
      const GreenOracle oracle = GreenSimpson(part, {Domain(part)});
      CHECK(Area(part) == doctest::Approx(oracle.Area).epsilon(1e-12));
      CHECK(Near(Centroid(part), oracle.Centroid, 1e-11));
    }
    CHECK(!Centroid(Ellipse{.Map = Scaling(2, 0), .Start = {}, .Sweep = Turn}).has_value());
  }

  TEST_CASE("area: Bézier curves and splines, exact or to an accuracy")
  {
    const Bezier::Quadratic arch = {
        .Start = {.X = 0, .Y = 0}, .Control = {.X = 1, .Y = 2}, .End = {.X = 2, .Y = 0}};
    CHECK(Area(arch) == doctest::Approx(-4.0 / 3.0).epsilon(1e-15));
    CHECK(Near(Centroid(arch), {.X = 1, .Y = 0.4}, 1e-15));

    const Bezier::Cubic hump   = {.Start = {.X = 0, .Y = 0},
                                  .Out   = {.X = 1, .Y = 3},
                                  .In    = {.X = 3, .Y = 2},
                                  .End   = {.X = 4, .Y = 0}};
    const Spline        humped = {.Degree  = 3,
                                  .Knots   = {0, 0, 0, 0, 1, 1, 1, 1},
                                  .Poles   = {hump.Start, hump.Out, hump.In, hump.End},
                                  .Weights = {}};
    const GreenOracle   oracle = GreenSimpson(hump, {{.Start = 0, .End = 1}});
    CHECK(Area(hump) == doctest::Approx(oracle.Area).epsilon(1e-12));
    CHECK(Near(Centroid(hump), oracle.Centroid, 1e-12));
    CHECK(Area(humped, {}) == doctest::Approx(Area(hump)).epsilon(1e-15));
    CHECK(Near(Centroid(humped, {}), oracle.Centroid, 1e-12));
    CHECK(Area(Tests::Wave) == doctest::Approx(0).epsilon(1e-15));
    CHECK(!Centroid(Tests::Wave).has_value());

    Spline even  = humped;
    even.Weights = {2, 2, 2, 2};
    CHECK(Area(even, {}) == doctest::Approx(Area(hump)).epsilon(1e-15));

    const double segment = (std::numbers::pi / 4) - 0.5;
    CHECK(Area(Tests::Quarter(), {.Error = 1e-13}) == doctest::Approx(segment).epsilon(1e-12));
    const double reach = 4 * std::pow(std::sqrt(0.5), 3) / (3 * ((std::numbers::pi / 2) - 1));
    CHECK(Near(Centroid(Tests::Quarter(), {.Error = 1e-13}),
               {.X = reach * std::sqrt(0.5), .Y = reach * std::sqrt(0.5)}, 1e-11));

    const Spline      uneven  = Tests::Uneven();
    const GreenOracle rounded = GreenSimpson(
        uneven, {{.Start = 0, .End = 0.3}, {.Start = 0.3, .End = 1.1}, {.Start = 1.1, .End = 2}});
    CHECK(Area(uneven, {.Error = 1e-12}) == doctest::Approx(rounded.Area).epsilon(1e-11));
    CHECK(Near(Centroid(uneven, {.Error = 1e-12}), rounded.Centroid, 1e-10));
  }

  /**
   *  An arc's area was its sector about its start plus its chord, which cancel for a small sweep or
   *  a far center: sweep 1e-6 gave the wrong sign and 78 times the area. The measures audit's rows,
   *  of radius 3, now come within a step of a double of references of 30 digits.
   */
  /**
   *  A rational spline's area was an adaptive quadrature over its own parameter, wrong where its
   *  weights lie far apart: 0.77 for 100.64 with weights 2^100 apart. Every one of the corpus now
   *  comes within a few steps of a double of its area, or of the scale its rate's length gives it
   *  where its parts cancel, 110 of the 117 within 16: an unclamped span far smaller than its
   *  poles' spread, its Bézier points sums of poles far larger than it, and a curve that runs back
   *  along a line, enclosing nothing, come within a step of a double of the poles' extent squared.
   *  A spline's reverse encloses the area negated exactly, and the centroid is the same to the bit.
   */
  TEST_CASE("area: rational splines and those of any degree, a reverse negated exactly")
  {
    for (const Tests::Rated &rated : Tests::Rateds())
    {
      const Spline spline = Tests::Splined(rated);
      CAPTURE(rated.Area);
      const double area   = Area(spline, {.Error = 0});
      double       extent = 0.0;
      for (const double coordinate : rated.Points)
      {
        extent = std::max(extent, std::abs(coordinate));
      }
      CHECK(std::abs(area - rated.Area)
            <= std::max(16 * Spacing(std::max(std::abs(rated.Area), rated.Scale)),
                        Spacing(extent * extent)));
      const Spline back = AreaReversed(spline);
      CHECK(Area(back, {})     == -area);
      CHECK(Centroid(back, {}) == Centroid(spline, {}));
    }

    /**
     *  The circle of radius 15 about the origin as four arcs of integer weights, and one of them
     *  with its weights a power of two times further apart: closed forms.
     */
    const Tests::Rated &circle   = Tests::Rateds().front();
    const auto          centered = Centroid(Tests::Splined(circle), {});
    REQUIRE(centered.has_value());
    CHECK(std::abs(centered->X) <= 1e-14);
    CHECK(std::abs(centered->Y) <= 1e-14);
    const double angle   = std::acos(0.6);
    const double segment = 225 * ((2 * angle) - std::sin(2 * angle)) / 2;
    const double reach   = 4 * 15 * 0.512 / (3 * ((2 * angle) - std::sin(2 * angle)));
    for (int power = 0; power <= 500; power += 25)
    {
      CAPTURE(power);
      const double heavy  = std::ldexp(1.0, power);
      const Spline spline = {.Degree  = 2,
                             .Knots   = {0, 0, 0, 1, 1, 1},
                             .Poles   = {{.X = 9, .Y = -12}, {.X = 25, .Y = 0}, {.X = 9, .Y = 12}},
                             .Weights = {25, 3 * heavy, heavy * heavy}};
      CHECK(Area(spline, {}) == doctest::Approx(segment).epsilon(4e-16));
      CHECK(Near(Centroid(spline, {}), {.X = reach, .Y = 0}, 1e-14));
    }
  }

  TEST_CASE("area: arcs, ellipses and bulges in closed forms that do not cancel")
  {
    constexpr std::array<ArcCase, 12> Arcs = {{
        {.Sweep =  1.0, .Offset = 0.0, .Area =    0.713380568364465720063739552664},
        {.Sweep =  1.0, .Offset = 1e4, .Area =    0.713380568364465720063739552664},
        {.Sweep =  1.0, .Offset = 1e8, .Area =    0.713380568364465720063739552664},
        {.Sweep = 1e-2, .Offset = 0.0, .Area =  7.49996250008928605864932591112e-7},
        {.Sweep = 1e-2, .Offset = 1e4, .Area =  7.49996250008928605864932591112e-7},
        {.Sweep = 1e-2, .Offset = 1e8, .Area =  7.49996250008928605864932591112e-7},
        {.Sweep = 1e-4, .Offset = 0.0, .Area = 7.49999999625000107913191678103e-13},
        {.Sweep = 1e-4, .Offset = 1e4, .Area = 7.49999999625000107913191678103e-13},
        {.Sweep = 1e-4, .Offset = 1e8, .Area = 7.49999999625000107913191678103e-13},
        {.Sweep = 1e-6, .Offset = 0.0, .Area = 7.49999999999962398183251609145e-19},
        {.Sweep = 1e-6, .Offset = 1e4, .Area = 7.49999999999962398183251609145e-19},
        {.Sweep = 1e-6, .Offset = 1e8, .Area = 7.49999999999962398183251609145e-19},
    }};
    for (const ArcCase &row : Arcs)
    {
      CAPTURE(row.Sweep);
      CAPTURE(row.Offset);
      const Arc arc = {.Center = {.X = row.Offset, .Y = row.Offset},
                       .Radius = 3,
                       .Start  = Radians(0.1),
                       .Sweep  = Radians(row.Sweep)};
      CHECK(Steps(Area(arc), row.Area) <= 1);
      const Arc back = {.Center = arc.Center,
                        .Radius = 3,
                        .Start  = Radians(0.1 + row.Sweep),
                        .Sweep  = Radians(-row.Sweep)};
      CHECK(Area(back) == -Area(arc));
    }

    /**
     *  An ellipse of 1e6 by 1 turned by 0.2, swept 2.5 and 1e-3.
     */
    const Matrix flat = {.X      = {.X = 1e6 * std::cos(0.2), .Y = 1e6 * std::sin(0.2)},
                         .Y      = {.X = -std::sin(0.2), .Y = std::cos(0.2)},
                         .Origin = {}};
    CHECK(Steps(Area(Ellipse{.Map = flat, .Start = Radians(0.1), .Sweep = Radians(2.5)}),
                950763.927948021731165882300543)
          <= 1);
    CHECK(Steps(Area(Ellipse{.Map = flat, .Start = Radians(0.1), .Sweep = Radians(1e-3)}),
                8.33333291666667691657225045598e-5)
          <= 1);
    CHECK(Area(Ellipse{.Map = flat, .Start = Radians(2.6), .Sweep = Radians(-2.5)})
          == -Area(Ellipse{.Map = flat, .Start = Radians(0.1), .Sweep = Radians(2.5)}));

    /**
     *  A bulge on a chord of length 1: its series up to 1/2, its sweep past it.
     */
    constexpr std::array<std::array<double, 2>, 11> Bulges = {{
        {1e-3, 3.3333339999999048313254885945e-4},
        {1e-5, 3.33333333340000027267589401551e-6},
        {0.25, 0.0843658878411229730916549771098},
        {0.5, 0.174724694531879778292387680829},
        {0.5000000000000001, 0.174724694531879820554218729268},
        {0.75, 0.276205028642189879992843548566},
        {1.0, 0.39269908169872415480783042291},
        {3.0, 2.06811912833090892476377371845},
        {1000.0, 196349.933215306725668928905169},
        {-0.3, -0.101777524332158179436830376437},
        {-2.0, -1.05245993577663320548208239076},
    }};
    for (const auto &[bulge, area] : Bulges)
    {
      CAPTURE(bulge);
      const Polyline polyline = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = bulge},
                                              {.Position = {.X = 1, .Y = 0}, .Bulge = 0}},
                                 .Closed   = false};
      CHECK(Steps(Area(polyline), area) <= 1);
    }

    /**
     *  The centroids, likewise from the center and from the chord: a far, thin arc and a bulge of
     *  1e-5 on a chord a million away.
     */
    const Arc thin = {
        .Center = {.X = 1e8, .Y = 1e8}, .Radius = 3, .Start = Radians(0.1), .Sweep = Radians(1e-4)};
    CHECK(Near(Centroid(thin),
               {.X = 100000002.984997514851572755702, .Y = 100000000.299649499966101552166}, 3e-8));
    const Polyline shallow = {.Vertices = {{.Position = {.X = 1e6, .Y = 0}, .Bulge = 1e-5},
                                           {.Position = {.X = 1e6 + 1, .Y = 0}, .Bulge = 0}},
                              .Closed   = false};
    const std::optional<Point> lobe = Centroid(shallow);
    REQUIRE(lobe.has_value());
    CHECK(lobe->X == 1e6 + 0.5);
    CHECK(lobe->Y == doctest::Approx(-0.0000020000000000171430207482411656).epsilon(1e-12));
    CHECK(!Centroid(Arc{.Center = {}, .Radius = 1, .Start = {}, .Sweep = {}}).has_value());
  }

  /**
   *  The Bézier area was called exact and was 210 steps off at a million: the closed form, its
   *  offsets exact and its sum in double-double, is the double nearest the true one.
   */
  TEST_CASE("area: Bézier curves in closed form, the double nearest")
  {
    const Bezier::Cubic hump = {.Start = {.X = 0, .Y = 0},
                                .Out   = {.X = 1, .Y = 3},
                                .In    = {.X = 3, .Y = 2},
                                .End   = {.X = 4, .Y = 0}};
    CHECK(Area(hump) == -5.25);
    const Bezier::Cubic far = {.Start = {.X = 1e6 + 0.5, .Y = 1e6 + 0.25},
                               .Out   = {.X = 1e6 + 3.125, .Y = 1e6 - 1.5},
                               .In    = {.X = 1e6 - 2.0625, .Y = 1e6 + 4.75},
                               .End   = {.X = 1e6 + 1.0, .Y = 1e6}};
    CHECK(Area(far) == 831.0 / 1280.0);
    const Bezier::Cubic reversed = {
        .Start = far.End, .Out = far.In, .In = far.Out, .End = far.Start};
    CHECK(Area(reversed) == -Area(far));
    const Bezier::Quadratic arch = {.Start   = {.X = 1e6, .Y = 0},
                                    .Control = {.X = 1e6 + 1, .Y = 2},
                                    .End     = {.X = 1e6 + 2, .Y = 0}};
    CHECK(Area(arch) == -4.0 / 3.0);
    const Bezier::Cubic dot = {.Start = {.X = 1, .Y = 2},
                               .Out   = {.X = 1, .Y = 2},
                               .In    = {.X = 1, .Y = 2},
                               .End   = {.X = 1, .Y = 2}};
    CHECK(Area(dot) == 0);
    const Bezier::Cubic large = {.Start = {.X = 0, .Y = 0},
                                 .Out   = {.X = 0x1p500, .Y = 0x1.8p501},
                                 .In    = {.X = 0x1.8p501, .Y = 0x1p501},
                                 .End   = {.X = 0x1p502, .Y = 0}};
    CHECK(Area(large) == -5.25 * 0x1p1000);
  }

  /**
   *  A polynomial spline's area was its Bézier pieces' worked out in doubles, off by up to ten
   *  steps of a double on the comparison's splines: rS9 by 9.6, S1 by 1.03. Each piece's control
   *  points are now worked out in double-double with their bounds, and the area in closed form from
   *  them: the double nearest, to the bit, as mpmath's 30 digits round.
   */
  TEST_CASE("area: polynomial splines, the double nearest")
  {
    const Spline bent     = {.Degree  = 3,
                             .Knots   = {0, 0, 0, 0, 0.7, 1.9, 3, 3, 3, 3},
                             .Poles   = {{.X = 0, .Y = 0},
                                         {.X =   1, .Y =  3},
                                         {.X = 2.5, .Y = -1},
                                         {.X =   4, .Y =  4},
                                         {.X =   6, .Y =  0},
                                         {.X = 7, .Y = 2}},
                             .Weights = {}};
    const Spline arched   = {.Degree  = 2,
                             .Knots   = {0.0, 0.0, 0.0, 0.3458243800811477, 1.0, 1.0, 1.0},
                             .Poles   = {{.X = 0.9106232176529652, .Y = 0.5657448127478695},
                                         {.X =   3.952008389128002, .Y =  -4.946252033066458},
                                         {.X = -3.1281197735486446, .Y = -0.5191628686621232},
                                         {.X = 4.43837469275355, .Y = -4.638525337169864}},
                             .Weights = {}};
    const Spline reckoned = {
        .Degree  = 3,
        .Knots   = {0.0, 0.0, 0.0, 0.0, 0.2003898059481518, 0.7057395194767506, 1.0, 1.0, 1.0, 1.0},
        .Poles   = {{.X = 2.3259675147015946, .Y = -3.7476680729407494},
                    {.X = 1.6638812112063253, .Y =  0.49737073446468827},
                    {.X = 3.2875272789740055, .Y = -0.19561704218467746},
                    {.X = 3.2610394724386857, .Y =   -1.737134402007403},
                    {.X = -0.543801421583229, .Y =  -2.6114327783979463},
                    {.X = -4.643496717099503, .Y = -0.7346559977535598}},
        .Weights = {}};
    CHECK(Area(bent, {})     == -0x1.625eb187e4a31p+1);
    CHECK(Area(reckoned, {}) == 0x1.bbad6588f7cdbp+0);
    CHECK(Area(AreaReversed(bent), {}) == 0x1.625eb187e4a31p+1);
    CHECK(Area(arched, {})             == 0x1.96875a9753820p+0);
    CHECK(Area(arched, {})             == -Area(AreaReversed(arched), {}));

    /**
     *  A spline of one span, clamped, is its Bézier curve, whose area is the double nearest in
     *  closed form: the two to the bit.
     */
    std::uint64_t state = 0x5EED5U;
    const auto    draw  = [&state]() -> double
    {
      state = (state * 6364136223846793005U) + 1442695040888963407U;
      return (static_cast<double>(state >> 11U) * 0x1p-53 * 8.0) - 4.0;
    };
    for (int round = 0; round < 2000; ++round)
    {
      const std::array<Point, 4> points    = {{{.X = draw(), .Y = draw()},
                                               {.X = draw(), .Y = draw()},
                                               {.X = draw(), .Y = draw()},
                                               {.X = draw(), .Y = draw()}}};
      const Spline               quadratic = {.Degree  = 2,
                                              .Knots   = {0, 0, 0, 1, 1, 1},
                                              .Poles   = {points[0], points[1], points[2]},
                                              .Weights = {}};
      const Spline               cubic     = {.Degree  = 3,
                                              .Knots   = {0, 0, 0, 0, 1, 1, 1, 1},
                                              .Poles   = {points.begin(), points.end()},
                                              .Weights = {}};
      CAPTURE(round);
      CHECK(Area(quadratic, {})
            == Area(Bezier::Quadratic{.Start = points[0], .Control = points[1], .End = points[2]}));
      CHECK(Area(cubic, {})
            == Area(Bezier::Cubic{
                .Start = points[0], .Out = points[1], .In = points[2], .End = points[3]}));
    }

    /**
     *  A spline that runs out along a line and back encloses nothing, which no bound settles: the
     *  pieces' sum decides it. One past the highest degree worked out so is summed as before.
     */
    const Spline back = {
        .Degree  = 2,
        .Knots   = {0, 0, 0, 1, 2, 2, 2},
        .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 1}, {.X = 2, .Y = 2}, {.X = 0, .Y = 0}},
        .Weights = {}};
    CHECK(Area(back, {}) == 0.0);
    Spline high = {.Degree = 27, .Knots = {}, .Poles = {}, .Weights = {}};
    for (int index = 0; index <= 27; ++index)
    {
      const double angle = std::numbers::pi * index / 27.0;
      high.Poles.push_back({.X = std::cos(angle), .Y = std::sin(angle)});
    }
    high.Knots.assign(28, 0.0);
    high.Knots.insert(high.Knots.end(), 28, 1.0);
    REQUIRE(Valid(high));
    CHECK(Area(high, {}) > 1.0);
    CHECK(Area(AreaReversed(high), {}) == -Area(high, {}));
  }

  TEST_CASE("area: any curve as its own kind")
  {
    CHECK(Area(Curve{Circle({}, 1)}, {}) == doctest::Approx(std::numbers::pi).epsilon(1e-15));
    CHECK(Area(Curve{Tests::Quarter()}, {.Error = 1e-13})
          == doctest::Approx((std::numbers::pi / 4) - 0.5).epsilon(1e-12));
    CHECK(Near(Centroid(Curve{Circle({.X = 1, .Y = 1}, 1)}, {}), {.X = 1, .Y = 1}, 1e-15));
    CHECK(!Centroid(Curve{Tests::Clamped()}, {}).has_value());
  }
} // namespace Geoxy
