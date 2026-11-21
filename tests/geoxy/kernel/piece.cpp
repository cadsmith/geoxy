/***************************************************************************************************
 *  @file      tests/geoxy/kernel/piece.cpp
 *  @brief     Tests of the rational Bézier kernel's pieces.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "kernel/piece.hpp"

#include "../curves.hpp"
#include "basis.hpp"
#include "scratch.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/box.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <span>
#include <vector>

namespace Geoxy::Kernel
{
  namespace
  {
    /**
     *  The point of a piece at a parameter, by de Casteljau.
     *
     *  @param  piece The piece.
     *  @param  along Where.
     *  @return       The point.
     */
    auto Evaluated(const Piece &piece, double along) -> Point
    {
      return End(Split(piece, along)[0]);
    }
  } // namespace

  TEST_CASE("piece: a polynomial piece is its curve, and its halves meet on it")
  {
    const Piece piece = Polynomial(Tests::Wave);
    CHECK(Start(piece) == Tests::Wave.Start);
    CHECK(End(piece)   == Tests::Wave.End);
    for (const double along : {0.25, 0.5, 0.75})
    {
      CHECK(Distance(Evaluated(piece, along), At(Tests::Wave, along)) <= 1e-14);
    }
    const std::array<Piece, 2> halves = Split(piece, 0.3);
    CHECK(End(halves[0])   == Start(halves[1]));
    CHECK(Start(halves[0]) == Start(piece));
    CHECK(End(halves[1])   == End(piece));
  }

  TEST_CASE("piece: a line is of degree 1")
  {
    const Piece piece = Straight({.X = 1, .Y = 2}, {.X = 3, .Y = 4});
    CHECK(piece.Degree == 1);
    CHECK(Start(piece) == Point{.X = 1, .Y = 2});
    CHECK(End(piece)   == Point{.X = 3, .Y = 4});
    CHECK(Flat(piece, 0.0));
  }

  TEST_CASE("piece: a spline's span is the Bézier curve of its blossom")
  {
    const Piece piece = Section(Tests::Clamped(), 3);
    CHECK(Distance(Projected(Poles(piece)[1]), Tests::Wave.Out) <= 1e-15);
    CHECK(Distance(Projected(Poles(piece)[2]), Tests::Wave.In) <= 1e-15);
    std::vector<Point> ends;
    Sections(Tests::Uneven(),
             [&ends](const Piece &each, Range /*span*/) -> void { ends.push_back(End(each)); });
    CHECK(ends.size() == 3);
    CHECK(Distance(ends[0], At(Tests::Uneven(), 0.3)) <= 1e-14);
    CHECK(Distance(ends[1], At(Tests::Uneven(), 1.1)) <= 1e-14);
    const Piece quarter = Section(Tests::Quarter(), 2);
    CHECK(Distance(Evaluated(quarter, 0.5), {}) == doctest::Approx(1.0).epsilon(1e-15));
  }

  TEST_CASE("piece: the hull holds the piece, and flatness is measured against the chord")
  {
    const Piece piece = Polynomial(Tests::Wave);
    const Box   hull  = Hull(piece);
    for (const double along : {0.1, 0.4, 0.6, 0.9})
    {
      CHECK(Contains(hull, Evaluated(piece, along)));
    }
    CHECK_FALSE(Flat(piece, 1.0));
    CHECK(Flat(piece, 10.0));
    const std::array loop = {Point{.X = 0, .Y = 0}, Point{.X = 1, .Y = 1}, Point{.X = 0, .Y = 0}};
    CHECK_FALSE(Flat(Polynomial(loop), 1.0));
    CHECK(Flat(Polynomial(loop), 1.5));
  }

  TEST_CASE(
      "piece: bound grows a box to hold a piece, within a few steps, or whole past its budget")
  {
    const Piece piece = Polynomial(Tests::Wave);
    const Box   ends  = Grown(Grown(Box{}, Start(piece)), End(piece));
    const Box   tight = Bound(piece, ends, Budget);
    const Box   loose = Bound(piece, ends, 0);
    CHECK(Contains(Padded(Hull(piece), 0.0), tight.Max));
    CHECK(loose == Hull(piece));
    CHECK(tight.Max.Y < loose.Max.Y);
    CHECK(Bound(piece, Hull(piece), Budget) == Hull(piece));
    const Box above = {.Min = {.X = -1, .Y = -1}, .Max = {.X = 5, .Y = 4}};
    CHECK(Bound(piece, above, Budget) == above);
    Box nearly   = Hull(piece);
    nearly.Max.Y = std::nextafter(nearly.Max.Y, 0.0);
    CHECK(Bound(piece, nearly, Budget) == Hull(piece));
  }

  TEST_CASE("piece: roots are where a piece crosses a level, to the last bit")
  {
    const Piece line = Straight({.X = 0, .Y = 0}, {.X = 4, .Y = 1});
    CHECK(Roots(line, {.Along = Axis::X, .Value = 1}, Budget).Roots   == List<double>{0.25});
    CHECK(Roots(line, {.Along = Axis::Y, .Value = 0.5}, Budget).Roots == List<double>{0.5});
    CHECK(Roots(line, {.Along = Axis::X, .Value = 5}, Budget).Roots.empty());
    CHECK(Roots(line, {.Along = Axis::X, .Value = 0}, Budget).Roots == List<double>{0.0});
    CHECK(Roots(line, {.Along = Axis::X, .Value = 4}, Budget).Roots == List<double>{1.0});

    const Piece upright = Straight({.X = 1, .Y = 0}, {.X = 1, .Y = 5});
    CHECK(Roots(upright, {.Along = Axis::X, .Value = 1}, Budget).Roots == List<double>{0.0, 1.0});

    const Level        half  = {.Along = Axis::Y, .Value = 0.5};
    const List<double> roots = Roots(Polynomial(Tests::Wave), half, Budget).Roots;
    REQUIRE(roots.size() == 3);
    for (const double root : roots)
    {
      CHECK(std::abs(At(Tests::Wave, root).Y - 0.5) < 1e-15);
    }
    CHECK(roots[0] < roots[1]);
    CHECK(roots[1] < roots[2]);
  }

  TEST_CASE("piece: a root a halving lands on is kept, and a spent budget keeps one of several")
  {
    const std::array odd   = {Point{.X = 0, .Y = -1}, Point{.X = 1, .Y = 1}, Point{.X = 2, .Y = -1},
                              Point{.X = 3, .Y = 1}};
    const Level      level = {.Along = Axis::Y, .Value = 0};
    CHECK(Roots(Polynomial(odd), level, Budget).Roots == List<double>{0.5});

    const Level half = {.Along = Axis::Y, .Value = 0.5};
    CHECK(Roots(Polynomial(Tests::Wave), half, 0).Roots.size() == 1);

    const std::array arch = {Point{.X = 0, .Y = 0}, Point{.X = 1, .Y = 2}, Point{.X = 2, .Y = 0}};
    CHECK(Roots(Polynomial(arch), half, Budget).Roots.size() == 2);
    CHECK(Roots(Polynomial(arch), half, 0).Roots.empty());
  }

  TEST_CASE("piece: a chain's seams hold the roots its pieces miss, a closed chain's last too")
  {
    const std::vector<Solution> open  = {{.Roots = {0.5}, .Start = -1, .End = 1e-17},
                                         {.Roots = {}, .Start = -1e-17, .End = -1},
                                         {.Roots = {}, .Start = -1, .End = 1}};
    const List<Spot>            spots = Chained(open, false);
    REQUIRE(spots.size()             == 2);
    CHECK(spots[0].Index             == 0);
    CHECK(spots[0].Along             == doctest::Approx(0.5));
    CHECK(spots[1].Index             == 1);
    CHECK(spots[1].Along             == doctest::Approx(0.0));
    CHECK(Chained(open, true).size() == 3);
    CHECK(Chained(open, true).back().Index == 0);
  }

  TEST_CASE("piece: a chain is looked across a piece that found nothing and tells nothing")
  {
    /**
     *  A piece of no length between two others, as a vertex given twice: the function turns across
     *  it, and the seam is the start of the piece after. Before, neither side saw the turn.
     */
    const std::vector<Solution> open  = {{.Roots = {}, .Start = -1, .End = -2},
                                         {.Roots = {}, .Start = 0, .End = 0},
                                         {.Roots = {}, .Start = 0, .End = 0},
                                         {.Roots = {}, .Start = 3, .End = 1}};
    const List<Spot>            spots = Chained(open, false);
    REQUIRE(spots.size() == 1);
    CHECK(spots[0].Index == 3);
    CHECK(spots[0].Along == 0.0);

    const std::vector<Solution> same = {{.Roots = {}, .Start = -1, .End = -2},
                                        {.Roots = {}, .Start = 0, .End = 0},
                                        {.Roots = {}, .Start = -3, .End = -1}};
    CHECK(Chained(same, false).empty());

    /**
     *  Round a closed chain the piece after the last may be past pieces of no length at its start,
     *  and a chain of nothing else has no seam at all.
     */
    const std::vector<Solution> round = {{.Roots = {}, .Start = 0, .End = 0},
                                         {.Roots = {}, .Start = 1, .End = 2},
                                         {.Roots = {}, .Start = -2, .End = -1},
                                         {.Roots = {}, .Start = 0, .End = 0}};
    const List<Spot>            wound = Chained(round, true);
    REQUIRE(wound.size() == 2);
    CHECK(wound[0].Index == 2);
    CHECK(wound[1].Index == 1);
    CHECK(Chained(round, false).size() == 1);
    const std::vector<Solution> blank = {{.Roots = {}, .Start = 0, .End = 0},
                                         {.Roots = {}, .Start = 0, .End = 0}};
    CHECK(Chained(blank, true).empty());

    /**
     *  A piece that tells of one end is not looked across.
     */
    const std::vector<Solution> ended = {{.Roots = {}, .Start = -1, .End = -2},
                                         {.Roots = {}, .Start = 0, .End = -1},
                                         {.Roots = {}, .Start = 3, .End = 1}};
    REQUIRE(Chained(ended, false).size() == 1);
    CHECK(Chained(ended, false)[0].Index == 2);
  }

  TEST_CASE("piece: signs are told apart where their product would come to nothing")
  {
    /**
     *  Two values of 2^-600 and opposite signs: their product, 2^-1200, is zero to a double, and
     *  the root between them and the seam across them were both let go.
     */
    const std::array<double, 2> tiny = {-0x1p-600, 0x1p-600};
    const List<double>          root = Roots(tiny, Budget);
    REQUIRE(root.size() == 1);
    CHECK(root[0]       == doctest::Approx(0.5));
    const std::array<double, 4> spent = {-0x1p-600, 0x1p-600, -0x1p-600, 0x1p-600};
    CHECK(Roots(spent, 0).size() == 1);
    const std::vector<Solution> seam = {{.Roots = {}, .Start = -0x1p-600, .End = -0x1p-600},
                                        {.Roots = {}, .Start = 0x1p-600, .End = 0x1p-600}};
    CHECK(Chained(seam, false).size() == 1);
  }

  TEST_CASE("piece: a derivative's numerator, pair by pair, keeps a light weight's end")
  {
    /**
     *  A quadratic whose middle weight is 2^-100 of its ends': it leaves its start toward the
     *  middle pole, at a rate 2^-100 of the rest. From its poles' slopes the rate at the end was
     *  rounding, and the end came out a foot of every point.
     */
    const Piece     light = {.Degree = 2,
                             .Poles  = {Homogeneous{.X = 0, .Y = 0, .W = 1},
                                        Homogeneous{.X = 0x1p-100, .Y = 0x1p-100, .W = 0x1p-100},
                                        Homogeneous{.X = 2, .Y = 0, .W = 1}}};
    const Hodograph speed = Derived(light);
    REQUIRE(speed.Degree  == 2);
    CHECK(speed.Across[0] == 0x1p-99);
    CHECK(speed.Rising[0] == 0x1p-99);
    CHECK(speed.Across[1] == 2);
    CHECK(speed.Rising[1] == 0);
    CHECK(speed.Across[2] == 0x1p-99);
    CHECK(speed.Rising[2] == -0x1p-99);

    const Point    above = {.X = 0.7, .Y = 0.9};
    const Solution feet  = Perpendiculars(light, above, Budget);
    REQUIRE(feet.Roots.size() == 1);
    CHECK(feet.Roots[0] > 0.3);
    CHECK(feet.Roots[0] < 0.6);
    CHECK(feet.Start < 0.0);
    CHECK(feet.End > 0.0);
    const Solution trod = Perpendiculars(Tread{.Shape = light, .Speed = speed}, above, Budget);
    REQUIRE(trod.Roots.size() == 1);
    CHECK(trod.Roots[0]       == doctest::Approx(feet.Roots[0]).epsilon(1e-12));

    /**
     *  Seen from above the piece turns one way about the point all along but for its last sliver,
     *  where it bends from the middle pole's way: one touch there, and none seen from below.
     */
    const Solution touches = Touches(light, above, Budget);
    REQUIRE(touches.Roots.size() == 1);
    CHECK(touches.Roots[0] > 0.999);
    CHECK(touches.Start > 0.0);
    CHECK(touches.End < 0.0);
    const Solution touched = Touches(Tread{.Shape = light, .Speed = speed}, above, Budget);
    REQUIRE(touched.Roots.size() == 1);
    CHECK(touched.Roots[0] > 0.999);
    CHECK(Touches(light, {.X = 1, .Y = -1}, Budget).Roots.empty());

    /**
     *  Halved, the numerator is the same polynomial over each half, and a span's, drawn from the
     *  spline's own poles, is its piece's.
     */
    const std::array<Hodograph, 2> halves = Halved(speed);
    CHECK(halves[0].Across[0] == speed.Across[0]);
    CHECK(halves[1].Across[2] == speed.Across[2]);
    CHECK(halves[0].Across[2] == halves[1].Across[0]);
    CHECK(halves[0].Rising[2] == halves[1].Rising[0]);
    const Spline    uneven = Tests::Uneven();
    const Hodograph whole  = Derived(uneven, 6);
    const Hodograph pieced = Derived(Section(uneven, 6));
    REQUIRE(whole.Degree == pieced.Degree);
    const double scale = whole.Across[0] / pieced.Across[0];
    CHECK(scale > 0.0);
    for (std::size_t index = 0; index <= whole.Degree; ++index)
    {
      CHECK(whole.Across.at(index) == doctest::Approx(scale * pieced.Across.at(index)));
      CHECK(whole.Rising.at(index) == doctest::Approx(scale * pieced.Rising.at(index)));
    }
    Spline plain  = uneven;
    plain.Weights = {};
    CHECK(Derived(plain, 6).Degree == 4);
    CHECK(Stepped(Spline{}, [](std::size_t /*span*/) -> Solution { return {}; }).empty());

    /**
     *  A span's touches, from its poles' own pairs, are its piece's.
     */
    const Point    aside  = {.X = 3, .Y = 3};
    const Solution polar  = Touches(uneven, 6, aside);
    const Solution homely = Touches(Section(uneven, 6), aside, Budget);
    REQUIRE(polar.Roots.size() == homely.Roots.size());
    for (std::size_t index = 0; index < polar.Roots.size(); ++index)
    {
      CHECK(polar.Roots[index] == doctest::Approx(homely.Roots[index]).epsilon(1e-12));
    }
    CHECK((polar.Start > 0.0) == (homely.Start > 0.0));
    CHECK((polar.End > 0.0)   == (homely.End > 0.0));

    /**
     *  A numerator past every double, as of control points near the largest of them, is left as it
     *  is: nothing is scaled by a power no number has.
     */
    Hodograph beyond = speed;
    beyond.Across[1] = std::numeric_limits<double>::infinity();
    CHECK(Perpendiculars(Tread{.Shape = light, .Speed = beyond}, above, Budget).Roots.size() <= 2);
  }

  TEST_CASE("piece: a quarter of an ellipse is an exact conic, its parameter turned to an angle")
  {
    const Matrix map     = Rotation(Degrees(30)) * Scaling(3, 1);
    const Piece  quarter = Quarter(map, 1);
    for (const double along : {0.0, 0.2, 0.5, 0.9, 1.0})
    {
      const double angle = Quartered(1, along);
      const Point  point = map * (Point{} + Direction(Radians(angle)));
      CHECK(Distance(Evaluated(quarter, along), point) < 1e-15);
    }
    CHECK(Quartered(0, 0.5) == doctest::Approx(std::numbers::pi / 4));
    CHECK(Quartered(2, 1.0) == Quartered(3, 0.0));
  }

  TEST_CASE("piece: feet of a piece are where its tangent is square to the point")
  {
    const Piece line = Straight({.X = 0, .Y = 0}, {.X = 4, .Y = 0});
    CHECK(Perpendiculars(line, {.X = 1, .Y = 2}, Budget).Roots == List<double>{0.25});
    CHECK(Perpendiculars(line, {.X = 1, .Y = 2}, Budget).Start < 0.0);
    CHECK(Perpendiculars(line, {.X = 1, .Y = 2}, Budget).End > 0.0);

    /**
     *  One weight throughout, negative: the same feet, by the full rational numerator.
     */
    Piece backward = line;
    for (Homogeneous &pole : std::span{backward.Poles}.first(backward.Degree + 1))
    {
      pole = {.X = -pole.X, .Y = -pole.Y, .W = -pole.W};
    }
    CHECK(Perpendiculars(backward, {.X = 1, .Y = 2}, Budget).Roots == List<double>{0.25});

    const Piece circle = Quarter({}, 0);
    CHECK(Perpendiculars(circle, {}, Budget).Roots.empty());
    CHECK(Perpendiculars(circle, {}, Budget).Start == 0.0);

    const Range span = {.Start = 1, .End = 3};
    CHECK(Across(span, 0.5) == 2);
    CHECK(Across(span, 1.0) == 3);
  }

  TEST_CASE("piece: touches are where a piece's tangent passes through the point")
  {
    const Piece line = Straight({.X = 0, .Y = 0}, {.X = 4, .Y = 0});
    CHECK(Touches(line, {.X = 1, .Y = 2}, Budget).Roots.empty());
    CHECK(Touches(line, {.X = 1, .Y = 2}, Budget).Start > 0.0);
    CHECK(Touches(line, {.X = 9, .Y = 0}, Budget).End == 0.0);

    const Piece arch = Polynomial(Bezier::Quadratic{
        .Start = {.X = 0, .Y = 0}, .Control = {.X = 1, .Y = 2}, .End = {.X = 2, .Y = 0}});
    CHECK(Touches(arch, {.X = 1, .Y = 1.5}, Budget).Roots.size() == 2);
  }
} // namespace Geoxy::Kernel
