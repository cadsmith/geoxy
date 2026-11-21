/***************************************************************************************************
 *  @file      tests/geoxy/fit.cpp
 *  @brief     Tests of splines fit through points, against mpmath solving the same system at 60
 *             digits, and of knots inserted, degrees raised, polylines framed and splines exploded.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/fit.hpp"

#include "curves.hpp"
#include "real.hpp"

#include "geoxy/at.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/error.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <expected>
#include <limits>
#include <optional>
#include <span>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  How many steps of a double at a scale two points are apart, the larger of their two
     *  coordinates' gaps.
     *
     *  @param  point The point.
     *  @param  near  The point it should be.
     *  @param  scale The scale whose steps count.
     *  @return       The steps.
     */
    auto Steps(Point point, Point near, double scale) -> double
    {
      return std::max(std::abs(point.X - near.X), std::abs(point.Y - near.Y)) / Spacing(scale);
    }

    /**
     *  Points whose chords are whole: 5, 13, 17, 25, 29 and 13.
     */
    constexpr std::array<Point, 7> Open = {{
        {.X =  0, .Y =  0},
        {.X =  3, .Y =  4},
        {.X =  8, .Y = 16},
        {.X = 16, .Y = 31},
        {.X = 23, .Y =  7},
        {.X = 43, .Y = 28},
        {.X = 31, .Y = 33},
    }};

    /**
     *  A loop whose chords are whole: 5, 13, 20, 15 and 11.
     */
    constexpr std::array<Point, 5> WholeChords = {{
        {.X =  0, .Y =   0},
        {.X =  3, .Y =   4},
        {.X = 15, .Y =   9},
        {.X = 15, .Y = -11},
        {.X =  0, .Y = -11},
    }};

    /**
     *  A fit's case: its points, its fitting, and the poles mpmath finds at 60 digits, rounded.
     */
    struct Solved
    {
      std::span<const Point> Points  = {}; /**< The points.  */
      Geoxy::Fitting         Fitting = {}; /**< The fitting. */
      std::vector<Point>     Poles   = {}; /**< The poles.   */
    };

    /**
     *  The cases worked out by the oracle, the poles as it printed them.
     *
     *  @return The cases.
     */
    auto SolvedCases() -> std::vector<Solved>
    {
      return {
          {.Points  = Open,
           .Fitting = {.Knots = Knotting::Chord},
           .Poles   = {{.X = 0.0, .Y = 0.0},
                       {.X = 1.0649960771153089, .Y =  1.3292862164045751},
                       {.X =  4.898981954730421, .Y =   6.114716595461045},
                       {.X =  7.440251064620232, .Y =   15.72568215010719},
                       {.X =  19.73424489722888, .Y =  44.230864453135126},
                       {.X = 15.524900495058391, .Y = -12.710099010962226},
                       {.X =  53.45176655931281, .Y =  29.084567928225628},
                       {.X =  36.30678118674666, .Y =    32.0745342375806},
                       {.X = 31.0, .Y = 33.0}}},
          {.Points  = Open,
           .Fitting = {.Knots = Knotting::Root},
           .Poles   = {{.X = 0.0, .Y = 0.0},
                       {.X = 1.0302697342692666, .Y =  1.2096905532515003},
                       {.X =   3.72179949812851, .Y =   4.369948513702094},
                       {.X = 7.2558497090682845, .Y =  14.079122926129807},
                       {.X = 18.790397576703004, .Y =   44.83334604070585},
                       {.X =  16.39874497410347, .Y = -11.321823776703152},
                       {.X =   52.0888058909492, .Y =   32.34535138377823},
                       {.X =  37.03645261065815, .Y =   32.81261360321251},
                       {.X = 31.0, .Y = 33.0}}},
          {.Points  = Open,
           .Fitting = {.Knots = Knotting::Uniform,
                       .Start = Vector{.X = 1, .Y = 0},
                       .End   = Vector{.X = 0, .Y = 1}},
           .Poles   = {{.X = 0.0, .Y = 0.0},
                       {.X = 1.6666666666666667, .Y =                 0.0},
                       {.X = 2.4854700854700855, .Y =   3.562393162393162},
                       {.X =  6.800854700854701, .Y =  11.531623931623932},
                       {.X =  18.31111111111111, .Y =   46.31111111111111},
                       {.X = 15.954700854700855, .Y = -10.776068376068377},
                       {.X =  55.87008547008547, .Y =   38.79316239316239},
                       {.X =               31.0, .Y =  28.666666666666668},
                       {.X = 31.0, .Y = 33.0}}},
          {.Points  = Open,
           .Fitting = {.Knots = Knotting::Chord,
                       .Start = Vector{.X = 1, .Y = 1},
                       .End   = Vector{.X = -1, .Y = 2}},
           .Poles   = {{.X = 0.0, .Y = 0.0},
                       {.X = 1.1785113019775793, .Y =  1.1785113019775793},
                       {.X = 4.7607165220295755, .Y =  6.3045288947146245},
                       {.X =  7.476757670456623, .Y =  15.609054216883026},
                       {.X = 19.883979747443597, .Y =   44.42943443173881},
                       {.X = 14.445726289631255, .Y = -13.680484984131418},
                       {.X =  57.05925766055522, .Y =  32.247221887478055},
                       {.X = 32.937925580499815, .Y =  29.124148839000366},
                       {.X = 31.0, .Y = 33.0}}},
          {.Points  = WholeChords,
           .Fitting = {.Knots = Knotting::Chord, .Closed = true},
           .Poles   = {{.X = -2.8123023264097076, .Y = -14.73933171443933},
                       {.X = -1.4113229870991906, .Y = -1.1884856132641826},
                       {.X =   4.271476231091928, .Y =  5.9881447478512495},
                       {.X =   19.82173010596626, .Y =  15.162601245376656},
                       {.X =  21.119658455275836, .Y = -15.169028995860323},
                       {.X = -2.8123023264097076, .Y =  -14.73933171443933},
                       {.X = -1.4113229870991906, .Y = -1.1884856132641826},
                       {.X = 4.271476231091928, .Y = 5.9881447478512495}}},
          {.Points  = WholeChords,
           .Fitting = {.Knots = Knotting::Root, .Closed = true},
           .Poles   = {{.X = -4.049024768108056, .Y = -13.917836735372324},
                       {.X = -0.12315038665075263, .Y =   0.798221744014723},
                       {.X =    2.180030470900711, .Y =  3.6073777876486837},
                       {.X =   18.804324813732922, .Y =  16.563738608431283},
                       {.X =   20.189670370553245, .Y = -16.505683758407905},
                       {.X =   -4.049024768108056, .Y = -13.917836735372324},
                       {.X = -0.12315038665075263, .Y =   0.798221744014723},
                       {.X = 2.180030470900711, .Y = 3.6073777876486837}}},
      };
    }

    /**
     *  Points moved by an offset.
     *
     *  @param  points The points.
     *  @param  offset The offset.
     *  @return        The points moved.
     */
    auto Moved(std::span<const Point> points, Vector offset) -> std::vector<Point>
    {
      std::vector<Point> moved;
      for (const Point point : points)
      {
        moved.push_back(point + offset);
      }
      return moved;
    }

    /**
     *  The most steps of a double at a scale the spline is from its points at their knots.
     *
     *  @param  spline The spline.
     *  @param  points The points.
     *  @param  scale  The scale.
     *  @return        The steps.
     */
    auto Missed(const Spline &spline, std::span<const Point> points, double scale) -> double
    {
      double worst = 0.0;
      for (std::size_t index = 0; index < points.size(); ++index)
      {
        worst = std::max(worst, Steps(At(spline, spline.Knots[index + 3]), points[index], scale));
      }
      return worst;
    }

    /**
     *  The most steps of a double at a scale one spline is from another, at a thousand places
     *  across the first one's domain.
     *
     *  @param  spline The spline.
     *  @param  other  The other.
     *  @param  scale  The scale.
     *  @return        The steps.
     */
    auto Strayed(const Spline &spline, const Spline &other, double scale) -> double
    {
      const Range domain = Domain(spline);
      double      worst  = 0.0;
      for (int index = 0; index <= 1000; ++index)
      {
        const double along = domain.Start + ((domain.End - domain.Start) * index / 1000.0);
        worst              = std::max(worst, Steps(At(other, along), At(spline, along), scale));
      }
      return worst;
    }

    /**
     *  A spline of uniform knots running past its domain at both ends.
     *
     *  @return The spline.
     */
    auto Unclamped() -> Spline
    {
      return {.Degree  = 3,
              .Knots   = {0, 1, 2, 3, 4, 5, 6, 7, 8},
              .Poles   = {{.X = 0, .Y = 0},
                          {.X = 1, .Y = 3},
                          {.X = 3, .Y = 4},
                          {.X = 5, .Y = 1},
                          {.X = 7, .Y = 2}},
              .Weights = {}};
    }
  } // namespace

  TEST_CASE("fit: poles as mpmath solves the same system")
  {
    for (const Vector offset : {Vector{}, Vector{.X = 5e6, .Y = -5e6}, Vector{.X = 1e8, .Y = 1e8}})
    {
      for (const Solved &solved : SolvedCases())
      {
        const std::vector<Point>           points = Moved(solved.Points, offset);
        const std::expected<Spline, Error> spline = Fit(points, solved.Fitting);
        REQUIRE(spline.has_value());
        REQUIRE(spline->Poles.size() == solved.Poles.size());
        const double scale = std::max(64.0, std::abs(offset.X));
        for (std::size_t index = 0; index < solved.Poles.size(); ++index)
        {
          CHECK(Steps(spline->Poles[index], solved.Poles[index] + offset, scale) <= 4);
        }
        CHECK(Missed(*spline, points, scale) <= 4);
      }
    }
  }

  TEST_CASE("fit: as many poles and knots as AutoCAD's, the ends exact")
  {
    const std::span<const Point>       six    = std::span{Open}.first(6);
    const std::expected<Spline, Error> spline = Fit(six, {});
    REQUIRE(spline.has_value());
    CHECK(spline->Degree        == 3);
    CHECK(spline->Poles.size()  == 8);
    CHECK(spline->Knots         == std::vector<double>{0, 0, 0, 0, 5, 18, 35, 60, 89, 89, 89, 89});
    CHECK(spline->Poles.front() == six.front());
    CHECK(spline->Poles.back()  == six.back());
    CHECK(At(*spline, 89)       == six.back());
    CHECK(Valid(*spline));

    const std::expected<Spline, Error> uniform = Fit(six, {.Knots = Knotting::Uniform});
    REQUIRE(uniform.has_value());
    CHECK(uniform->Knots == std::vector<double>{0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5});
  }

  TEST_CASE("fit: tangents leave at unit speed, natural ends do not bend")
  {
    const Fitting tangents = {.Start = Vector{.X = 0, .Y = 2}, .End = Vector{.X = 3, .Y = 0}};
    const std::expected<Spline, Error> leaving = Fit(Open, tangents);
    REQUIRE(leaving.has_value());
    CHECK(Steps(Point{} + Derivative(*leaving, 0), {.X = 0, .Y = 1}, 1) <= 4);
    CHECK(Steps(Point{} + Derivative(*leaving, leaving->Knots.back()), {.X = 1, .Y = 0}, 1) <= 4);

    const std::expected<Spline, Error> natural = Fit(Open, {});
    REQUIRE(natural.has_value());
    const std::vector<Point>  &poles  = natural->Poles;
    const std::vector<double> &knots  = natural->Knots;
    const Vector               first  = (poles[1] - poles[0]) / (knots[4] - knots[1]);
    const Vector               second = (poles[2] - poles[1]) / (knots[5] - knots[2]);
    CHECK(Steps(Point{} + first, Point{} + second, 1) <= 8);
    const std::size_t last = poles.size() - 1;
    const Vector      end  = (poles[last] - poles[last - 1]) / (knots[last + 3] - knots[last]);
    const Vector before = (poles[last - 1] - poles[last - 2]) / (knots[last + 2] - knots[last - 1]);
    CHECK(Steps(Point{} + end, Point{} + before, 1) <= 8);

    const std::array<Point, 2>         two  = {Point{.X = 1, .Y = 1}, Point{.X = 4, .Y = 5}};
    const std::expected<Spline, Error> line = Fit(two, {});
    REQUIRE(line.has_value());
    CHECK(Steps(At(*line, 2.5), {.X = 2.5, .Y = 3}, 4) <= 2);
  }

  TEST_CASE("fit: a closed fit runs on through its first point smoothly")
  {
    std::vector<Point> loop(WholeChords.begin(), WholeChords.end());
    loop.push_back(WholeChords.front());
    const std::expected<Spline, Error> spline = Fit(loop, {.Closed = true});
    REQUIRE(spline.has_value());
    CHECK(spline->Poles.size() == 8);
    CHECK(Valid(*spline));
    const Range domain = Domain(*spline);
    CHECK(domain == Range{.Start = 0, .End = 64});
    CHECK(Steps(At(*spline, domain.End), WholeChords.front(), 64) <= 4);
    CHECK(Steps(Point{} + Derivative(*spline, domain.Start),
                Point{} + Derivative(*spline, domain.End), 64)
          <= 8);

    const std::expected<Spline, Error> uniform =
        Fit(WholeChords, {.Knots = Knotting::Uniform, .Closed = true});
    REQUIRE(uniform.has_value());
    CHECK(Missed(*uniform, WholeChords, 64) <= 4);
  }

  TEST_CASE("fit: what cannot be fit")
  {
    const double             nan  = std::numeric_limits<double>::quiet_NaN();
    const double             huge = std::numeric_limits<double>::infinity();
    const std::vector<Point> one  = {{.X = 1, .Y = 1}};
    CHECK(Fit(one, {}).error() == Error::Degenerate);
    CHECK(Fit({}, {}).error()  == Error::Degenerate);
    const std::vector<Point> two = {{.X = 0, .Y = 0}, {.X = 1, .Y = 0}};
    CHECK(Fit(two, {.Closed = true}).error() == Error::Degenerate);
    const std::vector<Point> again = {{.X = 0, .Y = 0}, {.X = 1, .Y = 0}, {.X = 1, .Y = 0}};
    CHECK(Fit(again, {}).error() == Error::Degenerate);
    const std::vector<Point> near = {{.X = 0, .Y = 0}, {.X = 1e5, .Y = 0}, {.X = 1e5, .Y = 1e-12}};
    CHECK(Fit(near, {}).error() == Error::Degenerate);
    CHECK(Fit(two, {.Start = Vector{}}).error() == Error::Degenerate);
    CHECK(Fit(two, {.End = Vector{}}).error()   == Error::Degenerate);
    const std::vector<Point> lost = {{.X = 0, .Y = 0}, {.X = nan, .Y = 0}};
    CHECK(Fit(lost, {}).error() == Error::NotFinite);
    const std::vector<Point> high = {{.X = 0, .Y = 0}, {.X = 0, .Y = huge}};
    CHECK(Fit(high, {}).error() == Error::NotFinite);
    CHECK(Fit(two, {.Start = Vector{.X = huge, .Y = 0}}).error() == Error::NotFinite);
    CHECK(Fit(two, {.Start = Vector{.X = 0, .Y = nan}}).error()  == Error::NotFinite);
    CHECK(Fit(two, {.End = Vector{.X = nan, .Y = 0}}).error()    == Error::NotFinite);
    CHECK(Fit(two, {.End = Vector{.X = 0, .Y = huge}}).error()   == Error::NotFinite);
    const std::vector<Point> far = {{.X = -1.7e308, .Y = 0}, {.X = 1.7e308, .Y = 0}};
    CHECK(Fit(far, {}).error() == Error::NotFinite);
    const std::vector<Point> spread = {{.X = 0, .Y = 0}, {.X = 1.7e308, .Y = 0}, {.X = 0, .Y = 0}};
    CHECK(Fit(spread, {}).error() == Error::NotFinite);
    CHECK(Fit(spread, {.Knots = Knotting::Uniform}).error() == Error::NotFinite);
    CHECK(Fit(one, {.Closed = true}).error() == Error::Degenerate);
  }

  TEST_CASE("fit: knots inserted leave the curve as it was")
  {
    for (const Spline &spline : {Tests::Clamped(), Tests::Quarter(), Tests::Uneven(), Unclamped()})
    {
      const Range domain = Domain(spline);
      for (const double fraction : {0.0, 0.25, 0.5, 0.875, 1.0})
      {
        const double                knot = domain.Start + ((domain.End - domain.Start) * fraction);
        const std::optional<Spline> inserted = Inserted(spline, knot, 1);
        const bool clamped = std::ranges::count(spline.Knots, knot) > spline.Degree;
        REQUIRE(inserted.has_value() != clamped);
        if (inserted)
        {
          CHECK(Valid(*inserted));
          CHECK(inserted->Poles.size() == spline.Poles.size() + 1);
          CHECK(Domain(*inserted)      == domain);
          CHECK(Strayed(spline, *inserted, 8) <= 8);
        }
      }
    }
    const std::optional<Spline> twice = Inserted(Tests::Uneven(), 0.3, 1);
    REQUIRE(twice.has_value());
    CHECK(Strayed(Tests::Uneven(), *twice, 8) <= 4);
    const std::optional<Spline> thrice = Inserted(Unclamped(), 3.5, 3);
    REQUIRE(thrice.has_value());
    CHECK(Strayed(Unclamped(), *thrice, 8) <= 4);
    CHECK(Inserted(Unclamped(), 3.5, 0) == Unclamped());

    CHECK_FALSE(Inserted(Tests::Uneven(), 0.3, 2).has_value());
    CHECK_FALSE(Inserted(Tests::Clamped(), 0.5, 4).has_value());
    CHECK_FALSE(Inserted(Tests::Clamped(), 0.5, -1).has_value());
    CHECK_FALSE(Inserted(Tests::Clamped(), 1.5, 1).has_value());
    CHECK_FALSE(Inserted(Tests::Clamped(), -0.5, 1).has_value());
    CHECK_FALSE(
        Inserted(Tests::Clamped(), std::numeric_limits<double>::quiet_NaN(), 1).has_value());
    CHECK_FALSE(Inserted(Spline{}, 0.5, 1).has_value());
  }

  TEST_CASE("fit: degrees raised leave the curve as it was")
  {
    for (const Spline &spline : {Tests::Clamped(), Tests::Quarter(), Tests::Uneven(), Unclamped()})
    {
      for (const int times : {1, 2, 5})
      {
        const std::optional<Spline> raised = Elevated(spline, times);
        REQUIRE(raised.has_value());
        CHECK(raised->Degree == spline.Degree + times);
        CHECK(Valid(*raised));
        CHECK(Domain(*raised) == Domain(spline));
        CHECK(Strayed(spline, *raised, 8) <= 8);
      }
    }
    const std::optional<Spline> uneven = Elevated(Tests::Uneven(), 1);
    REQUIRE(uneven.has_value());
    CHECK(uneven->Knots
          == std::vector<double>{0, 0, 0, 0, 0, 0.3, 0.3, 0.3, 1.1, 1.1, 2, 2, 2, 2, 2});
    CHECK(Elevated(Tests::Clamped(), 0) == Tests::Clamped());
    CHECK_FALSE(Elevated(Tests::Clamped(), -1).has_value());
    CHECK_FALSE(Elevated(Tests::Clamped(), MaxDegree - 2).has_value());
    CHECK(Elevated(Tests::Clamped(), MaxDegree - 3).has_value());
    CHECK_FALSE(Elevated(Spline{}, 1).has_value());
  }

  TEST_CASE("fit: a polyline framing a spline, as PEDIT Spline makes it")
  {
    Polyline polyline = {};
    for (const Point point : std::span{Open}.first(5))
    {
      polyline.Vertices.push_back({.Position = point, .Bulge = 0.5});
    }
    const std::optional<Spline> cubic = Framed(polyline, 3);
    REQUIRE(cubic.has_value());
    CHECK(cubic->Knots        == std::vector<double>{0, 0, 0, 0, 1, 2, 2, 2, 2});
    CHECK(cubic->Poles.size() == 5);
    CHECK(Valid(*cubic));
    CHECK(At(*cubic, 2) == Open[4]);

    const std::optional<Spline> quadratic = Framed(polyline, 2);
    REQUIRE(quadratic.has_value());
    CHECK(quadratic->Knots == std::vector<double>{0, 0, 0, 1, 2, 3, 3, 3});

    polyline.Closed                  = true;
    const std::optional<Spline> loop = Framed(polyline, 3);
    REQUIRE(loop.has_value());
    CHECK(loop->Poles.size() == 8);
    CHECK(loop->Knots        == std::vector<double>{-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8});
    CHECK(Valid(*loop));
    CHECK(Steps(At(*loop, 0), At(*loop, 5), 32) <= 4);

    CHECK_FALSE(Framed(polyline, 0).has_value());
    CHECK_FALSE(Framed(polyline, MaxDegree + 1).has_value());
    CHECK_FALSE(Framed(polyline, 5).has_value());
    CHECK(Framed(polyline, 4).has_value());
  }

  TEST_CASE("fit: a spline exploded into chords")
  {
    const Polyline open = Explode(Tests::Clamped(), {.Error = 1e-3});
    CHECK_FALSE(open.Closed);
    REQUIRE(open.Vertices.size() > 8);
    CHECK(open.Vertices.front().Position == Tests::Wave.Start);
    CHECK(open.Vertices.back().Position  == Tests::Wave.End);

    Polyline                    square = {.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                                       {.Position = {.X = 4, .Y = 0}},
                                                       {.Position = {.X = 4, .Y = 4}},
                                                       {.Position = {.X = 0, .Y = 4}},
                                                       {.Position = {.X = 0, .Y = 0}}},
                                          .Closed   = false};
    const std::optional<Spline> framed = Framed(square, 2);
    REQUIRE(framed.has_value());
    const Polyline closed = Explode(*framed, {.Error = 1e-3});
    CHECK(closed.Closed);
    CHECK(closed.Vertices.front().Position      == Point{});
    CHECK_FALSE(closed.Vertices.back().Position == Point{});

    const Spline straight = {.Degree  = 1,
                             .Knots   = {0, 0, 1, 1},
                             .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 0}},
                             .Weights = {}};
    CHECK(Explode(straight, {.Error = 1e-3}).Vertices.size() == 2);
    CHECK(Explode(Spline{}, {.Error = 1e-3}).Vertices.empty());
  }
} // namespace Geoxy
