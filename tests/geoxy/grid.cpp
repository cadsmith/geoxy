/***************************************************************************************************
 *  @file      tests/geoxy/grid.cpp
 *  @brief     Tests that a grid's curves met all at once, or with one curve more, meet as every
 *             pair of them does, to the bit: on random grids of every kind at every scale, on the
 *             suites' curves, with curves of no box, and after curves are moved; met on pools of
 *             any size as alone, in the same order; and the curves a box picks.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/grid.hpp"

#include "../suites/cases.hpp"
#include "draws.hpp"
#include "gridwork.hpp"
#include "plain.hpp"
#include "sane.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/box.hpp"
#include "geoxy/contain.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/extents.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/loop.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/prepared.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/region.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/workers.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <thread>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  A pair that meets, as it was told.
     */
    struct Found
    {
      Grid::Pair     Curves  = {}; /**< The two curves' indices. */
      Geoxy::Marking Marking = {}; /**< What they meet in.       */
    };

    /**
     *  Whether two findings are the same, to the bit: the same curves, marks, stretches and budget.
     *
     *  @param  lhs One.
     *  @param  rhs The other.
     *  @return     Whether they are.
     */
    auto Same(const Found &lhs, const Found &rhs) -> bool
    {
      const auto overlap = [](const Overlap &one, const Overlap &other) -> bool
      {
        return one.First == other.First && one.Second == other.Second && one.Along == other.Along;
      };
      return lhs.Curves == rhs.Curves
          && lhs.Marking.Marks == rhs.Marking.Marks
          && std::ranges::equal(lhs.Marking.Overlaps, rhs.Marking.Overlaps, overlap)
          && lhs.Marking.Spent == rhs.Marking.Spent;
    }

    /**
     *  Whether two runs of findings are the same, one by one, in order.
     *
     *  @param  lhs One run.
     *  @param  rhs The other.
     *  @return     Whether they are.
     */
    auto Same(const std::vector<Found> &lhs, const std::vector<Found> &rhs) -> bool
    {
      return std::ranges::equal(lhs, rhs, [](const Found &one, const Found &other) -> bool
                                { return Same(one, other); });
    }

    /**
     *  A sink that keeps what it is told.
     *
     *  @param  into Where.
     *  @return      The sink.
     */
    auto Keeping(std::vector<Found> &into) -> Grid::Sink
    {
      return [&into](const Grid::Pair &curves, const Marking &marking) -> void
      {
        into.push_back({.Curves = curves, .Marking = marking});
      };
    }

    /**
     *  Every pair of a grid that meets, as IntersectAll tells them, in its order.
     *
     *  @param  grid      The grid.
     *  @param  tolerance The tolerance.
     *  @return           The pairs.
     */
    auto All(const Grid &grid, const Tolerance &tolerance) -> std::vector<Found>
    {
      std::vector<Found> found;
      IntersectAll(grid, tolerance, Keeping(found));
      return found;
    }

    /**
     *  Every pair of a grid that meets, met on the calling thread alone, in IntersectAll's order.
     *
     *  @param  grid      The grid.
     *  @param  tolerance The tolerance.
     *  @return           The pairs.
     */
    auto Alone(const Grid &grid, const Tolerance &tolerance) -> std::vector<Found>
    {
      std::vector<Found> found;
      Gridwork::Alone(grid, tolerance, Keeping(found));
      return found;
    }

    /**
     *  Every pair of a grid that meets, met on a pool, as IntersectAll tells them.
     *
     *  @param  by        The grid and the pool.
     *  @param  tolerance The tolerance.
     *  @return           The pairs.
     */
    auto Manned(const Grid::By &by, const Tolerance &tolerance) -> std::vector<Found>
    {
      std::vector<Found> found;
      IntersectAll(by, tolerance, Keeping(found));
      return found;
    }

    /**
     *  Every pair of a grid that meets, found by meeting every pair, the lower index first.
     *
     *  @param  grid      The grid.
     *  @param  tolerance The tolerance.
     *  @return           The pairs, in order of their indices.
     */
    auto Brute(const Grid &grid, const Tolerance &tolerance) -> std::vector<Found>
    {
      std::vector<Found> found;
      Marking            marking;
      for (std::size_t first = 0; first < grid.Size(); ++first)
      {
        for (std::size_t second = first + 1; second < grid.Size(); ++second)
        {
          Intersect({grid[first], grid[second]}, tolerance, marking);
          if (!marking.Marks.empty() || !marking.Overlaps.empty() || marking.Spent)
          {
            found.push_back({.Curves = {first, second}, .Marking = marking});
          }
        }
      }
      return found;
    }

    /**
     *  Findings put in order of their indices.
     *
     *  @param  found The findings.
     *  @return       Them, in order.
     */
    auto Sorted(std::vector<Found> found) -> std::vector<Found>
    {
      std::ranges::sort(found, {}, &Found::Curves);
      return found;
    }

    /**
     *  Whether IntersectAll tells a grid's pairs as meeting every pair does: the same pairs and
     *  marks, to the bit, and each kind's pairs in order of their indices.
     *
     *  @param  grid      The grid.
     *  @param  tolerance The tolerance.
     *  @return           Whether it does.
     */
    auto Agreed(const Grid &grid, const Tolerance &tolerance) -> bool
    {
      const std::vector<Found> all   = All(grid, tolerance);
      const std::vector<Found> brute = Brute(grid, tolerance);
      return Same(Sorted(all), brute);
    }

    /**
     *  Every curve of a grid a curve meets, as IntersectWith tells them.
     *
     *  @param  with      The curve and the grid.
     *  @param  tolerance The tolerance.
     *  @return           What it meets.
     */
    auto With(const Grid::With &with, const Tolerance &tolerance) -> std::vector<Found>
    {
      std::vector<Found> found;
      IntersectWith(with, tolerance, Keeping(found));
      return found;
    }

    /**
     *  Whether IntersectWith tells what a curve meets in a grid as meeting it with every curve
     *  does.
     *
     *  @param  with      The curve and the grid.
     *  @param  tolerance The tolerance.
     *  @return           Whether it does.
     */
    auto Probed(const Grid::With &with, const Tolerance &tolerance) -> bool
    {
      const Prepared    &curve = with.first;
      const Grid        &grid  = with.second;
      std::vector<Found> brute;
      Marking            marking;
      for (std::size_t index = 0; index < grid.Size(); ++index)
      {
        Intersect({curve, grid[index]}, tolerance, marking);
        if (!marking.Marks.empty() || !marking.Overlaps.empty() || marking.Spent)
        {
          brute.push_back({.Curves = {Grid::Outside, index}, .Marking = marking});
        }
      }
      return Same(With(with, tolerance), brute);
    }

    /**
     *  Random curves of every kind, about a place.
     */
    class Scatter
    {
    public:
      /**
       *  Curves drawn from a seed.
       *
       *  @param seed The seed.
       */
      explicit Scatter(std::uint64_t seed) : m_Draws(seed) {}

      /**
       *  A whole number below a count.
       *
       *  @param  count The count.
       *  @return       The number.
       */
      auto Pick(std::uint64_t count) -> std::uint64_t { return m_Draws.Next() % count; }

      /**
       *  A double about an offset.
       *
       *  @param  offset The offset.
       *  @param  spread How far either way.
       *  @return        The double.
       */
      auto Near(double offset, double spread) -> double { return m_Draws.Near(offset, spread); }

      /**
       *  A point within a box.
       *
       *  @param  box The box.
       *  @return     The point.
       */
      auto Within(const Box &box) -> Point
      {
        const Point  middle = Center(box);
        const Vector half   = (box.Max - box.Min) / 2;
        return {.X = Near(middle.X, half.X), .Y = Near(middle.Y, half.Y)};
      }

      /**
       *  A curve of any bounded kind, of about a size, about a point: lines, a fifth of them
       *  horizontal or upright, arcs and circles, ellipses, Bézier curves, splines and polylines.
       *
       *  @param  where About where.
       *  @param  size  About how large.
       *  @return       The curve.
       */
      auto Bounded(Point where, double size) -> Curve
      {
        const auto off = [this, where, size]() -> Point
        {
          return {.X = Near(where.X, size), .Y = Near(where.Y, size)};
        };
        switch (Pick(10))
        {
        case 0:
        {
          const Point start = off();
          const Point end   = off();
          return Pick(2) == 0 ? Line{.Start = start, .End = {.X = end.X, .Y = start.Y}}
                              : Line{.Start = start, .End = {.X = start.X, .Y = end.Y}};
        }
        case 1:
        case 2:
          return Line{.Start = off(), .End = off()};
        case 3:
          return Arc{.Center = where,
                     .Radius = size * Near(1, 0.9),
                     .Start  = Degrees(Near(0, 180)),
                     .Sweep  = Pick(3) == 0 ? Turn : Degrees(Near(0, 300))};
        case 4:
          return Ellipse{.Map   = Translation(where - Point{})
                                * Rotation(Degrees(Near(0, 90)))
                                * Scaling(size * Near(1, 0.8), size * Near(0.6, 0.5)),
                         .Start = Degrees(Near(0, 180)),
                         .Sweep = Pick(2) == 0 ? Turn : Degrees(Near(0, 300))};
        case 5:
          return Bezier::Quadratic{.Start = off(), .Control = off(), .End = off()};
        case 6:
          return Bezier::Cubic{.Start = off(), .Out = off(), .In = off(), .End = off()};
        case 7:
          return Spline{.Degree  = 2,
                        .Knots   = {0, 0, 0, 0.4, 1, 1, 1},
                        .Poles   = {off(), off(), off(), off()},
                        .Weights = Pick(2) == 0 ? std::vector<double>{}
                                                : std::vector<double>{1, Near(1, 0.5), 2, 1}};
        default:
        {
          Polyline polyline = {.Vertices = {}, .Closed = Pick(3) == 0};
          for (std::uint64_t index = 0; index < 2 + Pick(3); ++index)
          {
            polyline.Vertices.push_back(
                {.Position = off(), .Bulge = Pick(2) == 0 ? 0.0 : Near(0, 1)});
          }
          return polyline;
        }
        }
      }

      /**
       *  Lines as a grid's lines come: slanted at random, horizontal and upright, starting on
       *  another line or at its end, near parallel to another, and lines crossing a level line
       *  exactly halfway between two doubles.
       *
       *  @param  count  About how many.
       *  @param  region Where they lie.
       *  @param  size   About how long each is.
       *  @return        The lines.
       */
      auto Strokes(std::size_t count, const Box &region, double size) -> std::vector<Curve>
      {
        std::vector<Curve> lines;
        const auto         heading = [this, size](Point start) -> Line
        {
          const double angle = Near(0, std::numbers::pi);
          const double reach = size * Near(1, 0.5);
          return {.Start = start,
                  .End   = {.X = start.X + (reach * std::cos(angle)),
                            .Y = start.Y + (reach * std::sin(angle))}};
        };
        while (lines.size() < count)
        {
          const Point  where = Within(region);
          const double reach = size * Near(1, 0.5);
          const Line  *other =
              lines.empty() ? nullptr : &std::get<Line>(lines.at(Pick(lines.size())));
          switch (other == nullptr ? 0 : Pick(10))
          {
          case 0:
          case 1:
          case 2:
            lines.emplace_back(heading(where));
            break;
          case 3:
            lines.emplace_back(Line{.Start = where, .End = {.X = where.X + reach, .Y = where.Y}});
            break;
          case 4:
            lines.emplace_back(Line{.Start = where, .End = {.X = where.X, .Y = where.Y + reach}});
            break;
          case 5:
            lines.emplace_back(
                heading(other->Start + ((other->End - other->Start) * Near(0.5, 0.5))));
            break;
          case 6:
            lines.emplace_back(heading(Pick(2) == 0 ? other->Start : other->End));
            break;
          case 7:
          {
            const Vector way   = other->End - other->Start;
            const double angle = std::ldexp(Near(0, 1), -static_cast<int>(Pick(30)));
            const Vector turn  = {.X = way.X - (way.Y * angle), .Y = way.Y + (way.X * angle)};
            const Point  mid   = Midpoint(other->Start, other->End);
            lines.emplace_back(Line{.Start = mid - (turn / 2), .End = mid + (turn / 2)});
            break;
          }
          default:
          {
            const double next = std::nextafter(where.X, std::numeric_limits<double>::infinity());
            lines.emplace_back(Line{.Start = {.X = where.X, .Y = where.Y - reach},
                                    .End   = {.X = next, .Y = where.Y + reach}});
            lines.emplace_back(Line{.Start = {.X = where.X - reach, .Y = where.Y},
                                    .End   = {.X = where.X + reach, .Y = where.Y}});
          }
          }
        }
        return lines;
      }

      /**
       *  A curve with no box: a ray or an endless line.
       *
       *  @param  where Where it starts, or runs through.
       *  @return       The curve.
       */
      auto Endless(Point where) -> Curve
      {
        const Vector way = {.X = Near(0, 1), .Y = Near(0, 1)};
        if (Pick(2) == 0)
        {
          return Ray{.Start = where, .Way = way};
        }
        return Xline{.Through = where, .Way = way};
      }

      /**
       *  A grid's curves: some in clusters about a few places, some copies of others, some starting
       *  where others end, and a few with no box.
       *
       *  @param  count  How many.
       *  @param  region Where they lie.
       *  @param  size   About how large each is.
       *  @return        The curves.
       */
      auto Curves(std::size_t count, const Box &region, double size) -> std::vector<Curve>
      {
        std::vector<Curve> curves;
        std::vector<Point> places(1 + Pick(4));
        for (Point &place : places)
        {
          place = Within(region);
        }
        while (curves.size() < count)
        {
          const Point where = Pick(3) == 0 ? Within(region) : places.at(Pick(places.size()));
          switch (Pick(12))
          {
          case 0:
            if (!curves.empty())
            {
              curves.push_back(curves.at(Pick(curves.size())));
              break;
            }
            [[fallthrough]];
          case 1:
            if (const auto *line = curves.empty() ? nullptr : std::get_if<Line>(&curves.back()))
            {
              curves.emplace_back(Line{.Start = line->End, .End = Within(region)});
              break;
            }
            [[fallthrough]];
          default:
            curves.push_back(Bounded(where, size));
          }
        }
        return curves;
      }

    private:
      Draws m_Draws; /**< The draws. */
    };

    /**
     *  Rooms of a floor plan, as the benchmarks draw them: double walls, a door's swing and leaf, a
     *  column, a line and an arc of furniture; and a few lines across the whole plan.
     *
     *  @param  rooms How many rooms.
     *  @param  seed  The seed.
     *  @return       The curves, ten a room.
     */
    auto Plan(std::size_t rooms, std::uint64_t seed) -> std::vector<Curve>
    {
      Draws              draws{seed};
      const auto         side = static_cast<std::size_t>(std::ceil(std::sqrt(rooms)));
      std::vector<Curve> curves;
      const auto         unit = [&draws]() -> double
      {
        return draws.Near(0.5, 0.5);
      };
      for (std::size_t room = 0; room < rooms; ++room)
      {
        const std::size_t row = room / side;
        const double      x   = (10.0 * static_cast<double>(room - (row * side))) + (0.05 * unit());
        const double      y   = (10.0 * static_cast<double>(row)) + (0.05 * unit());
        curves.emplace_back(Line{.Start = {.X = x, .Y = y}, .End = {.X = x + 10, .Y = y}});
        curves.emplace_back(
            Line{.Start = {.X = x, .Y = y + 0.2}, .End = {.X = x + 3, .Y = y + 0.2}});
        curves.emplace_back(
            Line{.Start = {.X = x + 4, .Y = y + 0.2}, .End = {.X = x + 10, .Y = y + 0.2}});
        curves.emplace_back(Line{.Start = {.X = x, .Y = y}, .End = {.X = x, .Y = y + 10}});
        curves.emplace_back(
            Line{.Start = {.X = x + 0.2, .Y = y + 0.2}, .End = {.X = x + 0.2, .Y = y + 10}});
        curves.emplace_back(Arc{.Center = {.X = x + 3, .Y = y + 0.2},
                                .Radius = 1.0,
                                .Start  = Degrees(0),
                                .Sweep  = Degrees(90)});
        curves.emplace_back(
            Line{.Start = {.X = x + 3, .Y = y + 0.2}, .End = {.X = x + 3, .Y = y + 1.2}});
        curves.emplace_back(
            Circle({.X = x + 0.1 + (0.4 * unit()), .Y = y + 0.1 + (0.4 * unit())}, 0.35));
        const Point desk = {.X = x + 2 + (6 * unit()), .Y = y + 2 + (6 * unit())};
        curves.emplace_back(
            Line{.Start = desk,
                 .End   = {.X = desk.X + (3 * unit()) - 1.5, .Y = desk.Y + (3 * unit()) - 1.5}});
        curves.emplace_back(Arc{.Center = {.X = x + 2 + (6 * unit()), .Y = y + 2 + (6 * unit())},
                                .Radius = 0.5 + (2 * unit()),
                                .Start  = Degrees(360 * unit()),
                                .Sweep  = Degrees(30 + (300 * unit()))});
      }
      const double size = 10.0 * static_cast<double>(side);
      for (std::size_t across = 0; across < 1 + (rooms / 200); ++across)
      {
        const double where = size * unit();
        curves.emplace_back(
            Line{.Start = {.X = -5, .Y = where}, .End = {.X = size + 5, .Y = where}});
        curves.emplace_back(
            Line{.Start = {.X = where, .Y = -5}, .End = {.X = where, .Y = size + 5}});
      }
      return curves;
    }

    /**
     *  A box about the origin.
     *
     *  @param  half Half its side.
     *  @return      The box.
     */
    auto About(double half) -> Box
    {
      return {.Min = {.X = -half, .Y = -half}, .Max = {.X = half, .Y = half}};
    }

    /**
     *  The extents of some curves, those with a box.
     *
     *  @param  curves The curves.
     *  @return        Their extents.
     */
    auto Spanned(const std::vector<Curve> &curves) -> Box
    {
      Box all;
      for (const Curve &curve : curves)
      {
        const Box box = Extents(curve);
        if (std::isfinite(box.Min.X + box.Min.Y + box.Max.X + box.Max.Y))
        {
          all = Joined(all, box);
        }
      }
      return all;
    }
  } // namespace

  TEST_CASE("grid: an empty grid and a grid of one curve meet in nothing")
  {
    const Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
    const Grid      empty{{}, tolerance};
    CHECK(empty.Size() == 0);
    CHECK(All(empty, tolerance).empty());
    const Curve    cross = Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 1, .Y = 1}};
    const Prepared ready{cross};
    CHECK(With({ready, empty}, tolerance).empty());

    const std::vector<Curve> one = {Line{.Start = {.X = 0, .Y = 1}, .End = {.X = 1, .Y = 0}}};
    const Grid               alone{one, tolerance};
    CHECK(alone.Size() == 1);
    CHECK(All(alone, tolerance).empty());
    const std::vector<Found> met = With({ready, alone}, tolerance);
    REQUIRE(met.size()       == 1);
    CHECK(met.front().Curves == Grid::Pair{Grid::Outside, 0});
    CHECK(met.front().Marking.Marks == std::vector{Mark{.At = {.X = 0.5, .Y = 0.5}}});
    CHECK(&alone[0].Shape() != &one.front());
  }

  TEST_CASE("grid: random grids of every kind meet as every pair does, at every scale")
  {
    Scatter     scatter{0xD8A3};
    int         differ = 0;
    std::size_t met    = 0;
    std::size_t pairs  = 0;
    const int   rounds = Rounds(32);
    for (int round = 0; round < rounds; ++round)
    {
      const double             scale  = std::ldexp(1.0, static_cast<int>(scatter.Pick(80)) - 40);
      const double             offset = scatter.Pick(4) == 0 ? scale * 1e6 : 0.0;
      const Box                region = About(scale);
      const Box                placed = {.Min = region.Min + Vector{.X = offset, .Y = 0},
                                         .Max = region.Max + Vector{.X = offset, .Y = 0}};
      const std::size_t        count  = 20 + scatter.Pick(60);
      const std::vector<Curve> curves = scatter.Curves(count, placed, scale / 4);
      const Tolerance          sized  = Tolerance::Of(Spanned(curves));
      const std::array<Tolerance, 3> tolerances = {
          sized,
          Tolerance{.Length = scale * 1e-4 , .Radians = 1e-9},
          Tolerance{.Length = scale * 1e-13, .Radians = 1e-9},
      };
      const Tolerance &tolerance = tolerances.at(scatter.Pick(tolerances.size()));
      CAPTURE(round);
      const Grid grid{curves, tolerance};
      differ += Agreed(grid, tolerance) ? 0 : 1;
      met    += All(grid, tolerance).size();
      pairs  += count * (count - 1) / 2;
      /**
       *  Met within a finer tolerance than the tiling was made for, and a coarser one.
       */
      for (const double times : {0.25, 64.0})
      {
        const Tolerance other  = {.Length = tolerance.Length * times, .Radians = 1e-9};
        differ                += Agreed(grid, other) ? 0 : 1;
      }
    }
    CHECK(differ == 0);
    CHECK(met > pairs / 50);
  }

  /**
   *  Intersect tells two lines or arcs far apart past four lengths between their boxes, any other
   *  pair past sixteen, and the tiling grows each box by thirty-two: pairs set about each of those
   *  gaps, along x and along y, at sizes where the length is the tolerance's and where it is the
   *  grain's, are met as every pair is.
   */
  TEST_CASE("grid: pairs set about each gap that decides are met as every pair is")
  {
    int differ = 0;
    for (const int power : {-60, 0, 30, 60})
    {
      const double scale = std::ldexp(1.0, power);
      for (const double length : {scale * 1e-6, 0.0})
      {
        const Tolerance    tolerance = {.Length = length, .Radians = 1e-9};
        const double       grain     = 16 * std::ldexp(1.0, power - 52);
        const double       unit      = std::max(length, grain);
        std::vector<Curve> curves = {Line{.Start = {.X = 0, .Y = 0}, .End = {.X = scale, .Y = 0}},
                                     Bezier::Quadratic{.Start   = {.X = 0, .Y = 0},
                                                       .Control = {.X = scale / 2, .Y = scale / 4},
                                                       .End     = {.X = scale, .Y = 0}}};
        for (const double gap : {3.99, 4.0, 4.01, 15.99, 16.0, 16.01, 63.99, 64.0, 64.01})
        {
          const double off = scale + (gap * unit);
          curves.emplace_back(
              Line{.Start = {.X = off, .Y = -scale}, .End = {.X = off, .Y = scale}});
          curves.emplace_back(Line{.Start = {.X = 0, .Y = -off}, .End = {.X = scale, .Y = -off}});
        }
        const Grid grid{curves, tolerance};
        differ += Agreed(grid, tolerance) ? 0 : 1;
      }
    }
    CHECK(differ == 0);
  }

  /**
   *  Lines are met kind by kind through their exact pass, which tells most pairs and leaves the
   *  rest to Intersect: every pair it tells must be told as Intersect tells it, and IntersectAll as
   *  every pair.
   */
  TEST_CASE("grid: lines met by their exact pass meet as every pair does, at every scale")
  {
    Scatter     scatter{0x11E5};
    int         differ = 0;
    int         wrong  = 0;
    std::size_t told   = 0;
    std::size_t lines  = 0;
    const int   rounds = Rounds(48);
    for (int round = 0; round < rounds; ++round)
    {
      const double             scale  = std::ldexp(1.0, static_cast<int>(scatter.Pick(80)) - 40);
      const double             offset = scatter.Pick(3) == 0 ? scale * 1e6 : 0.0;
      const Box                region = {.Min = {.X = offset - scale, .Y = -scale},
                                         .Max = {.X = offset + scale, .Y = scale}};
      const std::vector<Curve> curves = scatter.Strokes(
          40 + scatter.Pick(80), region, scale / static_cast<double>(2 + scatter.Pick(6)));
      const std::array<Tolerance, 4> tolerances = {
          Tolerance::Of(Spanned(curves)),
          Tolerance{.Length = scale * 1e-4 , .Radians = 1e-9},
          Tolerance{.Length = scale * 1e-13, .Radians = 1e-9},
          Tolerance{.Length = scale * 1e-6 , .Radians =  0.1},
      };
      const Tolerance &tolerance = tolerances.at(scatter.Pick(tolerances.size()));
      CAPTURE(round);
      const Grid grid{curves, tolerance};
      differ += Agreed(grid, tolerance) ? 0 : 1;
      Marking marking;
      Marking plainly;
      for (std::size_t first = 0; first < grid.Size(); ++first)
      {
        for (std::size_t second = first + 1; second < grid.Size(); ++second)
        {
          const Prepared::Pair pair  = {grid[first], grid[second]};
          const Plain          plain = Plainly(pair, Sane(tolerance), plainly);
          Intersect(pair, tolerance, marking);
          if (!plain.Told)
          {
            wrong += Same({.Marking = marking}, {.Marking = plainly}) ? 0 : 1;
            continue;
          }
          const std::vector<Mark> marks = plain.Count == 0
                                            ? std::vector<Mark>{}
                                            : std::vector{Mark{.At = plain.At, .Kind = plain.Kind}};
          wrong += marking.Marks == marks && marking.Overlaps.empty() && !marking.Spent ? 0 : 1;
          told  += plain.Count;
        }
      }
      lines += curves.size();
    }
    CHECK(differ == 0);
    CHECK(wrong  == 0);
    CHECK(told > lines / 4);
  }

  /**
   *  What the pass leaves to Intersect: crossings exactly halfway between two doubles, a slanted
   *  line near parallel crossing another near its end, a tolerance whose angle is past what the
   *  pass allows for, and numbers too large or too small for it; and what the slanted pass tells of
   *  ends on a line along an axis, as Intersect does.
   */
  TEST_CASE("grid: lines the exact pass leaves are met as every pair is")
  {
    const Tolerance          tolerance = {.Length = 1e-9, .Radians = 1e-9};
    const Tolerance          sane      = Sane(tolerance);
    const double             next      = std::nextafter(1.0, 2.0);
    const std::vector<Curve> curves    = {
        Line{.Start = {.X = 1, .Y = -1}, .End = {.X = next, .Y = 1}},
        Line{.Start = {.X = -1, .Y = 0}, .End = {.X = 3, .Y = 0}},
        Line{.Start = {.X = -1, .Y = -1}, .End = {.X = 1, .Y = 1}},
        Line{.Start = {.X = -1, .Y = 1}, .End = {.X = 1, .Y = -1}},
        Line{.Start = {.X = 2, .Y = -1}, .End = {.X = 2, .Y = 0}},
        Line{.Start = {.X = 2.5, .Y = 0}, .End = {.X = 2.5, .Y = 2}},
        Line{.Start = {.X = -0.5, .Y = -0.01}, .End = {.X = 2.5, .Y = 0.02}},
        Line{.Start = {.X = 0, .Y = -1}, .End = {.X = 0, .Y = 2}},
        Line{.Start = {.X = 1.5, .Y = 1}, .End = {.X = 3.5, .Y = 1.25}},
    };
    const Grid grid{curves, tolerance};
    Marking    marking;
    const auto plainly = [&grid, &sane, &marking](std::size_t first, std::size_t second) -> Plain
    {
      return Plainly({grid[first], grid[second]}, sane, marking);
    };
    CHECK(!plainly(0, 1).Told);
    CHECK(plainly(2, 3).Count == 1);
    CHECK(plainly(1, 4).Count == 1);
    CHECK(plainly(1, 5).Count == 1);
    CHECK(plainly(1, 6).Count == 1);
    CHECK(plainly(5, 8).Count == 1);
    CHECK(plainly(1, 7).Count == 1);
    CHECK(plainly(1, 7).At    == Point{.X = 0, .Y = 0});
    CHECK(!Plainly({grid[5], grid[8]}, {.Length = 1e-9, .Radians = 0.1}, marking).Told);
    CHECK(marking.Marks.size() == 1);
    CHECK(Agreed(grid, tolerance));
    CHECK(Agreed(grid, {.Length = 1e-9, .Radians = 0.1}));

    /**
     *  Lines of numbers so large, or so small, that the length they are met within is past what the
     *  pass takes; and lines of numbers that are not.
     */
    for (const int power : {-1060, -1000, 400, 1000})
    {
      CAPTURE(power);
      std::vector<Curve> scaled;
      for (const Curve &curve : curves)
      {
        const Line &line = std::get<Line>(curve);
        const auto  far  = [power](Point point) -> Point
        {
          return {.X = std::ldexp(point.X, power), .Y = std::ldexp(point.Y, power)};
        };
        scaled.emplace_back(Line{.Start = far(line.Start), .End = far(line.End)});
      }
      for (const double length : {0.0, std::ldexp(1e-9, power)})
      {
        const Tolerance within = {.Length = length, .Radians = 1e-9};
        const Grid      apart{scaled, within};
        CHECK(!Plainly({apart[5], apart[8]}, Sane(within), marking).Told);
        CHECK(Agreed(apart, within));
      }
    }
    const double             none  = std::numeric_limits<double>::quiet_NaN();
    const std::vector<Curve> loose = {
        Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 2, .Y = 1}},
        Line{.Start = {.X = none, .Y = 0}, .End = {.X = 1, .Y = 1}},
        Line{.Start = {.X = 0, .Y = 1}, .End = {.X = 2, .Y = 0}},
    };
    const Grid lost{loose, tolerance};
    CHECK(!Plainly({lost[0], lost[1]}, sane, marking).Told);
    CHECK(!Plainly({lost[1], lost[2]}, sane, marking).Told);
    CHECK(Agreed(lost, tolerance));
  }

  TEST_CASE("grid: curves with no box are met with every other")
  {
    Scatter                  scatter{0xB0C5};
    const Tolerance          tolerance = {.Length = 1e-9, .Radians = 1e-9};
    std::vector<Curve>       curves    = scatter.Curves(30, About(10), 3);
    const double             none      = std::numeric_limits<double>::quiet_NaN();
    const std::vector<Curve> loose     = {
        scatter.Endless({.X = 0, .Y = 0}),
        scatter.Endless({.X = 5, .Y = -3}),
        Ray{.Start = {.X = -20, .Y = 1}, .Way = {.X = 1, .Y = 0}},
        Xline{.Through = {.X = 1, .Y = 1}, .Way = {.X = 0, .Y = 1}},
        Line{.Start = {.X = none, .Y = 0}, .End = {.X = 1, .Y = 1}},
        Arc{.Center = {.X = 1.7e308, .Y = -1.7e308}, .Radius = 1e308, .Start = {}, .Sweep = Turn},
        Polyline{},
    };
    for (std::size_t index = 0; index < loose.size(); ++index)
    {
      curves.insert(curves.begin() + static_cast<std::ptrdiff_t>(3 * index), loose.at(index));
    }
    const Grid grid{curves, tolerance};
    CHECK(Agreed(grid, tolerance));
    CHECK(Agreed(grid, {.Length = 1e-3, .Radians = 1e-9}));

    /**
     *  A grid of curves with no box only, and one met within a tolerance that is not a number,
     *  where every curve is.
     */
    const Grid endless{loose, tolerance};
    CHECK(Agreed(endless, tolerance));
    const Grid lost{curves, {.Length = none, .Radians = 1e-9}};
    CHECK(Agreed(lost, tolerance));
    CHECK(Agreed(grid, {.Length = none, .Radians = 1e-9}));
  }

  TEST_CASE("grid: clusters far apart, and curves of huge and tiny numbers")
  {
    Scatter scatter{0xFA2};
    for (const double apart : {1e3, 1e9, 1e15, 1e300})
    {
      CAPTURE(apart);
      std::vector<Curve> curves = scatter.Curves(20, About(1), 0.4);
      const Box          away   = {.Min = {.X = apart - 1, .Y = apart - 1},
                                   .Max = {.X = apart + 1, .Y = apart + 1}};
      for (const Curve &curve : scatter.Curves(20, away, 0.4))
      {
        curves.push_back(curve);
      }
      const Tolerance       tolerance = Tolerance::Of(Spanned(curves));
      const Grid            grid{curves, tolerance};
      const Gridwork::Pairs pairs = Gridwork::Candidates(grid, tolerance);
      std::size_t           count = 0;
      for (const std::vector<Grid::Pair> &kind : pairs)
      {
        count += kind.size();
      }
      CHECK(Agreed(grid, tolerance));
      CHECK(Agreed(grid, {.Length = 1e-9, .Radians = 1e-9}));
      CHECK(count < curves.size() * curves.size() / 4);
    }
    for (const int power : {-1000, -300, -60, 60, 300, 1000})
    {
      CAPTURE(power);
      const double             scale  = std::ldexp(1.0, power);
      const std::vector<Curve> curves = scatter.Curves(24, About(scale), scale / 4);
      for (const Tolerance &tolerance :
           {Tolerance::Of(Spanned(curves)), Tolerance{.Length = scale * 1e-6, .Radians = 1e-9},
            Tolerance{.Length = 0, .Radians = 0}})
      {
        const Grid grid{curves, tolerance};
        CHECK(Agreed(grid, tolerance));
      }
    }
    /**
     *  Points scattered thin over a wide plane, which would ask for more cells than curves.
     */
    std::vector<Curve> specks;
    for (int index = 0; index < 60; ++index)
    {
      const Point where = scatter.Within(About(1e6));
      specks.emplace_back(Line{.Start = where, .End = where + Vector{.X = 1e-3, .Y = 0}});
    }
    specks.emplace_back(Line{.Start = {.X = -1e6, .Y = 0}, .End = {.X = 1e6, .Y = 1}});
    const Tolerance fine = {.Length = 1e-9, .Radians = 1e-9};
    CHECK(Agreed(Grid{specks, fine}, fine));
  }

  TEST_CASE("grid: the suites' curves, met as grids, meet as every pair of them does")
  {
    const std::filesystem::path suites               = std::filesystem::path{GEOXY_FUZZ_DIR} / "..";
    const std::array<std::filesystem::path, 1> roots = {suites / "suites"};
    int                                        differ = 0;
    int                                        met    = 0;
    for (const std::filesystem::path &file : Suites::Files(roots))
    {
      CAPTURE(file);
      std::vector<Curve> curves;
      for (const Suites::Case &each : Suites::Read(file))
      {
        for (const std::optional<Curve> &curve : {each.First, each.Second})
        {
          if (curve)
          {
            curves.push_back(*curve);
          }
        }
      }
      /**
       *  The curves in runs of 96, about as many cases as lie in one place, each run a grid.
       */
      constexpr std::size_t Run   = 96;
      const std::size_t     runs  = (curves.size() + Run - 1) / Run;
      const auto            taken = static_cast<std::size_t>(Rounds(static_cast<int>(runs)));
      for (std::size_t run = 0; run < runs; run += runs / taken)
      {
        const std::size_t        first = run * Run;
        const std::vector<Curve> some(
            curves.begin() + static_cast<std::ptrdiff_t>(first),
            curves.begin() + static_cast<std::ptrdiff_t>(std::min(curves.size(), first + Run)));
        const Tolerance          tolerance = Tolerance::Of(Spanned(some));
        const Grid               grid{some, tolerance};
        const std::vector<Found> all  = All(grid, tolerance);
        differ                       += Same(Sorted(all), Brute(grid, tolerance)) ? 0 : 1;
        met                          += static_cast<int>(all.size());
      }
    }
    CHECK(differ == 0);
    CHECK(met > 0);
  }

  /**
   *  Rooms of a floor plan, their walls meeting at corners and joints, each door's swing ending on
   *  its wall and its leaf, at several scales and far from the origin: IntersectAll tells every
   *  pair as Intersect tells the two curves themselves, to the bit, the corners and joints as hits
   *  at an end, each told by the pass for ends with no marking made.
   */
  TEST_CASE("grid: a plan's corners, joints and door swings meet as Intersect meets their curves")
  {
    Draws draws{0xD0025};
    for (const int power : {-30, 0, 30})
    {
      for (const double offset : {0.0, 1e4})
      {
        CAPTURE(power);
        CAPTURE(offset);
        const double       scale = std::ldexp(1.0, power);
        std::vector<Curve> curves;
        for (std::size_t room = 0; room < 12; ++room)
        {
          const std::size_t row    = room / 4;
          const auto        column = static_cast<double>(room - (4 * row));
          const double      x      = offset + (10.0 * column) + draws.Near(0.5, 0.5);
          const double      y = offset + (10.0 * static_cast<double>(row)) + draws.Near(0.5, 0.5);
          const auto        spot = [x, y, scale](double right, double above) -> Point
          {
            return {.X = (x + right) * scale, .Y = (y + above) * scale};
          };
          curves.emplace_back(Line{.Start = spot(0, 0), .End = spot(10, 0)});
          curves.emplace_back(Line{.Start = spot(0, 0.25), .End = spot(3, 0.25)});
          curves.emplace_back(Line{.Start = spot(4, 0.25), .End = spot(10, 0.25)});
          curves.emplace_back(Line{.Start = spot(0, 0), .End = spot(0, 10)});
          curves.emplace_back(Line{.Start = spot(0.25, 0.25), .End = spot(0.25, 10)});
          curves.emplace_back(Arc{
              .Center = spot(3, 0.25), .Radius = scale, .Start = Degrees(0), .Sweep = Degrees(90)});
          curves.emplace_back(Line{.Start = spot(3, 0.25), .End = spot(3, 1.25)});
        }
        const Tolerance    tolerance = Tolerance::Of(Spanned(curves));
        const Grid         grid{curves, tolerance};
        std::vector<Found> expected;
        int                ends = 0;
        for (std::size_t first = 0; first < curves.size(); ++first)
        {
          for (std::size_t second = first + 1; second < curves.size(); ++second)
          {
            const Intersections meeting = Intersect(curves.at(first), curves.at(second), tolerance);
            Marking             marking;
            for (const Hit &hit : meeting.Hits)
            {
              marking.Marks.push_back({.At = hit.At, .Kind = hit.Kind});
              ends += hit.Kind == Contact::End ? 1 : 0;
            }
            marking.Overlaps.assign(meeting.Overlaps.begin(), meeting.Overlaps.end());
            marking.Spent = meeting.Spent;
            if (!marking.Marks.empty() || !marking.Overlaps.empty() || marking.Spent)
            {
              expected.push_back({.Curves = {first, second}, .Marking = marking});
            }
          }
        }
        CHECK(Same(Sorted(All(grid, tolerance)), expected));
        CHECK(ends > 48);
        Marking     marking;
        const Plain corner = Plainly({grid[0], grid[3]}, Sane(tolerance), marking);
        CHECK(corner.Told);
        CHECK(corner.Kind == Contact::End);
      }
    }
  }

  TEST_CASE("grid: pairs are sorted by kind, each kind in order of the curves' indices")
  {
    const Tolerance          tolerance = {.Length = 1e-9, .Radians = 1e-9};
    const std::vector<Curve> curves    = {
        Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 10, .Y = 0}},
        Arc{.Center = {.X = 5, .Y = 0}, .Radius = 2, .Start = {}, .Sweep = Turn},
        Line{.Start = {.X = 3, .Y = -5}, .End = {.X = 3, .Y = 5}},
        Line{.Start = {.X = 0, .Y = -1}, .End = {.X = 10, .Y = 1}},
        Bezier::Cubic{.Start = {.X = 0, .Y = 3},
                      .Out   = {.X = 3, .Y = -3},
                      .In    = {.X = 6, .Y = 6},
                      .End   = {.X = 9, .Y = -2}},
        Arc{.Center = {.X = 6, .Y = 0}, .Radius = 2, .Start = {}, .Sweep = Turn},
        Line{.Start = {.X = 4, .Y = 4}, .End = {.X = 4, .Y = 4}},
    };
    const Grid            grid{curves, tolerance};
    const Gridwork::Pairs pairs = Gridwork::Candidates(grid, tolerance);
    const auto            kind  = [&pairs](Gridwork::Kind which) -> const std::vector<Grid::Pair> &
    {
      return pairs.at(static_cast<std::size_t>(which));
    };
    CHECK(kind(Gridwork::Kind::Square) == std::vector<Grid::Pair>{{0, 2}});
    CHECK(kind(Gridwork::Kind::Lines)  == std::vector<Grid::Pair>{{0, 3}, {2, 3}});
    CHECK(kind(Gridwork::Kind::Lined)
          == std::vector<Grid::Pair>{{0, 1}, {0, 5}, {1, 2}, {1, 3}, {3, 5}});
    CHECK(kind(Gridwork::Kind::Rounds)       == std::vector<Grid::Pair>{{1, 5}});
    CHECK(kind(Gridwork::Kind::Other).size() == 6);
    std::vector<Grid::Pair> told;
    IntersectAll(grid, tolerance, [&told](const Grid::Pair &pair, const Marking &) -> void
                 { told.push_back(pair); });
    CHECK(told.front() == Grid::Pair{0, 2});
    CHECK(Agreed(grid, tolerance));
  }

  TEST_CASE("grid: one curve met with a grid, within its tolerance, a finer and a coarser")
  {
    Scatter                  scatter{0x9A7E};
    const std::vector<Curve> curves  = scatter.Curves(60, About(10), 3);
    int                      differ  = 0;
    const int                rounds  = Rounds(40);
    const double             none    = std::numeric_limits<double>::quiet_NaN();
    const Tolerance          drawn   = {.Length = 1e-6, .Radians = 1e-9};
    std::vector<Curve>       holding = curves;
    holding.emplace_back(Ray{.Start = {.X = 0, .Y = 0}, .Way = {.X = 1, .Y = 2}});
    const Grid grid{holding, drawn};
    for (int round = 0; round < rounds; ++round)
    {
      const Curve    probe = scatter.Pick(8) == 0 ? scatter.Endless(scatter.Within(About(10)))
                                                  : scatter.Bounded(scatter.Within(About(10)), 4);
      const Prepared ready{probe};
      for (const double length : {1e-6, 1e-9, 1e-3, 0.5, none})
      {
        differ += Probed({ready, grid}, {.Length = length, .Radians = 1e-9}) ? 0 : 1;
      }
    }
    CHECK(differ == 0);
    const Curve   &copy = curves.front();
    const Prepared ready{copy};
    const auto     met = With({ready, grid}, drawn);
    REQUIRE(!met.empty());
    CHECK(std::ranges::any_of(met, [](const Found &found) -> bool
                              { return found.Curves == Grid::Pair{Grid::Outside, 0}; }));
    const Grid lost{curves, {.Length = none, .Radians = 1e-9}};
    CHECK(Probed({ready, lost}, drawn));
  }

  TEST_CASE("grid: curves moved meet as a grid made again of where they are")
  {
    Scatter            scatter{0x40E5};
    const Tolerance    tolerance = {.Length = 1e-7, .Radians = 1e-9};
    std::vector<Curve> curves    = scatter.Curves(40, About(20), 3);
    Grid               grid{curves, tolerance};
    int                differ = 0;
    const int          rounds = Rounds(160);
    for (int round = 0; round < rounds; ++round)
    {
      const std::size_t   index = scatter.Pick(curves.size());
      const std::uint64_t how   = scatter.Pick(10);
      Curve               curve;
      if (how == 0)
      {
        curve = scatter.Endless(scatter.Within(About(20)));
      }
      else if (how == 1)
      {
        curve = scatter.Bounded(scatter.Within(About(400)), 3);
      }
      else
      {
        curve = scatter.Bounded(scatter.Within(About(20)), scatter.Pick(4) == 0 ? 15 : 3);
      }
      curves.at(index) = curve;
      grid.Moved(index, curve);
      if (round % 16 == 0 || round == rounds - 1)
      {
        const Grid again{curves, tolerance};
        differ               += Same(All(grid, tolerance), All(again, tolerance)) ? 0 : 1;
        differ               += Same(All(grid, {.Length = 1e-4, .Radians = 1e-9}),
                                     All(again, {.Length = 1e-4, .Radians = 1e-9}))
                                  ? 0
                                  : 1;
        const Curve    probe  = scatter.Bounded(scatter.Within(About(20)), 6);
        const Prepared ready{probe};
        differ += Same(With({ready, grid}, tolerance), With({ready, again}, tolerance)) ? 0 : 1;
        differ += Agreed(grid, tolerance) ? 0 : 1;
      }
    }
    CHECK(differ == 0);

    /**
     *  One curve moved again and again into cells that are full, which move to the end of the
     *  entries, until the room they leave is given back.
     */
    for (int step = 0; step < 400; ++step)
    {
      const Point where = {.X = -20 + (0.1 * step), .Y = 20 - (0.1 * step)};
      const Curve curve = Line{.Start = where, .End = where + Vector{.X = 1, .Y = 0.5}};
      curves.at(7)      = curve;
      grid.Moved(7, curve);
    }
    CHECK(Same(All(grid, tolerance), All(Grid{curves, tolerance}, tolerance)));

    /**
     *  A grid of curves with no box has no cell, and is gridded once one has a box.
     */
    std::vector<Curve> loose = {Ray{.Start = {.X = 0, .Y = 0}, .Way = {.X = 1, .Y = 1}},
                                Xline{.Through = {.X = 0, .Y = 1}, .Way = {.X = 1, .Y = 0}}};
    Grid               endless{loose, tolerance};
    loose.at(1) = Line{.Start = {.X = -1, .Y = 0}, .End = {.X = 1, .Y = 0.5}};
    endless.Moved(1, loose.at(1));
    loose.at(0) = Line{.Start = {.X = 0, .Y = -1}, .End = {.X = 0.2, .Y = 1}};
    endless.Moved(0, loose.at(0));
    CHECK(Same(All(endless, tolerance), All(Grid{loose, tolerance}, tolerance)));
    CHECK(All(endless, tolerance).size() == 1);
    CHECK_THROWS_AS(endless.Moved(2, loose.at(0)), std::out_of_range);
  }

  TEST_CASE("grid: copied and taken over, a grid meets as it did")
  {
    Scatter                  scatter{0xC0B7};
    const Tolerance          tolerance = {.Length = 1e-8, .Radians = 1e-9};
    const std::vector<Curve> curves    = scatter.Curves(40, About(5), 2);
    Grid                     grid{curves, tolerance};
    grid.Moved(3, scatter.Bounded({.X = 1, .Y = 1}, 2));
    const std::vector<Found> found = All(grid, tolerance);
    REQUIRE(!found.empty());

    const Grid copy{grid};
    CHECK(Same(All(copy, tolerance), found));
    Grid assigned{{}, tolerance};
    assigned = grid;
    CHECK(Same(All(assigned, tolerance), found));
    CHECK(&assigned[0].Shape() != &grid[0].Shape());

    Grid taken{std::move(assigned)};
    CHECK(Same(All(taken, tolerance), found));
    Grid over{{}, tolerance};
    over = std::move(taken);
    CHECK(Same(All(over, tolerance), found));
    CHECK(over[0].Shape() == curves.front());

    /**
     *  A grid taken over is left with no curve, and meets in nothing.
     */
    const Curve    line = Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 1, .Y = 1}};
    const Prepared ready{line};
    CHECK(taken.Size() == 0);
    CHECK(All(taken, tolerance).empty());
    CHECK(All(taken, {.Length = 1, .Radians = 1e-9}).empty());
    CHECK(With({ready, taken}, tolerance).empty());
  }

  /**
   *  Pools of one, two, three, eight and sixty-four threads meet random grids of every kind and
   *  size, lines alone, floor plans and grids moved, within their own tolerance, a finer and a
   *  coarser one, as the calling thread meets them alone: the same pairs, marks and order, to the
   *  bit.
   */
  TEST_CASE("grid: met on a pool of any size, a grid meets as it does alone, to the bit")
  {
    std::array<Workers, 5> pools = {Workers{1}, Workers{2}, Workers{3}, Workers{8}, Workers{64}};
    Scatter                scatter{0xC4E3};
    int                    differ = 0;
    std::size_t            met    = 0;
    const auto agreed = [&pools, &met](const Grid &grid, const Tolerance &tolerance) -> bool
    {
      const std::vector<Found> alone  = Alone(grid, tolerance);
      met                            += alone.size();
      return std::ranges::all_of(pools, [&grid, &tolerance, &alone](Workers &pool) -> bool
                                 { return Same(Manned({grid, pool}, tolerance), alone); });
    };
    const int rounds = Rounds(24);
    for (int round = 0; round < rounds; ++round)
    {
      CAPTURE(round);
      const double       scale  = std::ldexp(1.0, static_cast<int>(scatter.Pick(40)) - 20);
      const Box          region = About(scale);
      const std::size_t  count  = 1 + scatter.Pick(120);
      std::vector<Curve> curves = scatter.Pick(3) == 0 ? scatter.Strokes(count, region, scale / 4)
                                                       : scatter.Curves(count, region, scale / 4);
      for (std::size_t loose = scatter.Pick(4); loose > 2; --loose)
      {
        curves.at(scatter.Pick(curves.size())) = scatter.Endless(scatter.Within(region));
      }
      const Tolerance tolerance = {.Length = scale * 1e-9, .Radians = 1e-9};
      Grid            grid{curves, tolerance};
      for (const double times : {1.0, 0.25, 64.0})
      {
        differ += agreed(grid, {.Length = tolerance.Length * times, .Radians = 1e-9}) ? 0 : 1;
      }
      grid.Moved(scatter.Pick(curves.size()), scatter.Bounded(scatter.Within(region), scale));
      differ += agreed(grid, tolerance) ? 0 : 1;
    }
    for (const std::size_t rooms : {std::size_t{1}, std::size_t{12}, std::size_t{300}})
    {
      std::vector<Curve> curves = Plan(rooms, rooms);
      for (std::size_t index = 5; index < curves.size(); index += 700)
      {
        curves.at(index) = scatter.Endless({.X = 5, .Y = 5});
      }
      const Grid grid{curves, Tolerance::Of(Spanned(curves))};
      differ += agreed(grid, Tolerance::Of(Spanned(curves))) ? 0 : 1;
    }
    const double             none = std::numeric_limits<double>::quiet_NaN();
    const Grid               empty{{}, {.Length = 1e-9, .Radians = 1e-9}};
    const std::vector<Curve> endless  = {scatter.Endless({.X = 0, .Y = 0}),
                                         scatter.Endless({.X = 1, .Y = 0}),
                                         scatter.Endless({.X = 0, .Y = 1})};
    differ                           += agreed(empty, {.Length = 1e-9, .Radians = 1e-9}) ? 0 : 1;
    differ +=
        agreed(Grid{endless, {.Length = 1e-9, .Radians = 1e-9}}, {.Length = 1e-9, .Radians = 1e-9})
            ? 0
            : 1;
    differ += agreed(Grid{Plan(4, 3), {.Length = none, .Radians = 1e-9}},
                     {.Length = 1e-9, .Radians = 1e-9})
                ? 0
                : 1;
    CHECK(differ == 0);
    CHECK(met > 1000);
    CHECK(pools.at(4).Size() == 64);
  }

  /**
   *  A plan large enough for IntersectAll to make a pool of its own meets as alone; one too small
   *  is met alone.
   */
  TEST_CASE(
      "grid: IntersectAll given no pool makes one where the grid is large, and gives the same")
  {
    const std::vector<Curve> curves    = Plan(512, 0x9A5);
    const Tolerance          tolerance = Tolerance::Of(Spanned(curves));
    const Grid               grid{curves, tolerance};
    const std::vector<Found> alone = Alone(grid, tolerance);
    CHECK(alone.size() > curves.size());
    CHECK(Same(All(grid, tolerance), alone));
  }

  /**
   *  A pool working one call meets another, made from within what the first tells or on another
   *  thread at once, alone on the calling thread: each gives what it gives alone.
   */
  TEST_CASE("grid: a pool already working leaves a call to meet alone, as it gives the same")
  {
    Workers                  pool{4};
    const std::vector<Curve> one       = Plan(30, 1);
    const std::vector<Curve> other     = Plan(20, 2);
    const Tolerance          tolerance = Tolerance::Of(Spanned(one));
    const Grid               first{one, tolerance};
    const Grid               second{other, tolerance};
    const std::vector<Found> firstAlone  = Alone(first, tolerance);
    const std::vector<Found> secondAlone = Alone(second, tolerance);

    std::vector<Found> outer;
    std::vector<Found> inner;
    IntersectAll({first, pool}, tolerance,
                 [&outer, &inner, &second, &pool, &tolerance](const Grid::Pair &curves,
                                                              const Marking    &marking) -> void
                 {
                   if (outer.empty())
                   {
                     inner = Manned({second, pool}, tolerance);
                   }
                   outer.push_back({.Curves = curves, .Marking = marking});
                 });
    CHECK(Same(outer, firstAlone));
    CHECK(Same(inner, secondAlone));

    std::vector<Found> there;
    std::vector<Found> here;
    {
      const std::jthread thread{[&there, &second, &pool, &tolerance]() -> void
                                { there = Manned({second, pool}, tolerance); }};
      here = Manned({first, pool}, tolerance);
    }
    CHECK(Same(here, firstAlone));
    CHECK(Same(there, secondAlone));
  }

  TEST_CASE("grid: what a sink throws on a pool is thrown on once the pool is done, and it meets "
            "again")
  {
    Workers                  pool{3};
    const std::vector<Curve> curves    = Plan(40, 5);
    const Tolerance          tolerance = Tolerance::Of(Spanned(curves));
    const Grid               grid{curves, tolerance};
    int                      told = 0;
    CHECK_THROWS_AS(IntersectAll({grid, pool}, tolerance,
                                 [&told](const Grid::Pair &, const Marking &) -> void
                                 {
                                   ++told;
                                   if (told == 3)
                                   {
                                     throw std::runtime_error{"enough"};
                                   }
                                 }),
                    std::runtime_error);
    CHECK(told == 3);
    CHECK(Same(Manned({grid, pool}, tolerance), Alone(grid, tolerance)));
  }

  /**
   *  Boxes of every size about random grids, a point, an empty box, one of a number that is not a
   *  number and the whole plane: the curves near each are those whose boxes grown overlap it and
   *  those with none, and a window and a crossing window pick what Inside and Intersects pick of
   *  every curve, before and after curves are moved.
   */
  TEST_CASE("grid: the curves near a box, and those a window and a crossing window pick")
  {
    Scatter      scatter{0x5E1EC7};
    int          differ   = 0;
    std::size_t  windowed = 0;
    std::size_t  crossed  = 0;
    const double none     = std::numeric_limits<double>::quiet_NaN();
    const double endless  = std::numeric_limits<double>::infinity();
    const auto picked = [&differ, &windowed, &crossed](const Grid &grid, const Tolerance &tolerance,
                                                       const Box &box) -> void
    {
      std::vector<std::size_t> near;
      std::vector<std::size_t> window;
      std::vector<std::size_t> crossing;
      for (std::size_t index = 0; index < grid.Size(); ++index)
      {
        const std::optional<Box> grown = Gridwork::Grown(grid, index, tolerance);
        if (!grown || Overlaps(*grown, box))
        {
          near.push_back(index);
        }
        if (Inside(grid[index].Shape(), box))
        {
          window.push_back(index);
        }
        if (Intersects(grid[index].Shape(), box))
        {
          crossing.push_back(index);
        }
      }
      differ   += grid.Near(box) == near ? 0 : 1;
      differ   += grid.Window(box) == window ? 0 : 1;
      differ   += grid.Crossing(box) == crossing ? 0 : 1;
      windowed += window.size();
      crossed  += crossing.size();
    };
    const int rounds = Rounds(32);
    for (int round = 0; round < rounds; ++round)
    {
      CAPTURE(round);
      const double       scale     = std::ldexp(1.0, static_cast<int>(scatter.Pick(40)) - 20);
      const Box          region    = About(scale);
      std::vector<Curve> curves    = scatter.Curves(10 + scatter.Pick(80), region, scale / 6);
      const Tolerance    tolerance = {.Length = scale * 1e-9, .Radians = 1e-9};
      curves.push_back(scatter.Endless(scatter.Within(region)));
      Grid grid{curves, tolerance};
      for (int box = 0; box < 8; ++box)
      {
        const Point  one = scatter.Within(region);
        const double width =
            scale * std::ldexp(scatter.Near(1, 0.5), -static_cast<int>(scatter.Pick(6)));
        const double height =
            scale * std::ldexp(scatter.Near(1, 0.5), -static_cast<int>(scatter.Pick(6)));
        picked(grid, tolerance, {.Min = one, .Max = {.X = one.X + width, .Y = one.Y + height}});
      }
      const Point spot = scatter.Within(region);
      picked(grid, tolerance, {.Min = spot, .Max = spot});
      picked(grid, tolerance, {});
      picked(grid, tolerance, {.Min = {.X = none, .Y = 0}, .Max = {.X = 1, .Y = 1}});
      picked(grid, tolerance,
             {.Min = {.X = -endless, .Y = -endless}, .Max = {.X = endless, .Y = endless}});
      grid.Moved(scatter.Pick(curves.size()),
                 scatter.Bounded(scatter.Within(About(4 * scale)), scale));
      picked(grid, tolerance, About(2 * scale));
    }
    const std::vector<Curve> ray = {Ray{.Start = {.X = 0, .Y = 0}, .Way = {.X = 1, .Y = 0}}};
    const Grid               endlessly{ray, {.Length = 1e-9, .Radians = 1e-9}};
    CHECK(endlessly.Near(About(1))     == std::vector<std::size_t>{0});
    CHECK(endlessly.Crossing(About(1)) == std::vector<std::size_t>{0});
    CHECK(endlessly.Window(About(1)).empty());
    CHECK(differ == 0);
    CHECK(windowed > 0);
    CHECK(crossed > windowed);
  }

  TEST_CASE("grid: the curves a polygon window and a polygon crossing window pick")
  {
    const Tolerance          tolerance = {.Length = 1e-9, .Radians = 1e-9};
    const std::vector<Curve> curves    = {
        Line{.Start = {.X = 3, .Y = 1}, .End = {.X = 7, .Y = 1}},
        Line{.Start = {.X = 3, .Y = 1}, .End = {.X = 12, .Y = 1}},
        Line{.Start = {.X = 11, .Y = 5}, .End = {.X = 12, .Y = 6}},
        Circle({.X = 5, .Y = 3}, 1),
        Circle({.X = 5, .Y = 3}, 20),
        Ray{.Start = {.X = 5, .Y = 3}, .Way = {.X = 1, .Y = 0}},
        Line{.Start = {.X = 50, .Y = 50}, .End = {.X = 51, .Y = 51}},
    };
    const Grid   grid(curves, tolerance);
    const Region triangle = {
        .Loops = {Loop{.Edges = {Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                                       {.Position = {.X = 10, .Y = 0}},
                                                       {.Position = {.X = 5, .Y = 8}}},
                                          .Closed   = true}}}}};
    CHECK(grid.Window(triangle)   == std::vector<std::size_t>{0, 3});
    CHECK(grid.Crossing(triangle) == std::vector<std::size_t>{0, 1, 3, 5});
  }
} // namespace Geoxy
