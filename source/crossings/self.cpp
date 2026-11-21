/***************************************************************************************************
 *  @file      source/crossings/self.cpp
 *  @brief     A curve's crossings with itself: cut into pieces that cannot cross themselves, the
 *             pieces met in pairs through a grid.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "crossings/self.hpp"

#include "real.hpp"
#include "scratch.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/grid.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/orient.hpp"
#include "geoxy/part.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/prepared.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <numeric>
#include <optional>
#include <span>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy::Crossings
{
  namespace
  {
    /**
     *  How many times a piece is halved at most before it is taken as it is.
     */
    constexpr int Deepest = 24;

    /**
     *  A piece of the curve: a curve of its own whose parameter, times a scale and plus a base, is
     *  the curve's.
     */
    struct Piece
    {
      Curve  Shape = {};  /**< The piece.                                */
      double Base  = 0.0; /**< The curve's parameter where its own is 0. */
      double Scale = 1.0; /**< The curve's parameter per its own.        */
    };

    /**
     *  The pieces of a curve, and whether one was taken as it is, too deep to halve.
     */
    struct Pieces
    {
      std::vector<Piece> All   = {};    /**< The pieces, in order along the curve. */
      bool               Spent = false; /**< Whether one was too deep to halve.    */
    };

    /**
     *  Whether numbers run one way, not all the same.
     *
     *  @param  values The numbers.
     *  @return        Whether none falls, or none rises, and they are not all one.
     */
    auto Monotone(std::span<const double> values) -> bool
    {
      const bool rising  = std::ranges::is_sorted(values);
      const bool falling = std::ranges::is_sorted(values, std::greater<>{});
      return rising != falling;
    }

    /**
     *  Whether control points run one way along x, along y or along the chord from the first to the
     *  last, so that the curve they make cannot cross itself; or are all one point.
     *
     *  @param  points The control points.
     *  @return        Whether they do.
     */
    auto Untangled(std::span<const Point> points) -> bool
    {
      const Vector        chord = points.back() - points.front();
      std::vector<double> across;
      std::vector<double> upward;
      std::vector<double> along;
      for (const Point point : points)
      {
        across.push_back(point.X);
        upward.push_back(point.Y);
        along.push_back(Dot(point - points.front(), chord));
      }
      const bool one = std::ranges::all_of(points, [&points](Point point) -> bool
                                           { return point == points.front(); });
      return one || Monotone(across) || Monotone(upward) || Monotone(along);
    }

    /**
     *  A piece's control points.
     *
     *  @tparam ShapeT A Bézier curve or a spline.
     *  @param  shape  The piece.
     *  @return        Its control points.
     */
    template <typename ShapeT>
    auto Controls(const ShapeT &shape) -> std::vector<Point>
    {
      if constexpr (std::is_same_v<ShapeT, Bezier::Quadratic>)
      {
        return {shape.Start, shape.Control, shape.End};
      }
      else if constexpr (std::is_same_v<ShapeT, Bezier::Cubic>)
      {
        return {shape.Start, shape.Out, shape.In, shape.End};
      }
      else
      {
        return shape.Poles;
      }
    }

    /**
     *  A Bézier curve's or a spline's stretch halved until each half's control points run one way,
     *  each half added to the pieces.
     *
     *  @tparam ShapeT  A Bézier curve, whose parts run from 0 to 1, or a spline, whose parts keep
     *                  its parameter.
     *  @param  shape   The whole curve.
     *  @param  stretch The stretch, and how many times it was halved.
     *  @param  pieces  The pieces, the halves added.
     */
    template <typename ShapeT>
    void Halved(const ShapeT &shape, std::pair<Range, int> stretch, Pieces &pieces)
    {
      const auto [range, depth] = stretch;
      const ShapeT part         = Part(shape, range);
      if (Untangled(Controls(part)) || depth == Deepest)
      {
        pieces.Spent = std::max(pieces.Spent, depth == Deepest);
        if constexpr (std::is_same_v<ShapeT, Spline>)
        {
          pieces.All.push_back({.Shape = part, .Base = 0.0, .Scale = 1.0});
        }
        else
        {
          pieces.All.push_back(
              {.Shape = part, .Base = range.Start, .Scale = range.End - range.Start});
        }
        return;
      }
      const double middle = (range.Start + range.End) / 2.0;
      Halved(shape, {{.Start = range.Start, .End = middle}, depth + 1}, pieces);
      Halved(shape, {{.Start = middle, .End = range.End}, depth + 1}, pieces);
    }

    /**
     *  The pieces of a curve that cannot cross themselves.
     *
     *  @param  curve The curve.
     *  @return       The pieces; none for a line, a ray, an xline, or what has no segment or span.
     */
    auto Cut(const Curve &curve) -> Pieces
    {
      Pieces pieces;
      std::visit(
          [&pieces](const auto &shape) -> void
          {
            using ShapeT = std::decay_t<decltype(shape)>;
            if constexpr (std::is_same_v<ShapeT, Polyline>)
            {
              for (std::size_t index = 0; index < Count(shape); ++index)
              {
                pieces.All.push_back({.Shape = std::visit([](const auto &segment) -> Curve
                                                          { return segment; }, Nth(shape, index)),
                                      .Base  = static_cast<double>(index),
                                      .Scale = 1.0});
              }
            }
            else if constexpr (std::is_same_v<ShapeT, Spline>)
            {
              if (!Valid(shape))
              {
                return;
              }
              for (auto span = static_cast<std::size_t>(shape.Degree); span < shape.Poles.size();
                   ++span)
              {
                const Range range = {.Start = shape.Knots[span], .End = shape.Knots[span + 1]};
                if (range.Start < range.End)
                {
                  Halved(shape, {range, 0}, pieces);
                }
              }
            }
            else if constexpr (std::is_same_v<ShapeT, Bezier::Quadratic>
                               || std::is_same_v<ShapeT, Bezier::Cubic>)
            {
              Halved(shape, {{.Start = 0.0, .End = 1.0}, 0}, pieces);
            }
            else if constexpr (std::is_same_v<ShapeT, Arc> || std::is_same_v<ShapeT, Ellipse>)
            {
              const double turns = std::abs(shape.Sweep.InRadians()) / Turn.InRadians();
              const auto   count = static_cast<std::size_t>(std::max(1.0, std::ceil(turns)));
              for (std::size_t index = 0; index < count; ++index)
              {
                const auto  whole = static_cast<double>(count);
                const Range range = {.Start = static_cast<double>(index) / whole,
                                     .End   = static_cast<double>(index + 1) / whole};
                pieces.All.push_back({.Shape = Part(shape, range),
                                      .Base  = range.Start,
                                      .Scale = range.End - range.Start});
              }
            }
          },
          curve);
      return pieces;
    }

    /**
     *  Where each piece starts in the curve's parameter.
     *
     *  @param  pieces The pieces.
     *  @return        Their starts, rising, and the last one's end after them.
     */
    auto Starts(std::span<const Piece> pieces) -> std::vector<double>
    {
      std::vector<double> starts;
      starts.reserve(pieces.size() + 1);
      for (const Piece &piece : pieces)
      {
        starts.push_back(piece.Base + (Domain(piece.Shape).Start * piece.Scale));
      }
      const Piece &last = pieces.back();
      starts.push_back(last.Base + (Domain(last.Shape).End * last.Scale));
      return starts;
    }

    /**
     *  The curve cut, and where its pieces start, and whether it closes.
     */
    struct Course
    {
      const Pieces       *Cut    = nullptr; /**< The pieces.                     */
      std::vector<double> Starts = {};      /**< Where each starts, and the end. */
      bool                Closed = false;   /**< Whether its end is its start.   */
    };

    /**
     *  Where a pass of the curve through a place comes from and goes to, as ways out of it: back
     *  along the way it came in, and on along the way it leaves, from the pieces either side of the
     *  place, or the one it is within.
     *
     *  @param  course The pieces, where they start, and whether the curve closes.
     *  @param  along  The place, in the curve's parameter; a closed curve's end taken as its start.
     *  @return        The two ways; none at an open curve's own end.
     */
    auto Ways(const Course &course, double along) -> std::optional<std::pair<Vector, Vector>>
    {
      const std::vector<Piece> &pieces = course.Cut->All;
      const std::size_t         count  = pieces.size();
      if (Exactly(along, course.Starts.back()))
      {
        return std::nullopt;
      }
      const auto   after = std::ranges::upper_bound(course.Starts, along);
      const auto   index = static_cast<std::size_t>(after - course.Starts.begin()) - 1;
      const Piece &piece = pieces[index];
      const double local = (along - piece.Base) / piece.Scale;
      if (!Exactly(along, course.Starts[index]))
      {
        return std::pair{-Derivative(piece.Shape, local), Derivative(piece.Shape, local)};
      }
      if (index == 0 && !course.Closed)
      {
        return std::nullopt;
      }
      const Piece &before = pieces[(index + count - 1) % count];
      return std::pair{-Derivative(before.Shape, Domain(before.Shape).End),
                       Derivative(piece.Shape, Domain(piece.Shape).Start)};
    }

    /**
     *  Whether a way lies strictly within the counterclockwise sweep from one way to another:
     *  within a half turn after the first and before the last where the sweep is under a half turn,
     *  within either where it is not.
     *
     *  @param  from The sweep's first way.
     *  @param  last Its last.
     *  @param  way  The way.
     *  @return      Whether it does.
     */
    auto Swept(Vector from, Vector last, Vector way) -> bool
    {
      const Point origin = {};
      const bool  narrow =
          Orient(origin, origin + from, origin + last) == Orientation::Counterclockwise;
      const int sides = static_cast<int>(Orient(origin, origin + from, origin + way)
                                         == Orientation::Counterclockwise)
                      + static_cast<int>(Orient(origin, origin + way, origin + last)
                                         == Orientation::Counterclockwise);
      return narrow ? sides == 2 : sides > 0;
    }

    /**
     *  Whether one pass's ways lie either side of another's: a crossing, not a touch.
     *
     *  @param  first  The first pass's ways.
     *  @param  second The second's.
     *  @return        Whether they do; not where two ways run along one another.
     */
    auto Separated(const std::pair<Vector, Vector> &first, const std::pair<Vector, Vector> &second)
        -> bool
    {
      const Point origin = {};
      for (const Vector way : {first.first, first.second})
      {
        for (const Vector other : {second.first, second.second})
        {
          if (Orient(origin, origin + way, origin + other) == Orientation::Collinear
              && Dot(way, other)                           > 0.0)
          {
            return false;
          }
        }
      }
      return Swept(first.first, first.second, second.first)
          != Swept(first.first, first.second, second.second);
    }

    /**
     *  A piece's place as the curve's, a closed curve's end taken as its start.
     *
     *  @param  course The pieces, and whether the curve closes.
     *  @param  piece  The piece's index.
     *  @param  place  The place on it.
     *  @return        The place on the curve.
     */
    auto Transformed(const Course &course, std::size_t piece, Place place) -> Place
    {
      const Piece &held  = course.Cut->All[piece];
      const double along = held.Base + (place.Along * held.Scale);
      const bool   wraps = course.Closed && Exactly(along, course.Starts.back());
      return {.At = place.At, .Along = wraps ? course.Starts.front() : along};
    }

    /**
     *  The hits and stretches two pieces share, as the curve's: the end two pieces that follow one
     *  another share left out, each hit's earlier place first.
     *
     *  @param  course    The pieces, and whether the curve closes.
     *  @param  pair      The two pieces, and their meeting.
     *  @param  tolerance Within which two points are one.
     *  @return           What they share.
     */
    auto Shared(const Course &course, std::pair<Grid::Pair, const Intersections *> pair,
                const Tolerance &tolerance) -> Intersections
    {
      const auto [first, second]     = pair.first;
      const std::vector<Piece> &all  = course.Cut->All;
      const bool                next = second == first + 1;
      const bool    adjacent = next || (course.Closed && first == 0 && second + 1 == all.size());
      const Curve  &later    = next ? all[second].Shape : all[first].Shape;
      const Point   joint    = At(later, Domain(later).Start);
      Intersections shared   = {.Hits = {}, .Overlaps = {}, .Spent = pair.second->Spent};
      for (const Hit &hit : pair.second->Hits)
      {
        if (adjacent && hit.Kind == Contact::End && Same(hit.At, joint, tolerance))
        {
          continue;
        }
        const Place one = Transformed(course, first, hit.First);
        const Place two = Transformed(course, second, hit.Second);
        shared.Hits.push_back({.At     = hit.At,
                               .First  = one.Along < two.Along ? one : two,
                               .Second = one.Along < two.Along ? two : one,
                               .Kind   = hit.Kind});
      }
      for (const Overlap &overlap : pair.second->Overlaps)
      {
        const Piece &one = all[first];
        const Piece &two = all[second];
        shared.Overlaps.push_back(
            {.First  = {.Start = one.Base + (overlap.First.Start * one.Scale),
                        .End   = one.Base + (overlap.First.End * one.Scale)},
             .Second = {.Start = two.Base + (overlap.Second.Start * two.Scale),
                        .End   = two.Base + (overlap.Second.End * two.Scale)},
             .Along  = overlap.Along});
      }
      return shared;
    }

    /**
     *  Hits at one place made one, in order along the curve, and those Intersect told only as ends
     *  judged a crossing or a touch by the ways the curve passes there.
     *
     *  @param  course The pieces, and whether the curve closes.
     *  @param  hits   The hits.
     *  @return        The hits merged.
     */
    auto Merged(const Course &course, const std::vector<Hit> &hits) -> std::vector<Hit>
    {
      /**
       *  The hits' indices sorted by their places and then by index, the order a stable sort gives,
       *  with no room taken but the indices' own.
       */
      List<std::size_t> order(hits.size());
      std::ranges::iota(order, std::size_t{0});
      std::ranges::sort(order,
                        [&hits](std::size_t lhs, std::size_t rhs) -> bool
                        {
                          return std::tuple{hits[lhs].First.Along, hits[lhs].Second.Along, lhs}
                               < std::tuple{hits[rhs].First.Along, hits[rhs].Second.Along, rhs};
                        });
      std::vector<Hit> merged;
      for (const std::size_t index : order)
      {
        Hit        hit   = hits[index];
        const bool again = !merged.empty()
                        && Exactly(merged.back().First.Along, hit.First.Along)
                        && Exactly(merged.back().Second.Along, hit.Second.Along);
        if (again)
        {
          continue;
        }
        const auto first  = Ways(course, hit.First.Along);
        const auto second = Ways(course, hit.Second.Along);
        if (hit.Kind == Contact::End && first && second)
        {
          hit.Kind = Separated(*first, *second) ? Contact::Cross : Contact::Touch;
        }
        merged.push_back(hit);
      }
      return merged;
    }
  } // namespace

  auto Self(const Curve &curve, const Tolerance &tolerance) -> Intersections
  {
    const Pieces       cut = Cut(curve);
    std::vector<Curve> shapes;
    shapes.reserve(cut.All.size());
    for (const Piece &piece : cut.All)
    {
      shapes.push_back(piece.Shape);
    }
    Intersections meeting = {.Hits = {}, .Overlaps = {}, .Spent = cut.Spent};
    if (shapes.size() < 2)
    {
      return meeting;
    }
    const Range      domain = Domain(curve);
    const Course     course = {.Cut    = &cut,
                               .Starts = Starts(cut.All),
                               .Closed =
                                   Same(At(curve, domain.Start), At(curve, domain.End), tolerance)};
    const Grid       grid(shapes, tolerance);
    std::vector<Hit> hits;
    IntersectAll(grid, tolerance,
                 [&](const Grid::Pair &pair, const Marking & /*marking*/) -> void
                 {
                   const auto [first, second] = std::minmax(pair[0], pair[1]);
                   Intersections met;
                   Intersect(Prepared::Pair{grid[first], grid[second]}, tolerance, met);
                   const Intersections shared = Shared(course, {{first, second}, &met}, tolerance);
                   meeting.Spent              = std::max(meeting.Spent, shared.Spent);
                   hits.insert(hits.end(), shared.Hits.begin(), shared.Hits.end());
                   meeting.Overlaps.insert(meeting.Overlaps.end(), shared.Overlaps.begin(),
                                           shared.Overlaps.end());
                 });
    meeting.Hits = Merged(course, hits);
    std::ranges::sort(meeting.Overlaps, [](const Overlap &lhs, const Overlap &rhs) -> bool
                      { return lhs.First.Start < rhs.First.Start; });
    return meeting;
  }
} // namespace Geoxy::Crossings
