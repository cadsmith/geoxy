/***************************************************************************************************
 *  @file      source/swift/bulge.cpp
 *  @brief     A line, a ray or an endless line met with a polyline: each straight segment by its
 *             cross products, each bulge by its exact rational circle, every crossing worked in
 *             double-double and proved the double nearest the true point.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "real.hpp"
#include "straight.hpp"
#include "sure.hpp"
#include "swift/lane.hpp"

#include "geoxy/domain.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <utility>

namespace Geoxy::Swift
{
  namespace
  {
    using Sure::Number;

    /**
     *  A segment of a polyline as the kernel takes it: its vertices and its bulge, and where it is.
     */
    struct Leg
    {
      Point       Start = {};  /**< Its first vertex, P0.              */
      Point       End   = {};  /**< Its second vertex, P1.             */
      double      Bend  = 0.0; /**< Its bulge, b; zero where straight. */
      std::size_t Index = 0;   /**< Which segment of the polyline.     */
    };

    /**
     *  How far a point is from a bulge's circle, at the least. The circle is the points X where
     *  G(X) = 4 b (X - P0).(X - P1) - 2 (1 - b^2) (X - P0).n is zero, n the chord c turned a
     *  quarter, and G(X) = 4 b (|X - C|^2 - r^2) for its center C and radius r = |c| (1 + b^2) / (4
     *  |b|): the distance | |X - C| - r | is |G| / (4 |b| (|X - C| + r)), and |X - C| + r is no
     *  more than |X - P0| + 2 r. G's few roundings are each within 2^-50 of its terms.
     *
     *  @param  leg   The bulge.
     *  @param  point The point.
     *  @return       The distance, at the least.
     */
    auto Bowed(const Leg &leg, Point point) -> double
    {
      const Vector from   = point - leg.Start;
      const Vector toEnd  = point - leg.End;
      const Vector chord  = leg.End - leg.Start;
      const Vector across = Perpendicular(chord);
      const double bend   = leg.Bend;
      const double rest   = 1.0 - (bend * bend);
      const double value  = (4.0 * bend * Dot(from, toEnd)) - (2.0 * rest * Dot(from, across));
      const double slack =
          0x1p-46
          * ((4.0 * std::abs(bend) * (std::abs(from.X * toEnd.X) + std::abs(from.Y * toEnd.Y)))
             + (2.0
                * std::abs(rest)
                * (std::abs(from.X * across.X) + std::abs(from.Y * across.Y))));
      const double below =
          ((4.0 * std::abs(bend) * Length(from)) + (2.0 * Length(chord) * (1.0 + (bend * bend))))
          * (1.0 + 0x1p-46);
      return (std::abs(value) - slack) / below;
    }

    /**
     *  How far a point is from a straight segment, at the least: from its carrier, or past an end
     *  from that end, as Away weighs a lane.
     *
     *  @param  leg   The segment.
     *  @param  point The point.
     *  @return       The distance, at the least.
     */
    auto Straightly(const Leg &leg, Point point) -> double
    {
      const Vector way   = leg.End - leg.Start;
      const Vector gap   = point - leg.Start;
      const double speed = Length(way);
      const double slack = 0x1p-50 * (std::abs(gap.X * way.Y) + std::abs(gap.Y * way.X));
      const double off   = (std::abs(Cross(gap, way)) - slack) / (speed * (1.0 + 0x1p-50));
      const double along = Dot(gap, way) / (speed * speed);
      if (along < -0x1p-40)
      {
        return std::max(off, Distance(point, leg.Start) * (1.0 - 0x1p-50));
      }
      if (along > 1.0 + 0x1p-40)
      {
        return std::max(off, Distance(point, leg.End) * (1.0 - 0x1p-50));
      }
      return off;
    }

    /**
     *  A segment of a polyline: straight where its bulge is below Straight or its ends are one.
     *
     *  @param  polyline The polyline.
     *  @param  index    Which segment.
     *  @return          The segment.
     */
    auto Legging(const Polyline &polyline, std::size_t index) -> Leg
    {
      const Vertex &start = polyline.Vertices[index];
      const Vertex &end   = polyline.Vertices[(index + 1) % polyline.Vertices.size()];
      const bool    plain = !(std::abs(start.Bulge) >= Straight) || start.Position == end.Position;
      return {.Start = start.Position,
              .End   = end.Position,
              .Bend  = plain ? 0.0 : start.Bulge,
              .Index = index};
    }

    /**
     *  What meeting a lane with a polyline's segments gathers: the lane, the tolerance, and the
     *  crossings found.
     */
    struct Harvest
    {
      const Lane      *Straight  = nullptr; /**< The lane.            */
      Geoxy::Tolerance Tolerance = {};      /**< The tolerance.       */
      Spots            Found     = {};      /**< The crossings found. */
    };

    /**
     *  Which side of a lane's carrier a point lies on, where rounding leaves it clear: the cross
     *  product of its gap from P with d, d rounded once at the most and each operation once, within
     *  2^-50 of its terms.
     *
     *  @param  lane  The lane.
     *  @param  point The point.
     *  @return       1 to the right, -1 to the left, 0 where it is not clear.
     */
    auto Sided(const Lane &lane, Point point) -> int
    {
      const Vector gap   = point - lane.Start;
      const double cross = Cross(gap, lane.Rough);
      const double slack =
          0x1p-50 * (std::abs(gap.X * lane.Rough.Y) + std::abs(gap.Y * lane.Rough.X));
      return static_cast<int>(cross > slack) - static_cast<int>(cross < -slack);
    }

    /**
     *  A crossing's point and the lane's place there, P + t d with each coordinate the double
     *  nearest the true one, a crossing that is a spot of the polyline's segment.
     *
     *  @param  lane  The lane.
     *  @param  along The lane's place, t, with its bound.
     *  @param  index The segment.
     *  @return       The spot; nothing where a coordinate is not settled.
     */
    auto Spotted(const Lane &lane, const Number &along, std::size_t index) -> std::optional<Spot>
    {
      const auto point =
          Sure::Rounded(std::array{Sure::Exact(lane.Start.X) + (lane.Way[0] * along),
                                   Sure::Exact(lane.Start.Y) + (lane.Way[1] * along)});
      return point.transform(
          [&along, index](const std::array<double, 2> &rounded) -> Spot
          {
            return {.At    = {.X = rounded[0], .Y = rounded[1]},
                    .Along = {along.Value.High, static_cast<double>(index) + 0.5},
                    .Own   = false};
          });
    }

    /**
     *  A lane and a straight segment: where the segment's ends lie on either side of the lane's
     *  carrier, their crossing, by the cross products of their ways and of the gap between them,
     *  exact in two doubles; else none. A crossing within twice the reach over its sine of the
     *  lane's end, or not steep, or at a slant at which the thorough way may take the two to lie
     *  along each other, is not sure.
     *
     *  @param  leg       The segment.
     *  @param  gathering The lane, the tolerance, and where its crossing goes.
     *  @return           Whether it was decided.
     */
    auto OnSegment(const Leg &leg, Harvest &gathering) -> bool
    {
      const Lane      &lane      = *gathering.Straight;
      const Tolerance &tolerance = gathering.Tolerance;
      const int        sides     = Sided(lane, leg.Start) * Sided(lane, leg.End);
      if (sides == 0)
      {
        return false;
      }
      if (sides > 0)
      {
        return true;
      }
      const Number wayX  = Sure::Apart(leg.End.X, leg.Start.X);
      const Number wayY  = Sure::Apart(leg.End.Y, leg.Start.Y);
      const Number gapX  = Sure::Apart(leg.Start.X, lane.Start.X);
      const Number gapY  = Sure::Apart(leg.Start.Y, lane.Start.Y);
      const Number turn  = (lane.Way[0] * wayY) - (lane.Way[1] * wayX);
      const Number along = Sure::Quotient((gapX * wayY) - (gapY * wayX), turn);
      const Vector way   = leg.End - leg.Start;
      const double sine = std::max(Sure::Least(turn) / (lane.Speed * Length(way) * (1.0 + 0x1p-48)),
                                   std::numeric_limits<double>::min());
      const double margin = (2.0 * Reach(tolerance) / (lane.Speed * sine * Sure::Shrink))
                          + (along.Bound * Sure::Grow)
                          + (0x1p-50 * std::abs(along.Value.High));
      const Side   place  = Within(lane.Along, along.Value.High, margin);
      if (place == Side::Out)
      {
        return true;
      }
      const Run run  = {.Start = lane.Start, .Way = lane.Rough, .Along = lane.Along};
      const Run side = {.Start = leg.Start, .Way = way, .Along = {.Start = 0.0, .End = 1.0}};
      if (place == Side::Near
          || !Steep(sine, tolerance)
          || !Unaligned(run, side, tolerance)
          || !Unaligned(side, run, tolerance))
      {
        return false;
      }
      const std::optional<Spot> spot = Spotted(lane, along, leg.Index);
      if (spot)
      {
        gathering.Found.push_back(*spot);
      }
      return spot.has_value();
    }

    /**
     *  A lane and a bulge, by the bulge's exact circle: the line's quadratic on it times 4 b, a t^2
     *  + 2 c t + e with a = 4 b d.d, c = 2 b d.(h + j) - (1 - b^2) d.n and e = 4 b h.j - 2 (1 -
     *  b^2) h.n, h = P - P0, j = P - P1, n the chord turned a quarter to the left; its discriminant
     *  D = c^2 - a e = 16 b^2 |d|^2 (r^2 - s^2), s the carrier's distance from the center. Where it
     *  cuts the circle, the shallower cap is r - s = (r^2 - s^2) / (r + s) >= D / (8 |b| |d|^2 |c|
     *  (1 + b^2)) high, and the crossings' sine with the circle is sqrt(D) / (|d| |c| (1 + b^2)). A
     *  crossing lies on the arc where it lies on the bulge's side of the chord.
     *
     *  @param  leg       The bulge.
     *  @param  gathering The lane, the tolerance, and where its crossings go.
     *  @return           Whether it was decided.
     */
    auto OnBulge(const Leg &leg, Harvest &gathering) -> bool
    {
      const Lane      &lane      = *gathering.Straight;
      const Tolerance &tolerance = gathering.Tolerance;
      const double     reach     = Reach(tolerance);
      const Number     dx        = lane.Way[0];
      const Number     dy        = lane.Way[1];
      const Number     fromX     = Sure::Apart(lane.Start.X, leg.Start.X);
      const Number     fromY     = Sure::Apart(lane.Start.Y, leg.Start.Y);
      const Number     pastX     = Sure::Apart(lane.Start.X, leg.End.X);
      const Number     pastY     = Sure::Apart(lane.Start.Y, leg.End.Y);
      const Number     chordX    = Sure::Apart(leg.End.X, leg.Start.X);
      const Number     chordY    = Sure::Apart(leg.End.Y, leg.Start.Y);
      const Number     acrossX   = -chordY;
      const Number     acrossY   = chordX;
      const double     bend      = leg.Bend;
      const Number     rest      = Sure::Exact(1.0) - Sure::Product(bend, bend);
      const Number     square    = (dx * dx) + (dy * dy);
      const Number     a         = Sure::Exact(4.0 * bend) * square;
      const Number c = (Sure::Exact(2.0 * bend) * ((dx * (fromX + pastX)) + (dy * (fromY + pastY))))
                     - (rest * ((dx * acrossX) + (dy * acrossY)));
      const Number constant = (Sure::Exact(4.0 * bend) * ((fromX * pastX) + (fromY * pastY)))
                            - (Sure::Exact(2.0) * rest * ((fromX * acrossX) + (fromY * acrossY)));
      const Number disc     = (c * c) - (a * constant);
      const double chord    = Length(leg.End - leg.Start) * (1.0 + 0x1p-50);
      const double size     = 1.0 + (bend * bend);
      if (disc.Value.High < 0.0)
      {
        /**
         *  The carrier misses the circle: by s - r = (s^2 - r^2) / (s + r), s + r no more than |h|
         *  + 2 r.
         */
        const double from = Length(lane.Start - leg.Start) * (1.0 + 0x1p-50);
        const double gap  = Sure::Least(disc)
                          / (4.0
                             * std::abs(bend)
                             * Sure::Most(square)
                             * ((4.0 * std::abs(bend) * from) + (2.0 * chord * size)))
                          * Sure::Shrink;
        return gap > reach;
      }
      const Number root  = Sure::Root(disc);
      const double speed = std::sqrt(Sure::Most(square)) * Sure::Grow;
      const double cap   = Sure::Least(disc)
                         / (8.0 * std::abs(bend) * Sure::Most(square) * chord * size)
                         * Sure::Shrink;
      const double sine  = std::max(
          std::sqrt(std::max(Sure::Least(disc), 0.0)) / (speed * chord * size) * Sure::Shrink,
          std::numeric_limits<double>::min());
      const Number skew = (chordX * fromY) - (chordY * fromX);
      const Number lean = (chordX * dy) - (chordY * dx);
      for (const double sign : {-1.0, 1.0})
      {
        const Number t      = Sure::Quotient(sign < 0.0 ? -c - root : -c + root, a);
        const double margin = (2.0 * reach / (lane.Speed * sine * Sure::Shrink))
                            + (t.Bound * Sure::Grow)
                            + (0x1p-50 * std::abs(t.Value.High));
        const Side   place  = Within(lane.Along, t.Value.High, margin);
        if (place == Side::Out)
        {
          continue;
        }
        /**
         *  The side of the chord: c x (X - P0) = c x h + t (c x d), the arc right of the chord for
         *  a bulge above zero; near the chord's line, a crossing is near a vertex.
         */
        const Number side = skew + (t * lean);
        if (place == Side::Near || !(Sure::Least(side) > 0.0))
        {
          return false;
        }
        if ((side.Value.High < 0.0) != (bend > 0.0))
        {
          continue;
        }
        const std::optional<Spot> spot =
            cap > reach && Steep(sine, tolerance) ? Spotted(lane, t, leg.Index) : std::nullopt;
        if (!spot)
        {
          return false;
        }
        gathering.Found.push_back(*spot);
      }
      return true;
    }
  } // namespace

  auto Shunned(const Polyline &polyline, Point point) -> double
  {
    double least = std::numeric_limits<double>::infinity();
    for (std::size_t index = 0; index < Count(polyline); ++index)
    {
      const Leg segment = Legging(polyline, index);
      if (segment.Start != segment.End)
      {
        least = std::min(least,
                         Zero(segment.Bend) ? Straightly(segment, point) : Bowed(segment, point));
      }
    }
    return least;
  }

  auto Bulged(const Lane &lane, const Polyline &polyline, const Tolerance &tolerance)
      -> std::optional<Spots>
  {
    const std::size_t count = Count(polyline);
    const double      reach = Reach(tolerance);
    const bool        held =
        std::ranges::all_of(polyline.Vertices, [](const Vertex &vertex) -> bool
                            { return Held({vertex.Position.X, vertex.Position.Y, vertex.Bulge}); });
    if (count == 0 || !held)
    {
      return std::nullopt;
    }
    Harvest    gathering = {.Straight = &lane, .Tolerance = tolerance, .Found = {}};
    const auto leg       = [&polyline](std::size_t index) -> Leg
    {
      return Legging(polyline, index);
    };
    /**
     *  Room for two crossings a segment, taken once: the arena is asked once a call at the most.
     */
    gathering.Found.reserve(2 * count);
    for (std::size_t index = 0; index < count; ++index)
    {
      const Leg segment = leg(index);
      if (segment.Start == segment.End)
      {
        continue;
      }
      if (!(Zero(segment.Bend) ? OnSegment(segment, gathering) : OnBulge(segment, gathering)))
      {
        return std::nullopt;
      }
    }
    /**
     *  Every vertex far from the lane, and every end of the lane far from every segment.
     */
    const bool vertices =
        std::ranges::all_of(polyline.Vertices, [&lane, reach](const Vertex &vertex) -> bool
                            { return Away(lane, vertex.Position) > reach; });
    const auto far = [&leg, count, reach](Point end) -> bool
    {
      for (std::size_t index = 0; index < count; ++index)
      {
        const Leg segment = leg(index);
        if (segment.Start != segment.End
            && !((Zero(segment.Bend) ? Straightly(segment, end) : Bowed(segment, end)) > reach))
        {
          return false;
        }
      }
      return true;
    };
    const bool ends = std::ranges::all_of(std::span{lane.Ends}.first(lane.Count), far);
    if (!vertices || !ends)
    {
      return std::nullopt;
    }
    return std::move(gathering.Found);
  }
} // namespace Geoxy::Swift
