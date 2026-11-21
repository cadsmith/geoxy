/***************************************************************************************************
 *  @file      source/swift/rim.cpp
 *  @brief     The circle an arc or an ellipse lies on, as its kernel sees it: its frame, its
 *             quarters and the side of it a point is on.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "swift/rim.hpp"

#include "elliptic.hpp"
#include "real.hpp"
#include "sure.hpp"
#include "swift/lane.hpp"
#include "swift/ring.hpp"
#include "swift/terms.hpp"
#include "twofold.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <optional>
#include <utility>
#include <variant>

namespace Geoxy::Swift
{
  auto Rimmed(const Curve &curve) -> std::optional<Ring>
  {
    const auto *arc = std::get_if<Arc>(&curve);
    if (arc != nullptr && !(arc->Radius > 0.0))
    {
      return std::nullopt;
    }
    const Ellipse shape = arc != nullptr ? Elliptic(*arc) : std::get<Ellipse>(curve);
    const Vector  axis  = shape.Map.X;
    const Vector  other = shape.Map.Y;
    const Point   home  = shape.Map.Origin;
    const double  sweep = shape.Sweep.InRadians();
    if (!Held({axis.X, axis.Y, other.X, other.Y, home.X, home.Y})
        || !(std::abs(Cross(axis, other)) > Narrowest)
        || Zero(sweep)
        || Coiled(shape))
    {
      return std::nullopt;
    }
    const Number det       = Sure::Product(axis.X, other.Y) - Sure::Product(axis.Y, other.X);
    const double frobenius = std::sqrt(Dot(axis, axis) + Dot(other, other)) * (1.0 + 0x1p-50);
    return Ring{.Map    = shape.Map,
                .Rows   = {Vector{.X = other.Y, .Y = -other.X}, Vector{.X = -axis.Y, .Y = axis.X}},
                .Level  = det * det,
                .High   = Sure::Most(det),
                .Low    = Sure::Least(det),
                .Least  = Shrink / frobenius,
                .Most   = frobenius,
                .Facing = det.Value.High > 0.0 ? 1.0 : -1.0,
                .From   = Direction(shape.Start),
                .To     = Direction(shape.Start + shape.Sweep),
                .Start  = shape.Start.InRadians(),
                .Sweep  = sweep,
                .Whole  = !(std::abs(sweep) < Turn.InRadians()),
                .Shape  = &curve,
                .Hub    = {Sure::Exact(home.X), Sure::Exact(home.Y)}};
  }

  auto Encircled(const Hump &hump) -> std::optional<Ring>
  {
    const Point  start   = hump.Start;
    const Point  end     = hump.End;
    const double bend    = hump.Bend;
    const Number across  = Sure::Apart(end.X, start.X);
    const Number upward  = Sure::Apart(end.Y, start.Y);
    const Number square  = Sure::Product(bend, bend);
    const Number rest    = Sure::Exact(1.0) - square;
    const Number heavy   = Sure::Exact(1.0) + square;
    const Number lean    = Sure::Exact(4.0 * bend);
    const Number middleX = (Sure::Exact(start.X) + Sure::Exact(end.X)) * Sure::Exact(0.5);
    const Number middleY = (Sure::Exact(start.Y) + Sure::Exact(end.Y)) * Sure::Exact(0.5);
    const Number centerX = middleX - Sure::Quotient(upward * rest, lean);
    const Number centerY = middleY + Sure::Quotient(across * rest, lean);
    const Number level   = Sure::Quotient(((across * across) + (upward * upward)) * heavy * heavy,
                                          Sure::Exact(16.0) * square);
    const double drift   = (std::max(std::abs(centerX.Value.Low) + centerX.Bound,
                                     std::abs(centerY.Value.Low) + centerY.Bound)
                            * Grow)
                         + Sure::Tiny;
    const double high    = std::sqrt(Sure::Most(level)) * Grow;
    const double low     = std::sqrt(Sure::Least(level)) * Shrink;
    const Point  center  = {.X = centerX.Value.High, .Y = centerY.Value.High};
    if (!Held({center.X, center.Y, high}) || !(low > 0.0))
    {
      return std::nullopt;
    }
    const double chord = Length(end - start);
    return Ring{
        .Map    = {.X = {.X = high, .Y = 0.0}, .Y = {.X = 0.0, .Y = high}, .Origin = center},
        .Rows   = {Vector{.X = 1.0, .Y = 0.0}, Vector{.X = 0.0, .Y = 1.0}},
        .Level  = level,
        .High   = high,
        .Low    = low,
        .Least  = Shrink,
        .Most   = Grow,
        .Facing = 1.0,
        .From   = {},
        .To     = {},
        .Start  = 0.0,
        .Sweep  = 0.0,
        .Whole  = false,
        .Shape  = nullptr,
        .Hub    = {centerX, centerY},
        .Drift  = drift,
        .Stride =
            chord * (1.0 + (bend * bend)) * Above(2.0 * static_cast<double>(hump.Index + 1)) * Grow,
        .Chord = &hump};
  }

  auto Seen(const Ring &ring, Point point) -> std::pair<Vector, double>
  {
    const Vector gap   = point - ring.Map.Origin;
    const Vector first = ring.Rows[0];
    const Vector other = ring.Rows[1];
    const double slack =
        (0x1p-50
         * (std::abs(first.X * gap.X)
            + std::abs(first.Y * gap.Y)
            + std::abs(other.X * gap.X)
            + std::abs(other.Y * gap.Y)))
        + (ring.Drift
           * (std::abs(first.X) + std::abs(first.Y) + std::abs(other.X) + std::abs(other.Y)));
    return {{.X = Dot(first, gap), .Y = Dot(other, gap)}, slack};
  }

  auto Hubbed(const Ring &ring, std::size_t axis) -> Number
  {
    return ring.Hub.at(axis);
  }

  auto Sidled(const Hump &hump, const std::array<Number, 2> &crossing, double margin) -> Side
  {
    const Vector chord = hump.End - hump.Start;
    const Vector gap = Point{.X = crossing[0].Value.High, .Y = crossing[1].Value.High} - hump.Start;
    const double cross = Cross(chord, gap);
    const double clear = (0x1p-48 * (std::abs(chord.X * gap.Y) + std::abs(chord.Y * gap.X)))
                       + (margin * Length(chord) * Grow);
    if (!(std::abs(cross) > clear))
    {
      return Side::Near;
    }
    return (cross < 0.0) == (hump.Bend > 0.0) ? Side::In : Side::Out;
  }

  auto RingGap(const Ring &ring, Point point) -> double
  {
    const auto [seen, slack] = Seen(ring, point);
    const double size        = Length(seen);
    const double off         = std::max(size - ring.High, ring.Low - size);
    return ring.Least * (off - (2.0 * slack) - (0x1p-50 * size)) * Shrink;
  }

  auto Swept(const Ring &ring, const std::pair<Vector, double> &seen, double angle) -> Side
  {
    const bool   forward = ring.Sweep > 0.0;
    const Vector start   = forward ? ring.From : ring.To;
    const Vector end     = forward ? ring.To : ring.From;
    const double size    = Length(seen.first);
    const double clear   = (((angle * Grow) + 0x1p-48) * size) + (2.0 * seen.second);
    const double enter   = Cross(start, seen.first);
    const double leave   = Cross(seen.first, end);
    if (!(std::abs(enter) > clear) || !(std::abs(leave) > clear))
    {
      return Side::Near;
    }
    const bool inside = ring.Whole
                     || (std::abs(ring.Sweep) <= std::numbers::pi ? enter > 0.0 && leave > 0.0
                                                                  : enter > 0.0 || leave > 0.0);
    return inside ? Side::In : Side::Out;
  }

  auto Bendings(const Ring &bent, const Ring &ring) -> Bend
  {
    const bool   hubbed = ring.Chord != nullptr;
    const Number gapX   = hubbed ? Sure::Exact(bent.Map.Origin.X) - Hubbed(ring, 0)
                                 : Sure::Apart(bent.Map.Origin.X, ring.Map.Origin.X);
    const Number gapY   = hubbed ? Sure::Exact(bent.Map.Origin.Y) - Hubbed(ring, 1)
                                 : Sure::Apart(bent.Map.Origin.Y, ring.Map.Origin.Y);
    const auto   row    = [](Vector way, const std::array<Number, 2> &column) -> Number
    {
      return (Sure::Exact(way.X) * column[0]) + (Sure::Exact(way.Y) * column[1]);
    };
    const std::array<Number, 2> axis  = {Sure::Exact(bent.Map.X.X), Sure::Exact(bent.Map.X.Y)};
    const std::array<Number, 2> other = {Sure::Exact(bent.Map.Y.X), Sure::Exact(bent.Map.Y.Y)};
    const Vector                first = ring.Rows[0];
    const Vector                next  = ring.Rows[1];
    return {.Bent   = &bent,
            .Turned = {row(first, axis), row(first, other), row(next, axis), row(next, other)},
            .Offset = {row(first, {gapX, gapY}), row(next, {gapX, gapY})},
            .Gap    = hubbed ? Vector{.X = gapX.Value.High, .Y = gapY.Value.High}
                             : bent.Map.Origin - ring.Map.Origin,
            .Drift  = hubbed ? (std::max(std::abs(gapX.Value.Low) + gapX.Bound,
                                         std::abs(gapY.Value.Low) + gapY.Bound)
                                * Grow)
                                   + Sure::Tiny
                             : 0.0};
  }
} // namespace Geoxy::Swift
