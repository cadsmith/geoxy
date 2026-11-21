/***************************************************************************************************
 *  @file      source/met/floor.cpp
 *  @brief     How near a point may come to a part, bounded below: the floor of each part, by which
 *             a point is told to lie farther from a curve than a length without Nearest.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "met/floor.hpp"

#include "basis.hpp"
#include "hull.hpp"
#include "kernel/piece.hpp"
#include "met/parts.hpp"
#include "scratch.hpp"
#include "straight.hpp"

#include "geoxy/box.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <variant>

namespace Geoxy::Met
{
  auto Magnitude(Point point) -> double
  {
    return std::max(std::abs(point.X), std::abs(point.Y));
  }

  auto Stretches(const Matrix &map) -> std::array<double, 2>
  {
    const double largest =
        std::max({std::abs(map.X.X), std::abs(map.X.Y), std::abs(map.Y.X), std::abs(map.Y.Y)});
    const int  power  = -std::ilogb(largest);
    const auto scaled = [power](Vector axis) -> Vector
    {
      return {.X = std::ldexp(axis.X, power), .Y = std::ldexp(axis.Y, power)};
    };
    const Vector x      = scaled(map.X);
    const Vector y      = scaled(map.Y);
    const double square = Dot(x, x) + Dot(y, y);
    const double area   = std::abs(Cross(x, y));
    const double gap    = std::sqrt(std::max(0.0, (square * square) - (4.0 * area * area)));
    const double most   = std::sqrt((square + gap) / 2.0);
    return {std::ldexp(area / most, -power), std::ldexp(most, -power)};
  }

  namespace
  {
    /**
     *  A part as a bound on how near a point comes to it is worked out from: a run as it is, none
     *  where it has no length; an ellipse by its map's inverse and stretches, none where the map
     *  has no inverse, as a flat ellipse's; a piece by the box and the strip across its chord its
     *  control points span, the box alone where its ends are one point.
     *
     *  @param  part The part.
     *  @return      Its floor.
     */
    auto Floored(const Part &part) -> Floor
    {
      if (const Run *run = std::get_if<Run>(&part.Shape))
      {
        const double size = Magnitude(run->Start) + Magnitude(Point{} + run->Way);
        if (!(Dot(run->Way, run->Way) > 0.0))
        {
          /**
           *  A run whose way a double cannot square, as a polyline's lone vertex, is bounded by
           *  nothing: an ellipse of no stretch.
           */
          return {.Shape = Oval{}, .Size = size};
        }
        return {.Shape = *run, .Size = size};
      }
      if (const auto *ellipse = std::get_if<Ellipse>(&part.Shape))
      {
        const Matrix &map = ellipse->Map;
        const double  size =
            Magnitude(map.Origin) + Magnitude(Point{} + map.X) + Magnitude(Point{} + map.Y);
        const std::optional<Matrix> unit = Inverse(map);
        if (!unit)
        {
          return {.Shape = Oval{}, .Size = size};
        }
        const auto [least, most] = Stretches(map);
        return {.Shape = Oval{.Unit = *unit, .Least = least, .Most = most}, .Size = size};
      }
      const auto  &piece  = std::get<Piece>(part.Shape);
      const Point  from   = Kernel::Start(piece);
      const Vector chord  = Kernel::End(piece) - from;
      const double reach  = Length(chord);
      const Vector across = reach > 0.0 ? Perpendicular(chord) / reach : Vector{};
      Shell        shell  = {.Bounds = Kernel::Hull(piece),
                             .From   = from,
                             .Across = across,
                             .Low    = std::numeric_limits<double>::infinity(),
                             .High   = -std::numeric_limits<double>::infinity()};
      for (const Homogeneous &pole : Kernel::Poles(piece))
      {
        const double offset = Dot(Projected(pole) - from, across);
        shell.Low           = std::min(shell.Low, offset);
        shell.High          = std::max(shell.High, offset);
      }
      const Box &bounds = shell.Bounds;
      return {.Shape = shell, .Size = std::max(Magnitude(bounds.Min), Magnitude(bounds.Max))};
    }

    /**
     *  How near a point comes to a run at the least: to its point nearest, worked out.
     *
     *  @param  run   The run.
     *  @param  point The point.
     *  @return       The distance.
     */
    auto Beneath(const Run &run, Point point) -> double
    {
      const double along = Dot(point - run.Start, run.Way) / Dot(run.Way, run.Way);
      return Distance(point,
                      run.Start + (run.Way * std::clamp(along, run.Along.Start, run.Along.End)));
    }

    /**
     *  How near a point comes to an ellipse at the least: the least stretch of its map times how
     *  far the point lies off the unit circle in its frame, as no point of the circle is nearer
     *  than that to the point's image there; less what the inverse's rounding may move the image
     *  by, as the map's numbers cancel in its determinant, some steps of a double times the most
     *  stretch and how far out the image is.
     *
     *  @param  oval  The ellipse.
     *  @param  point The point.
     *  @return       The distance; 0 or less for a flat ellipse.
     */
    auto Beneath(const Oval &oval, Point point) -> double
    {
      const double radius = Length((oval.Unit * point) - Point{});
      return (oval.Least * std::abs(radius - 1.0)) - (Absolute * oval.Most * (radius + 1.0));
    }

    /**
     *  How near a point comes to a piece at the least: to the box its control points span, or to
     *  the strip they span across its chord, whichever is farther, as each holds the piece.
     *
     *  @param  shell The piece.
     *  @param  point The point.
     *  @return       The distance.
     */
    auto Beneath(const Shell &shell, Point point) -> double
    {
      const double offset = Dot(point - shell.From, shell.Across);
      return std::max({Outlying(shell.Bounds, point), shell.Low - offset, offset - shell.High});
    }

  } // namespace

  auto Remote(std::span<const Floor> floors, Point point, double length) -> bool
  {
    const double size = Magnitude(point);
    return std::ranges::all_of(
        floors,
        [point, length, size](const Floor &floor) -> bool
        {
          const double bound = std::visit([point](const auto &shape) -> double
                                          { return Beneath(shape, point); }, floor.Shape);
          return bound * (1.0 - Relative) > length + (Absolute * (floor.Size + size));
        });
  }

  auto Floors(const List<Part> &parts) -> List<Floor>
  {
    List<Floor> floors;
    floors.reserve(parts.size());
    for (const Part &part : parts)
    {
      floors.push_back(Floored(part));
    }
    return floors;
  }
} // namespace Geoxy::Met
