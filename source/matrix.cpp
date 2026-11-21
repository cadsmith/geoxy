/***************************************************************************************************
 *  @file      source/matrix.cpp
 *  @brief     Affine maps: turns, mirrors, composition, inverses and stretch.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/matrix.hpp"

#include "finite.hpp"
#include "twofold.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/box.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>

namespace Geoxy
{
  namespace
  {
    auto Origin(Point point) -> Vector
    {
      return point - Point{};
    }

    auto Extent(const Matrix &map) -> double
    {
      return std::max(Largest(map.X), Largest(map.Y));
    }
  } // namespace

  auto Rotation(Angle angle) -> Matrix
  {
    const Vector along = Direction(angle);
    return {.X = along, .Y = Perpendicular(along)};
  }

  auto Rotation(Angle angle, Point about) -> Matrix
  {
    Matrix turn = Rotation(angle);
    turn.Origin = about - (turn * Origin(about));
    return turn;
  }

  auto Mirror(Point first, Point second) -> std::optional<Matrix>
  {
    const std::optional<Vector> along = Unit(second - first);
    if (!along)
    {
      return std::nullopt;
    }
    const double cosine = (along->X * along->X) - (along->Y * along->Y);
    const double sine   = 2.0 * along->X * along->Y;

    Matrix mirror = {.X = {.X = cosine, .Y = sine}, .Y = {.X = sine, .Y = -cosine}};
    mirror.Origin = first - (mirror * Origin(first));
    return mirror;
  }

  auto operator*(const Matrix &lhs, const Matrix &rhs) -> Matrix
  {
    return {.X = lhs * rhs.X, .Y = lhs * rhs.Y, .Origin = lhs * rhs.Origin};
  }

  auto operator*(const Matrix &lhs, Point rhs) -> Point
  {
    return lhs.Origin + ((lhs.X * rhs.X) + (lhs.Y * rhs.Y));
  }

  auto operator*(const Matrix &lhs, Vector rhs) -> Vector
  {
    return (lhs.X * rhs.X) + (lhs.Y * rhs.Y);
  }

  auto operator*(const Matrix &lhs, const Box &rhs) -> Box
  {
    if (IsEmpty(rhs))
    {
      return rhs;
    }
    const std::array corners{rhs.Min, Point{.X = rhs.Max.X, .Y = rhs.Min.Y}, rhs.Max,
                             Point{.X = rhs.Min.X, .Y = rhs.Max.Y}};

    Box mapped;
    for (const Point corner : corners)
    {
      mapped = Grown(mapped, lhs * corner);
    }
    return mapped;
  }

  auto Determinant(const Matrix &map) -> double
  {
    return Cross(map.X, map.Y);
  }

  auto Inverse(const Matrix &map) -> std::optional<Matrix>
  {
    const double extent = Extent(map);
    if (!(extent > 0.0) || !std::isfinite(extent))
    {
      return std::nullopt;
    }
    /**
     *  Brought to a size near 1 by a power of two, which rounds nothing.
     */
    const double scale = std::ldexp(1.0, std::ilogb(extent));
    const Vector x     = map.X / scale;
    const Vector y     = map.Y / scale;
    /**
     *  The determinant from its two products exactly, in two doubles a number: it is 0 only where
     *  the axes are parallel, and however nearly parallel they are it is right to its last bit.
     */
    const Twofold::Two along       = Twofold::TwoProduct(x.X, y.Y);
    const Twofold::Two across      = Twofold::TwoProduct(x.Y, y.X);
    const double       determinant = Twofold::Difference({.High = along.High, .Low = along.Low},
                                                         {.High = across.High, .Low = across.Low})
                                         .High;
    if (!(std::abs(determinant) > 0.0))
    {
      return std::nullopt;
    }
    const double factor = determinant * scale;

    Matrix inverse = {.X = Vector{.X = y.Y, .Y = -x.Y} / factor,
                      .Y = Vector{.X = -y.X, .Y = x.X} / factor};
    inverse.Origin = Point{} - (inverse * Origin(map.Origin));
    if (!Finite(inverse.X) || !Finite(inverse.Y) || !Finite(Origin(inverse.Origin)))
    {
      return std::nullopt;
    }
    return inverse;
  }

  auto Stretch(const Matrix &map) -> double
  {
    const double scale = Extent(map);
    if (!(scale > 0.0) || !std::isfinite(scale))
    {
      return scale;
    }
    const Vector x          = map.X / scale;
    const Vector y          = map.Y / scale;
    const double sum        = Length({.X = x.X + y.Y, .Y = y.X - x.Y});
    const double difference = Length({.X = x.X - y.Y, .Y = y.X + x.Y});
    return scale * ((sum * 0.5) + (difference * 0.5));
  }
} // namespace Geoxy
