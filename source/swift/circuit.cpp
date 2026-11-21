/***************************************************************************************************
 *  @file      source/swift/circuit.cpp
 *  @brief     A ring's equation along a piece, its roots, and each root's crossing certified.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "swift/circuit.hpp"

#include "sure.hpp"
#include "swift/casteljau.hpp"
#include "swift/lane.hpp"
#include "swift/piece.hpp"
#include "swift/rim.hpp"
#include "swift/roots.hpp"
#include "swift/section.hpp"
#include "swift/terms.hpp"
#include "twofold.hpp"

#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

namespace Geoxy::Swift
{
  namespace
  {
    /**
     *  A quarter's direction at a place, (1 - t^2, 2 t) turned by its quarters, its weight 1 + t^2
     *  first.
     *
     *  @param  quarter Which quarter.
     *  @param  along   t.
     *  @return         The weight and the direction's coordinates.
     */
    auto Facing(std::size_t quarter, const Number &along) -> std::array<Number, 3>
    {
      const Number square = along * along;
      const auto   way    = Quartered(
          std::array<Number, 2>{Sure::Exact(1.0) - square, Sure::Exact(2.0) * along}, quarter);
      return {Sure::Exact(1.0) + square, way[0], way[1]};
    }
  } // namespace

  auto Squaring(const Ring &ring, const std::array<Number, 2> &seen, const Number &weight) -> Number
  {
    return (seen[0] * seen[0]) + (seen[1] * seen[1]) - (ring.Level * weight * weight);
  }

  auto Raising(const Ring &ring, const Broad &section, double along)
      -> std::pair<Number, std::array<Number, 3>>
  {
    const std::array<Number, 3> lifted = Homogeneously(section, Sure::Exact(along));
    const Number                weight = lifted[2];
    const Number                gapX   = lifted[0] - (Hubbed(ring, 0) * weight);
    const Number                gapY   = lifted[1] - (Hubbed(ring, 1) * weight);
    const std::array<Number, 2> seen   = {
        (Sure::Exact(ring.Rows[0].X) * gapX) + (Sure::Exact(ring.Rows[0].Y) * gapY),
        (Sure::Exact(ring.Rows[1].X) * gapX) + (Sure::Exact(ring.Rows[1].Y) * gapY)};
    return {Squaring(ring, seen, weight), lifted};
  }

  auto Circuited(const Ring &ring, const Piece<Thin> &piece, double along)
      -> std::pair<Number, std::array<Number, 3>>
  {
    const Bend                 &bend   = *piece.Arc;
    const std::array<Number, 3> facing = Facing(piece.Quarter, Sure::Exact(along));
    const Number                weight = facing[0];
    std::array<Number, 2>       seen   = {};
    for (std::size_t axis = 0; axis < 2; ++axis)
    {
      seen.at(axis) = (bend.Offset.at(axis) * weight)
                    + (bend.Turned.at(2 * axis) * facing[1])
                    + (bend.Turned.at((2 * axis) + 1) * facing[2]);
    }
    return {Squaring(ring, seen, weight), {}};
  }

  auto Taylored(const Broad &section, const std::array<Number, 3> &point,
                std::array<double, 3> step) -> std::array<Number, 2>
  {
    const auto [along, ahead, off] = step;
    const std::size_t degree       = section.Degree;
    const auto        order        = static_cast<double>(degree);
    const Point       base         = section.Base;
    const auto        lifted       = [&section, base](std::size_t index) -> std::array<double, 3>
    {
      const Point  pole   = section.Poles.at(index);
      const double weight = section.Weights.at(index);
      return {pole.X + (base.X * weight), pole.Y + (base.Y * weight), weight};
    };
    /**
     *  The control points w P taken back from the base, w (P - B) + B w, each coordinate within two
     *  steps u of its size past the slack, and the slack of the weight times |B|: so each
     *  difference within four of them and the slack of its own and its weight's.
     */
    const double heavy   = 1.0 + std::abs(base.X) + std::abs(base.Y);
    double       largest = 0.0;
    for (std::size_t index = 0; index <= degree; ++index)
    {
      const auto cell = lifted(index);
      largest         = std::max({largest, std::abs(cell[0]), std::abs(cell[1])});
    }
    const double                        carried = 4.0 * Unit * largest;
    const double                        slack   = (section.Slack * heavy) + carried;
    std::array<Poled<double, Ample>, 3> ways    = {};
    double                              speed   = 0.0;
    double                              bend    = 0.0;
    for (std::size_t index = 0; index < degree; ++index)
    {
      double sum = 0.0;
      for (std::size_t axis = 0; axis < 3; ++axis)
      {
        const double change      = lifted(index + 1).at(axis) - lifted(index).at(axis);
        ways.at(axis).at(index)  = order * change;
        sum                     += std::abs(change);
      }
      speed = std::max(speed, sum);
    }
    for (std::size_t index = 0; index + 2 <= degree; ++index)
    {
      for (std::size_t axis = 0; axis < 3; ++axis)
      {
        bend = std::max(bend, std::abs((lifted(index + 2).at(axis) - lifted(index + 1).at(axis))
                                       - (lifted(index + 1).at(axis) - lifted(index).at(axis))));
      }
    }
    speed = (speed + (6.0 * slack * Grow)) * order * (1.0 + 0x1p-48);
    bend  = (bend + (4.0 * slack)) * order * (order - 1.0) * 2.0 * (1.0 + 0x1p-48);
    const double          reach = std::abs(ahead) + off;
    const double          tail  = ((speed * off)
                                   + (speed * 0x1p-48 * std::abs(ahead))
                                   + (2.0 * order * slack * std::abs(ahead))
                                   + (bend * reach * reach * 0.5))
                                * Grow;
    std::array<Number, 3> moved = {};
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      const Twofold::Two lean =
          Twofold::TwoProduct(Valued(ways.at(axis), degree - 1, along), -ahead);
      moved.at(axis) =
          point.at(axis) + Number{.Value = {.High = lean.High, .Low = lean.Low}, .Bound = tail};
    }
    const std::array<Number, 2> taylor = {Sure::Quotient(moved[0], moved[2]),
                                          Sure::Quotient(moved[1], moved[2])};
    if (Sure::Rounded(taylor))
    {
      return taylor;
    }
    const std::array<Number, 3> exact = Homogeneously(
        section,
        {.Value = Twofold::Add(Twofold::Lifted(along), Twofold::Lifted(-ahead)), .Bound = off});
    return {Sure::Quotient(exact[0], exact[2]), Sure::Quotient(exact[1], exact[2])};
  }

  auto Veered(const Piece<Thin> &piece, std::array<double, 3> step) -> std::array<Number, 2>
  {
    const auto [along, ahead, off]     = step;
    const Matrix               &map    = piece.Arc->Bent->Map;
    const std::array<Number, 3> facing = Facing(
        piece.Quarter,
        {.Value = Twofold::Add(Twofold::Lifted(along), Twofold::Lifted(-ahead)), .Bound = off});
    const auto coordinate = [&facing](double start, std::array<double, 2> axes) -> Number
    {
      return Sure::Exact(start)
           + Sure::Quotient((Sure::Exact(axes[0]) * facing[1]) + (Sure::Exact(axes[1]) * facing[2]),
                            facing[0]);
    };
    return {coordinate(map.Origin.X, {map.X.X, map.Y.X}),
            coordinate(map.Origin.Y, {map.X.Y, map.Y.Y})};
  }

  auto Sloping(std::size_t degree, double first, double largest) -> double
  {
    const auto order = static_cast<double>(degree);
    return (order <= 10.0 ? 2.0 * order * (first + (32.0 * Unit * largest))
                          : ((2.0 * order * first)
                             + (((6.0 * order * order) + (4.0 * order)) * Unit * largest))
                                * Grow)
         + std::numeric_limits<double>::min();
  }

  auto Swerving(const Piece<Thin> &piece, const Number &place, double apart) -> Side
  {
    const double t      = place.Value.High;
    const double wobble = (place.Bound + std::abs(place.Value.Low)) * Grow;
    const auto   way    = Quartered(std::array{(1.0 - t) * (1.0 + t), 2.0 * t}, piece.Quarter);
    return Swept(*piece.Arc->Bent, {Vector{.X = way[0], .Y = way[1]}, 0.0}, 2.0 * (apart + wobble));
  }
} // namespace Geoxy::Swift
