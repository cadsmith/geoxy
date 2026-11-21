/***************************************************************************************************
 *  @file      source/swift/piece.cpp
 *  @brief     A rational Bézier piece in a ring's frame: its poles shifted, scaled and weighed.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "swift/piece.hpp"

#include "real.hpp"
#include "swift/rim.hpp"
#include "swift/roots.hpp"
#include "swift/section.hpp"
#include "swift/terms.hpp"

#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <array>
#include <cmath>
#include <cstddef>

namespace Geoxy::Swift
{
  void Shifting(const Section &section, Point origin, Piece<Thin> &piece)
  {
    for (std::size_t index = 0; index <= section.Degree; ++index)
    {
      const Point  pole   = section.Poles.at(index);
      const double weight = section.Weights.at(index);
      const Vector moved  = {.X = pole.X - (origin.X * weight), .Y = pole.Y - (origin.Y * weight)};
      piece.Moved.at(index)   = moved;
      piece.Offs.at(index)    = Vector{.X = (Unit * std::abs(moved.X)) + section.Slack,
                                       .Y = (Unit * std::abs(moved.Y)) + section.Slack}
                              * Grow;
      piece.Weights.at(index) = weight;
    }
  }

  void Shifting(const Broad &section, Point center, Piece<Ample> &piece)
  {
    const bool  based = !Zero(section.Base.X) || !Zero(section.Base.Y);
    const Point origin =
        based ? Point{.X = center.X - section.Base.X, .Y = center.Y - section.Base.Y} : center;
    for (std::size_t index = 0; index <= section.Degree; ++index)
    {
      const Point  pole   = section.Poles.at(index);
      const double weight = section.Weights.at(index);
      /**
       *  A control point is w (P - B) about the section's base B, the ring's center less it, o - B,
       *  rounded once where B is not nothing, within u |o - B| (w + h), and a rational span's share
       *  (o - B) w rounded once, both off by the slack h: (1 + |o - B|) times it a coordinate.
       */
      const Vector moved = {.X = pole.X - (origin.X * weight), .Y = pole.Y - (origin.Y * weight)};
      const Vector rated = section.Rated ? Vector{.X = (Unit * std::abs(origin.X * weight))
                                                     + (piece.Heft * std::abs(origin.X)),
                                                  .Y = (Unit * std::abs(origin.Y * weight))
                                                     + (piece.Heft * std::abs(origin.Y))}
                                         : Vector{};
      const double carried =
          based ? Unit * (std::abs(origin.X) + std::abs(origin.Y)) * (weight + piece.Heft) : 0.0;
      piece.Moved.at(index) = moved;
      piece.Offs.at(index) =
          Vector{.X = (Unit * std::abs(moved.X)) + section.Slack + rated.X + carried,
                 .Y = (Unit * std::abs(moved.Y)) + section.Slack + rated.Y + carried}
          * Grow;
      piece.Weights.at(index) = weight;
    }
  }

  auto QuarterPiece(const Bend &bend, std::size_t quarter) -> Piece<Thin>
  {
    constexpr std::array<std::array<double, 2>, 3> Corners = {{{1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}}};
    constexpr std::array<double, 3>                Weights = {1.0, 1.0, 2.0};
    const Matrix                                  &map     = bend.Bent->Map;
    const Vector                                   gap     = bend.Gap;
    Piece<Thin>                                    piece   = {
        .Degree = 2, .Speed = 2.0 * bend.Bent->Most, .Arc = &bend, .Quarter = quarter};
    for (std::size_t index = 0; index < Corners.size(); ++index)
    {
      const auto   corner   = Quartered(Corners.at(index), quarter);
      const double weight   = Weights.at(index);
      const Vector moved    = gap + (map.X * corner[0]) + (map.Y * corner[1]);
      piece.Moved.at(index) = moved * weight;
      piece.Offs.at(index)  = (Vector{.X = std::abs(gap.X) + std::abs(map.X.X) + std::abs(map.Y.X),
                                      .Y = std::abs(gap.Y) + std::abs(map.X.Y) + std::abs(map.Y.Y)}
                               * (4.0 * Unit * weight * Grow))
                            + (Vector{.X = bend.Drift, .Y = bend.Drift} * (weight * Grow));
      piece.Weights.at(index) = weight;
    }
    return piece;
  }
} // namespace Geoxy::Swift
