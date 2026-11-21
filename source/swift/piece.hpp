/***************************************************************************************************
 *  @file      source/swift/piece.hpp
 *  @brief     A rational Bézier piece in a ring's frame: its poles shifted, scaled and weighed.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SWIFT_PIECE_HPP
#define GEOXY_SWIFT_PIECE_HPP

#include "sure.hpp"
#include "swift/lane.hpp"
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
#include <numbers>
#include <optional>
#include <span>

namespace Geoxy::Swift
{
  /**
   *  A piece of the other curve as the ring is met with it: a rational Bézier piece of degree two
   *  to Broadest piece, its homogeneous control points w_i (P_i - o) about the ring's center
   *  rounded, how far each coordinate may be from the true one's, and its weights, exact and above
   *  zero; a polynomial piece's are all 1. A Bézier curve's piece and a quarter are held in a
   *  cubic's room, a spline's span in a broad one.
   *
   *  @tparam CountT The room: the degree and one at the most.
   */
  template <std::size_t CountT>
  struct Piece
  {
    Poled<Vector, CountT>    Moved   = {};      /**< w_i (P_i - o), rounded.  */
    Poled<Vector, CountT>    Offs    = {};      /**< How far each may be off. */
    Poled<double, CountT>    Weights = {};      /**< w_i.                     */
    std::size_t              Degree  = 0;       /**< Two to Broadest.         */
    double                   Speed   = 0.0;     /**< |p'|, at the most.       */
    const Sectional<CountT> *Source  = nullptr; /**< The polynomial piece it is;
                                                     none for a quarter.             */
    const Bend *Arc     = nullptr;              /**< What it is a quarter of. */
    std::size_t Quarter = 0;                    /**< Which quarter.           */
    double      Heft    = 0.0;                  /**< How far a weight may be
                                                     off: a rational span's slack.   */
  };

  /**
   *  A polynomial piece's speed and bend, at the most: n times its longest control step, and n (n
   *  - 1) times its largest second difference, with what their rounding adds and the slack, which
   *    moves a step's coordinates by twice itself and a second difference's by four times.
   *
   *  @tparam CountT  The piece's room.
   *  @param  section The piece.
   *  @return         The most |B'| and |B''| may be.
   */
  template <std::size_t CountT>
  auto Paced(const Sectional<CountT> &section) -> std::array<double, 2>
  {
    const std::span<const Point> poles = std::span{section.Poles}.first(section.Degree + 1);
    const auto                   order = static_cast<double>(section.Degree);
    double                       speed = 0.0;
    double                       bend  = 0.0;
    for (std::size_t index = 0; index < section.Degree; ++index)
    {
      speed = std::max(speed, Length(poles[index + 1] - poles[index]));
    }
    for (std::size_t index = 0; index + 2 <= section.Degree; ++index)
    {
      const Vector second =
          (poles[index + 2] - poles[index + 1]) - (poles[index + 1] - poles[index]);
      bend = std::max({bend, std::abs(second.X), std::abs(second.Y)});
    }
    speed += 2.0 * std::numbers::sqrt2 * section.Slack * Grow;
    bend  += 4.0 * section.Slack;
    return {speed * order * (1.0 + 0x1p-48), bend * order * (order - 1.0) * 2.0 * (1.0 + 0x1p-48)};
  }

  /**
   *  A piece moved to a bulge's circle's true center C, in double-double: each control point w P
   *  less (C - B) w about the section's base B, rounded, within the double-double's low part and
   *  bound; a rational span's weight within the slack h of the true one moves it by h |C - B|, a
   *  polynomial piece's heft nothing. The rounded center's drift, a step of the doubles at the
   *  center's size, is never charged: met far from nothing, a bulge and a curve a few tolerances
   *  across are held to their own size.
   *
   *  @tparam CountT  The piece's room.
   *  @param  section The piece.
   *  @param  ring    The bulge's circle.
   *  @param  piece   The piece moved, its control points, slacks and weights filled in.
   */
  template <std::size_t CountT>
  void Centered(const Sectional<CountT> &section, const Ring &ring, Piece<CountT> &piece)
  {
    const Number hubX  = Hubbed(ring, 0) - Sure::Exact(section.Base.X);
    const Number hubY  = Hubbed(ring, 1) - Sure::Exact(section.Base.Y);
    const Vector rated = {.X = piece.Heft * Sure::Most(hubX), .Y = piece.Heft * Sure::Most(hubY)};
    for (std::size_t index = 0; index <= section.Degree; ++index)
    {
      const Point  pole       = section.Poles.at(index);
      const double weight     = section.Weights.at(index);
      const Number movedX     = Sure::Exact(pole.X) - (hubX * Sure::Exact(weight));
      const Number movedY     = Sure::Exact(pole.Y) - (hubY * Sure::Exact(weight));
      const Vector moved      = {.X = movedX.Value.High, .Y = movedY.Value.High};
      piece.Moved.at(index)   = moved;
      piece.Offs.at(index)    = Vector{.X = (Unit * std::abs(moved.X))
                                          + std::abs(movedX.Value.Low)
                                          + movedX.Bound
                                          + section.Slack
                                          + rated.X,
                                       .Y = (Unit * std::abs(moved.Y))
                                          + std::abs(movedY.Value.Low)
                                          + movedY.Bound
                                          + section.Slack
                                          + rated.Y}
                              * Grow;
      piece.Weights.at(index) = weight;
    }
  }

  /**
   *  A Bézier curve's piece moved to a point, a ring's center: each control point's gap from it,
   *  rounded once, and the slack.
   *
   *  @param section The piece.
   *  @param origin  The point.
   *  @param piece   The piece moved, its control points, slacks and weights filled in.
   */
  void Shifting(const Section &section, Point origin, Piece<Thin> &piece);

  /**
   *  A spline's span moved to a point, a ring's center or a parabola's middle control point, about
   *  its base: each control point's gap from it, rounded once, and the slack.
   *
   *  @param section The span.
   *  @param center  The point.
   *  @param piece   The piece moved, its control points, slacks and weights filled in.
   */
  void Shifting(const Broad &section, Point center, Piece<Ample> &piece);

  /**
   *  A polynomial or rational piece readied to be moved: its degree, speed and a rational span's
   *  heft, its control points and weights not yet filled in.
   *
   *  @tparam CountT  The piece's room.
   *  @param  section The piece.
   *  @return         The piece; nothing where a control point lies outside the range the kernels
   *                  take.
   */
  template <std::size_t CountT>
  auto Ready(const Sectional<CountT> &section) -> std::optional<Piece<CountT>>
  {
    if (!std::ranges::all_of(std::span{section.Poles}.first(section.Degree + 1),
                             [](Point pole) -> bool { return Held({pole.X, pole.Y}); }))
    {
      return std::nullopt;
    }
    Piece<CountT> piece = {.Degree = section.Degree, .Source = &section};
    if constexpr (CountT == Ample)
    {
      piece.Speed = section.Rated ? Fastest(section) : Paced(section)[0];
      piece.Heft  = section.Rated ? section.Slack : 0.0;
    }
    else
    {
      piece.Speed = Paced(section)[0];
    }
    return piece;
  }

  /**
   *  A polynomial or rational piece moved to a ring's center, as Shifting moves it; about a bulge's
   *  circle's true center, as Centered moves it.
   *
   *  @tparam CountT  The piece's room.
   *  @param  section The piece.
   *  @param  ring    The ring.
   *  @return         The piece; nothing where a control point lies outside the range the kernels
   *                  take.
   */
  template <std::size_t CountT>
  auto Plain(const Sectional<CountT> &section, const Ring &ring) -> std::optional<Piece<CountT>>
  {
    std::optional<Piece<CountT>> piece = Ready(section);
    if (!piece)
    {
      return std::nullopt;
    }
    if (ring.Chord != nullptr)
    {
      Centered(section, ring, *piece);
      return piece;
    }
    Shifting(section, ring.Map.Origin, *piece);
    return piece;
  }

  /**
   *  A quarter of an arc or an ellipse as a rational quadratic piece about the ring's center: the
   *  unit circle from angle q pi / 2 to (q + 1) pi / 2 is (1 - t^2, 2 t) / (1 + t^2), whose control
   *  points are (1, 0), (1, 1) and (0, 1) turned by q quarters, of weights 1, 1 and 2, every one
   *  exact (tests/proof/ring.py's quarter); mapped, each is o' + M c, about o w (o' - o + M c),
   *  worked out in doubles: the gap rounded and two sums, each within u of the sizes of its terms,
   *  and a bulge's circle's gap off by its drift past that, w times it. At t its angle is q pi / 2
   *  + 2 atan t, which moves no faster than 2, so the piece's speed is no more than twice
   *  sigma_max(M).
   *
   *  @param  bend    What it is a quarter of, about the ring.
   *  @param  quarter Which quarter, from 0 to 3.
   *  @return         The piece.
   */
  auto QuarterPiece(const Bend &bend, std::size_t quarter) -> Piece<Thin>;
} // namespace Geoxy::Swift

#endif /* GEOXY_SWIFT_PIECE_HPP */
