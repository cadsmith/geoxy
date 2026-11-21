/***************************************************************************************************
 *  @file      source/met/floor.hpp
 *  @brief     How near a point may come to a part, bounded below: the floor of each part, by which
 *             a point is told to lie farther from a curve than a length without Nearest.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_MET_FLOOR_HPP
#define GEOXY_MET_FLOOR_HPP

#include "met/parts.hpp"
#include "scratch.hpp"
#include "straight.hpp"

#include "geoxy/box.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <array>
#include <cstddef>
#include <span>
#include <variant>

namespace Geoxy::Met
{
  /**
   *  The rounding a bound on how near a point comes to a part allows for, past the width of a
   *  double: Remote holds a bound to be past a length by 2^-20 of itself and by 2^-36 of the size
   *  of the numbers about it, the part's and the point's.
   *
   *  A curve's point, as At works it out at a place, lies within some 2^8 steps of a double, at the
   *  size of the part's numbers and of the point, of the point of the part at that place worked out
   *  exactly: a run's start and way, an ellipse's origin and map, as a polyline's bulge or an arc
   *  is held, a piece's control points cut from a spline. Every place Nearest gives is such a
   *  point, and its distance to the point weighed is worked out within two steps more; so a point
   *  whose part lies farther than a length by both allowances, as the bound tells it, is nowhere
   *  within that length of the curve as Nearest finds it.
   */
  inline constexpr double Relative = 0x1p-20;
  inline constexpr double Absolute = 0x1p-36;

  /**
   *  The largest magnitude of a point's coordinates.
   *
   *  @param  point The point.
   *  @return       The magnitude.
   */
  auto Magnitude(Point point) -> double;

  /**
   *  How much a linear map stretches, at the least and at the most: its singular values, worked out
   *  with its entries scaled by one power of two to the size of 1, so that their squares neither
   *  run past every double nor come to nothing, and scaled back.
   *
   *  @pre        The map is invertible.
   *  @param  map The map, its origin aside.
   *  @return     The least and the most.
   */
  auto Stretches(const Matrix &map) -> std::array<double, 2>;

  /**
   *  An ellipse as a bound on how near a point comes to it is worked out: the inverse of its map,
   *  and how much the map stretches, at the least and at the most.
   */
  struct Oval
  {
    Matrix Unit  = {};  /**< The map's inverse. */
    double Least = 0.0; /**< Its least stretch. */
    double Most  = 0.0; /**< Its most.          */
  };

  /**
   *  A piece as a bound on how near a point comes to it is worked out: the box its control points
   *  span, and the strip they span across its chord, each holding it.
   */
  struct Shell
  {
    Box    Bounds = {};  /**< The box.                             */
    Point  From   = {};  /**< Where the chord starts.              */
    Vector Across = {};  /**< The chord's unit normal.             */
    double Low    = 0.0; /**< The least offset of a control point. */
    double High   = 0.0; /**< The greatest.                        */
  };

  /**
   *  A part as a bound on how near a point comes to it is worked out, once for every part: a run as
   *  it is, an ellipse by its unit circle's frame, a piece by its control points' box and strip;
   *  and the size of the numbers its points are worked out from.
   */
  struct Floor
  {
    std::variant<Run, Oval, Shell> Shape = Run{}; /**< The run, the ellipse or the piece. */
    double                         Size  = 0.0;   /**< The size of its numbers.           */
  };

  /**
   *  Whether a point lies certainly farther than a length from every part of a curve, so that
   *  Nearest would find none of it within the length: each part's bound past the length by the
   *  rounding it allows for.
   *
   *  @param  floors The curve's parts, as bounds are worked out from them.
   *  @param  point  The point.
   *  @param  length The length.
   *  @return        Whether it does; not where the point is not finite.
   */
  auto Remote(std::span<const Floor> floors, Point point, double length) -> bool;

  /**
   *  A curve's parts as bounds are worked out from them.
   *
   *  @param  parts The parts.
   *  @return       Their floors, in order.
   */
  auto Floors(const List<Part> &parts) -> List<Floor>;
} // namespace Geoxy::Met

#endif /* GEOXY_MET_FLOOR_HPP */
