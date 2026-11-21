/***************************************************************************************************
 *  @file      source/belt/carrier.hpp
 *  @brief     Where two pieces' lines or circles meet, in doubles; and the same lines and circles
 *             held in two doubles each, for the meetings and the distances a double is too coarse
 *             to decide.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_BELT_CARRIER_HPP
#define GEOXY_BELT_CARRIER_HPP

#include "belt/piece.hpp"
#include "twofold.hpp"

#include "geoxy/point.hpp"

#include <array>
#include <cstddef>

namespace Geoxy::Belt
{
  /**
   *  The points two lines or circles share.
   */
  struct Hits
  {
    std::array<Point, 2> At    = {};    /**< The points, the first Count. */
    std::size_t          Count = 0;     /**< How many: 0, 1 or 2.         */
    bool                 Close = false; /**< Whether they come so near touching that doubles
                                             cannot tell how many there are.               */
    Point Touch = {};                   /**< Where they come nearest touching, when close. */
    bool  Graze = false;                /**< Whether they come within the rounding of the numbers
                                             that place them of touching.                        */
  };

  /**
   *  The points the lines or circles of two pieces share.
   *
   *  @param  first  One piece.
   *  @param  second The other.
   *  @return        The points, wherever along the line or the circle they fall; none for two lines
   *                 within 1e-14 of one direction, or two circles about one center.
   */
  auto Met(const Piece &first, const Piece &second) -> Hits;

  /**
   *  A whole line, through a point along a way, or a whole circle, about a point, in two doubles a
   *  number.
   */
  struct Carrier
  {
    bool            Round = false; /**< Whether it is a circle.                  */
    Twofold::Double X     = {};    /**< The x of the line's point or the center. */
    Twofold::Double Y     = {};    /**< Its y.                                   */
    Twofold::Double WayX  = {};    /**< A line's way, its x: its end's less its
                                        start's, exactly.                        */
    Twofold::Double WayY   = {};   /**< Its y.             */
    Twofold::Double Radius = {};   /**< A circle's radius. */
  };

  /**
   *  A shape's piece as its carrier: a line's way, and an arc's center and radius worked out from
   *  its bulge, in two doubles.
   *
   *  @param  piece The piece.
   *  @param  bulge The bulge of the vertex it starts at.
   *  @return       The carrier.
   */
  auto Carried(const Piece &piece, double bulge) -> Carrier;

  /**
   *  A piece's carrier moved to one side.
   *
   *  @param  carrier The carrier.
   *  @param  left    Whether the piece, an arc, runs counterclockwise.
   *  @param  by      How far, to the piece's left when positive.
   *  @return         The carrier moved.
   */
  auto Shifted(const Carrier &carrier, bool left, double by) -> Carrier;

  /**
   *  The points two carriers share, decided in two doubles.
   *
   *  @param  first  One carrier.
   *  @param  second The other.
   *  @return        The points; none for two lines along one direction.
   */
  auto Met(const Carrier &first, const Carrier &second) -> Hits;

  /**
   *  A point of a carrier, and the distance it is measured against.
   */
  struct Probe
  {
    const Carrier *On   = nullptr; /**< The carrier. */
    Point          Near = {};      /**< A point near it: the carrier's point nearest
                                        this one is the point measured.              */
    double Distance = 0.0;         /**< The distance. */
  };

  /**
   *  The square of how far a carrier's point is from a shape's piece, less the distance squared, in
   *  two doubles. Which part of the piece is nearest is read in doubles: where that changes, the
   *  two readings agree.
   *
   *  @param  probe The point and the distance.
   *  @param  piece The shape's piece.
   *  @param  wide  The piece's carrier.
   *  @return       The difference: below 0 nearer than the distance.
   */
  auto Clear(const Probe &probe, const Piece &piece, const Carrier &wide) -> Twofold::Double;

  /**
   *  A point two carriers share, made over from the doubles that found it: stepped once onto both,
   *  the miss from each read in two doubles a number.
   *
   *  @param  first  One carrier.
   *  @param  second The other.
   *  @param  near   The point as doubles found it.
   *  @return        The doubles nearest the point; the point given where the two all but touch
   *                 there, and stepping leads nowhere.
   */
  auto Sharpened(const Carrier &first, const Carrier &second, Point near) -> Point;

  /**
   *  The square of how far a carrier's point is from a point, less the distance squared, in two
   *  doubles.
   *
   *  @param  probe The carrier's point and the distance.
   *  @param  point The point.
   *  @return       The difference: below 0 nearer than the distance.
   */
  auto Clear(const Probe &probe, Point point) -> Twofold::Double;

  /**
   *  How a carrier's point stands to another carrier, in two doubles a number.
   *
   *  @param  probe The point; its distance is not read.
   *  @param  edge  The other carrier.
   *  @return       For a line, the cross of its way with the way from its point to the point: above
   *                0 to the left of its way. For a circle, the square of how far the point is from
   *                its center less its radius squared: below 0 inside it.
   */
  auto Beside(const Probe &probe, const Carrier &edge) -> Twofold::Double;
} // namespace Geoxy::Belt

#endif /* GEOXY_BELT_CARRIER_HPP */
