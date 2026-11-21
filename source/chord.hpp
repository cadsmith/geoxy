/***************************************************************************************************
 *  @file      source/chord.hpp
 *  @brief     The two ends M + H and M - H of a common chord of a line's or a circle's carrier and
 *             a circle, each coordinate one exact sum rounded once and certified.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_CHORD_HPP
#define GEOXY_CHORD_HPP

#include "lean.hpp"
#include "twofold.hpp"

#include "geoxy/point.hpp"

#include <array>
#include <cmath>
#include <cstddef>

namespace Geoxy::Chord
{
  /**
   *  What the two crossings M + H and M - H of a chord share: M, each coordinate as a double, what
   *  it leaves and its terms' size, as Moved gives it; H in double-double; and how far each
   *  coordinate may be off from the error of the places M and H are worked from.
   */
  struct Halves
  {
    std::array<Lean::Rounded, 2> Middle = {}; /**< M, each coordinate.                     */
    std::array<Lean::Double, 2>  Half   = {}; /**< H, each coordinate.                     */
    std::array<double, 2>        Off    = {}; /**< How far each may be off for the places. */
  };

  /**
   *  A crossing worked out, and whether its point is certain.
   */
  struct Settle
  {
    Point               At     = {};    /**< The point, rounded.                     */
    bool                Sure   = false; /**< Whether it is the double nearest it.    */
    std::array<bool, 2> Held   = {};    /**< Whether each coordinate is the nearest. */
    std::array<bool, 2> Naught = {};    /**< Whether each, not held, may be zero: within
                                             its bound of it.                             */
  };

  /**
   *  The two ends M + H and M - H of a chord, each coordinate one exact sum and the certificate.
   *
   *  M's coordinate is within 50 u^2 of its terms' size besides what its place moves it, as Moved
   *  leaves it, and H's within 17.3 u^2 |H|, as a product in double-double leaves it; the sum is
   *  exact, its low part rounds twice, within 4.1 u^2 (|M| + |H|), and the double and its rest
   *  twice more, within 6 u^2 (|M| + |H|). In all within a quarter of Grain (|M| + |H|), with the
   *  places' share, Off, at most half what it is given here: under 0.6 of the bound.
   *
   *  @param  halves M, H, and how far the places move them.
   *  @return        M + H and M - H.
   */
  inline auto Ends(const Halves &halves) -> std::array<Settle, 2>
  {
    std::array<Settle, 2> ends = {};
    for (std::size_t side = 0; side < 2; ++side)
    {
      const double          sign   = side == 0 ? 1.0 : -1.0;
      std::array<double, 2> spot   = {};
      std::array<bool, 2>   held   = {};
      std::array<bool, 2>   naught = {};
      for (std::size_t axis = 0; axis < 2; ++axis)
      {
        const Lean::Rounded &middle = halves.Middle.at(axis);
        const Lean::Double  &half   = halves.Half.at(axis);
        const Twofold::Two   sum    = Twofold::TwoSum(middle.Value, sign * half.High);
        const double         low    = sum.Low + (middle.Rest + (sign * half.Low));
        const double         value  = sum.High + low;
        const double         bound =
            halves.Off.at(axis) + (Lean::Grain * (middle.Size + std::abs(half.High)));
        const bool certain = Lean::Certain(value, (sum.High - value) + low, bound);
        held.at(axis)      = certain;
        naught.at(axis)    = !certain && !(std::abs(value) > bound);
        spot.at(axis)      = value;
      }
      ends.at(side) = {.At     = {.X = spot[0], .Y = spot[1]},
                       .Sure   = static_cast<int>(held[0]) + static_cast<int>(held[1]) == 2,
                       .Held   = held,
                       .Naught = naught};
    }
    return ends;
  }

  /**
   *  A chord's two ends, a coordinate of one not certain but zero exactly made so: there the two
   *  halves of the sum, M's and H's of opposite signs, cancel exactly, which Intersect's quick
   *  roots never settle and its short sums and exact way round to +0, its numerator's square equal
   *  to the root's share's; the end then sure. Only the end whose bound reaches zero is asked,
   *  where the other's does not, as M not zero puts one root only on the axis, and only where its
   *  other coordinate is certain: the point of zero and that coordinate is the crossing where it
   *  lies on both carriers exactly, as the test given tells it.
   *
   *  @param  ends   M + H and M - H, as Ends works them out.
   *  @param  halves M and H.
   *  @param  zero   Whether a point lies on both carriers exactly, asked with the point and M's
   *                 coordinate where it is zero.
   *  @return        The ends.
   */
  template <typename ZeroT>
  auto Zeroed(std::array<Settle, 2> ends, const Halves &halves, const ZeroT &zero)
      -> std::array<Settle, 2>
  {
    if (ends[0].Sure && ends[1].Sure)
    {
      return ends;
    }
    for (std::size_t side = 0; side < 2; ++side)
    {
      Settle               &end   = ends.at(side);
      const Settle         &other = ends.at(1 - side);
      std::array<double, 2> spot  = {end.At.X, end.At.Y};
      for (std::size_t axis = 0; axis < 2; ++axis)
      {
        const double          middle = halves.Middle.at(axis).Value;
        const double          half   = (side == 0 ? 1.0 : -1.0) * halves.Half.at(axis).High;
        const int             asked  = static_cast<int>(end.Naught.at(axis))
                                     * static_cast<int>(!other.Naught.at(axis))
                                     * static_cast<int>(end.Held.at(1 - axis))
                                     * static_cast<int>(middle * half < 0.0);
        std::array<double, 2> naught = spot;
        naught.at(axis)              = 0.0;
        if (asked > 0 && zero(Point{.X = naught[0], .Y = naught[1]}, middle))
        {
          spot.at(axis)     = 0.0;
          end.Held.at(axis) = true;
        }
      }
      end.At   = {.X = spot[0], .Y = spot[1]};
      end.Sure = static_cast<int>(end.Held[0]) + static_cast<int>(end.Held[1]) == 2;
    }
    return ends;
  }
} // namespace Geoxy::Chord

#endif /* GEOXY_CHORD_HPP */
