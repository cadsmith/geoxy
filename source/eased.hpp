/***************************************************************************************************
 *  @file      source/eased.hpp
 *  @brief     A spline cut finer where its weights make it stiff, to be met in its place.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_EASED_HPP
#define GEOXY_EASED_HPP

#include "kernel/piece.hpp"
#include "twofold.hpp"

#include "geoxy/spline.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace Geoxy
{
  /**
   *  A place on a spline told from both ends of its span: how far past the span's start and how far
   *  short of its end, each in two doubles. Beside a heavy weight a place lies nearer a span's end
   *  than doubles are apart there, and its parameter alone, even in two doubles, tells that sliver
   *  to few of its digits; told from the end it lies by, it is told to all of them.
   */
  struct Sited
  {
    std::size_t     Span  = 0;  /**< The index of the knot that starts the span.        */
    Twofold::Double Ahead = {}; /**< How far past the span's start, in the parameter.   */
    Twofold::Double Short = {}; /**< How far short of the span's end, in the parameter. */
  };

  /**
   *  A piece of a spline's span, as cut for a stand-in: the span, its knots, and where the piece
   *  lies within it from 0 to 1, told from both ends.
   */
  struct Slice
  {
    std::size_t     Span  = 0;   /**< The index of the knot that starts the span.        */
    double          From  = 0.0; /**< The span's first knot.                             */
    double          To    = 0.0; /**< The span's last knot.                              */
    Twofold::Double Start = {};  /**< How far past the span's start the piece starts.    */
    Twofold::Double Rest  = {};  /**< How far short of the span's end the piece ends.    */
    double          Width = 1.0; /**< How much of the span the piece is: a power of two. */
  };

  /**
   *  A spline that stands for another, the same curve to a step of a double, cut finer where the
   *  other's weights make it stiff.
   *
   *  Beside a weight far heavier than its neighbours a spline rushes along its poles' legs within a
   *  sliver of its parameter, a sliver thinner than doubles are apart for a weight 2^52 times the
   *  rest: no double then tells a place on the leg. Each piece of the spline is halved, and its
   *  halves again, until the weights of each are within so many times the lighter of its ends, and
   *  each of those is a span of length 1 of the stand-in, whose parameter so tells every place of
   *  the curve as finely as its coordinates.
   */
  struct Easing
  {
    Spline             Shape  = {}; /**< The stand-in: a span of length 1 a piece.   */
    std::vector<Slice> Pieces = {}; /**< Where each of its spans lies on the spline. */
  };

  /**
   *  Within how many times the lighter of its end weights a piece's weights are, for the piece to
   *  stand whole where crossings are sought on it: each crossing found is then worked out again on
   *  the spline itself, in two doubles, so a piece need only be even enough to be solved on.
   */
  inline constexpr double Mild = 16.0;

  /**
   *  The same where a place found on a piece is told as it was found, its point worked out on the
   *  piece at the parameter found. A piece leaves its lighter end as many times faster as the
   *  weight beside that end is heavier, and a parameter near 1 is told to a step of a double and no
   *  finer: within twice, a step of the parameter moves the point a few steps of its coordinates.
   */
  inline constexpr double Gentle = 2.0;

  /**
   *  A spline's stand-in, where one is called for.
   *
   *  @param  spline The spline.
   *  @param  mild   Within how many times the lighter of its end weights each piece's weights are
   *                 to be: Mild or Gentle.
   *  @return        The stand-in; nothing for a spline that is not valid, has no weights, or has no
   *                 piece stiff enough to be cut. One of no pieces for a spline whose weights lie
   *                 farther apart than doubles hold beside each other, some 2^1070 times, and which
   *                 so cannot be stood for.
   */
  auto Eased(const Spline &spline, double mild) -> std::optional<Easing>;

  /**
   *  A spline as pieces even enough to be asked for places as they are found: each piece's weights
   *  within Gentle times the lighter of its ends, its derivative's numerator beside it, cut from
   *  the whole span's so that a piece far smaller than its coordinates still tells which way it
   *  runs.
   */
  struct Paving
  {
    std::vector<Kernel::Tread> Treads = {}; /**< The pieces, in the spline's order. */
    std::vector<Slice>         Pieces = {}; /**< Where each lies on the spline.     */
  };

  /**
   *  A spline's paving, where one is called for.
   *
   *  @param  spline The spline.
   *  @return        The paving; nothing for a spline that is not valid, has no weights, has no
   *                 piece stiff enough to be cut, or whose weights lie farther apart than doubles
   *                 hold beside each other, some 2^1070 times: that one is answered as its own
   *                 parameter tells it.
   */
  auto Paved(const Spline &spline) -> std::optional<Paving>;

  /**
   *  A place on a piece cut from a spline's span as the spline's parameter, the double nearest it:
   *  beside a weight 2^52 times its neighbours' or more, that may be the end of the span, for
   *  places apart on the curve.
   *
   *  @param  piece Where the piece lies on the spline.
   *  @param  along The piece's parameter there, from 0 to 1.
   *  @return       The spline's parameter there.
   */
  auto Homed(const Slice &piece, double along) -> double;

  /**
   *  A place of a stand-in as the place of the spline it stands for, told from both ends of its
   *  span.
   *
   *  @param  easing The stand-in.
   *  @param  along  The place on it.
   *  @return        The place on the spline.
   */
  auto Back(const Easing &easing, Twofold::Double along) -> Sited;

  /**
   *  A place of a stand-in as the place of the spline it stands for, the double nearest it: beside
   *  a weight 2^52 times its neighbours' or more, that may be the end of the span.
   *
   *  @param  easing The stand-in.
   *  @param  along  The place on it.
   *  @return        The spline's parameter there.
   */
  auto Homed(const Easing &easing, double along) -> double;
} // namespace Geoxy

#endif /* GEOXY_EASED_HPP */
