/***************************************************************************************************
 *  @file      source/swift/roots.hpp
 *  @brief     The roots of a polynomial held by its Bernstein coefficients on a stretch from 0 to
 *             1, isolated in doubles under one running bound, each alone where the polynomial is
 *             monotone, and a place near each polished in doubles: what the swift kernels that meet
 *             a curve's piece by an equation on it share.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SWIFT_ROOTS_HPP
#define GEOXY_SWIFT_ROOTS_HPP

#include "real.hpp"
#include "sure.hpp"
#include "swift/casteljau.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>

namespace Geoxy::Swift
{
  /**
   *  The unit of a double's rounding, u = 2^-53.
   */
  inline constexpr double Unit = 0x1p-53;

  /**
   *  How many halvings a stretch may take before a kernel gives up on it: 2^-14 of the curve.
   */
  inline constexpr int Deepest = 14;

  /**
   *  How many steps of Newton's method, kept within the stretch, a root takes in doubles before its
   *  one step in double-double: enough to bring it within about 2^-30 of the stretch, as it
   *  converges quadratically once near; a root left farther off has a larger bound on that one
   *  step, and is only less often settled.
   */
  inline constexpr int Polishing = 6;

  /**
   *  A stretch of a piece's parameter and an equation's Bernstein coefficients on it, in doubles.
   *
   *  @tparam CountT How many coefficients it holds room for.
   */
  template <std::size_t CountT>
  struct Stretch
  {
    std::array<double, CountT> C     = {};  /**< The coefficients on it.    */
    double                     Lo    = 0.0; /**< Where it starts.           */
    double                     Hi    = 1.0; /**< Where it ends.             */
    int                        Depth = 0;   /**< How many halvings made it. */
  };

  /**
   *  A stretch holding one simple root, the equation monotone on it.
   */
  struct Root
  {
    double Lo    = 0.0; /**< Where it starts.                                */
    double Hi    = 1.0; /**< Where it ends.                                  */
    double Slope = 0.0; /**< The least |f'| on it, in the piece's parameter. */
  };

  /**
   *  The roots found, held in place: no more than the degree, one fewer than the coefficients.
   *
   *  @tparam CountT How many coefficients the equation holds room for.
   */
  template <std::size_t CountT>
  struct Roots
  {
    std::array<Root, CountT - 1> At    = {}; /**< The roots.          */
    std::size_t                  Count = 0;  /**< How many there are. */
  };

  /**
   *  A slope to divide by: itself, or the least normal double of its sign where it is zero, so that
   *  a step taken by it is huge, not infinite, and its bound settles nothing.
   *
   *  @param  slope The slope.
   *  @return       The slope to divide by.
   */
  inline auto Steady(double slope) -> double
  {
    return std::copysign(std::max(std::abs(slope), std::numeric_limits<double>::min()), slope);
  }

  /**
   *  What isolating an equation's roots carries from stretch to stretch: the coefficients' bound,
   *  the largest's size and the reach T, their degree, and the roots found.
   *
   *  @tparam CountT How many coefficients the equation holds room for.
   */
  template <std::size_t CountT>
  struct Isolation
  {
    std::array<double, 3> Bounds = {}; /**< The bound, the largest's size, and T. */
    std::size_t           Degree = 0;  /**< The degree.                           */
    Roots<CountT>         Found  = {}; /**< The roots found.                      */
  };

  /**
   *  A stretch of the isolation: far, a root's, or halved and each half isolated, the left first,
   *  so that the roots come in order along the piece.
   *
   *  @tparam CountT    How many coefficients the equation holds room for.
   *  @param  stretch   The stretch.
   *  @param  isolation What the isolation carries.
   *  @return           Whether every root on it was isolated within Deepest halvings.
   */
  template <std::size_t CountT>
  auto Isolating(const Stretch<CountT> &stretch, Isolation<CountT> &isolation) -> bool
  {
    const auto [first, largest, reach] = isolation.Bounds;
    const std::size_t degree           = isolation.Degree;
    const auto        order            = static_cast<double>(degree);
    const double      bound =
        (first + (order * static_cast<double>(stretch.Depth) * 2.0 * Unit * largest)) * Sure::Grow;
    const auto   coefficients = std::span{stretch.C}.first(degree + 1);
    const double head         = coefficients.front();
    const double tail         = coefficients.back();
    const bool   above        = std::ranges::all_of(coefficients, [&](double value) -> bool
                                                    { return value > reach + bound; });
    const bool   below        = std::ranges::all_of(coefficients, [&](double value) -> bool
                                                    { return value < -(reach + bound); });
    if (above || below)
    {
      return true;
    }
    /**
     *  The differences of the coefficients are n / (Hi - Lo) times f' 's own: each within twice the
     *  bound and its own rounding.
     */
    double      steepest = std::numeric_limits<double>::infinity();
    std::size_t rising   = 0;
    std::size_t falling  = 0;
    for (std::size_t index = 0; index < degree; ++index)
    {
      const double step   = coefficients[index + 1] - coefficients[index];
      const double slack  = (2.0 * bound) + (2.0 * Unit * std::abs(step));
      rising             += static_cast<std::size_t>(step > slack);
      falling            += static_cast<std::size_t>(step < -slack);
      steepest            = std::min(steepest, std::abs(step) - slack);
    }
    const bool monotone = rising == degree || falling == degree;
    const bool clear    = std::abs(head) > reach + bound && std::abs(tail) > reach + bound;
    /**
     *  Monotone with both ends clear, the ends lie on either side of zero: on one side, every
     *  coefficient between them would be clear as well, and the stretch far. Each such stretch
     *  holds a root of its own, so no more than the degree.
     */
    if (monotone && clear)
    {
      Roots<CountT> &roots       = isolation.Found;
      roots.At.at(roots.Count++) = {
          .Lo    = stretch.Lo,
          .Hi    = stretch.Hi,
          .Slope = steepest * order / (stretch.Hi - stretch.Lo) * Sure::Shrink};
      return true;
    }
    if (stretch.Depth >= Deepest)
    {
      return false;
    }
    const double    middle = stretch.Lo + ((stretch.Hi - stretch.Lo) / 2.0);
    Stretch<CountT> left   = {.C = {}, .Lo = stretch.Lo, .Hi = middle, .Depth = stretch.Depth + 1};
    Stretch<CountT> right  = {.C = {}, .Lo = middle, .Hi = stretch.Hi, .Depth = stretch.Depth + 1};
    std::array<double, CountT> work = stretch.C;
    left.C[0]                       = work[0];
    right.C.at(degree)              = work.at(degree);
    for (std::size_t round = 1; round <= degree; ++round)
    {
      for (std::size_t index = 0; index + round <= degree; ++index)
      {
        work.at(index) = (work.at(index) + work.at(index + 1)) * 0.5;
      }
      left.C.at(round)           = work[0];
      right.C.at(degree - round) = work.at(degree - round);
    }
    return Isolating(left, isolation) && Isolating(right, isolation);
  }

  /**
   *  The roots of an equation f on a piece, each alone in a stretch where f is monotone, isolated
   *  in doubles: every coefficient within one bound of the exact one, which each halving, a convex
   *  combination rounded once per round, grows by u times the largest coefficient per round. A
   *  stretch is far where every coefficient is clear of the reach T beyond the bound, so that f, in
   *  their hull, is too; it holds one simple root where f's differences are all of one sign, f
   *  monotone there, and its ends are clear of T on either side. Anything else is halved, and past
   *  Deepest halvings the kernel is not sure.
   *
   *  @tparam CountT How many coefficients the equation holds room for.
   *  @param  start  The coefficients on the whole piece, rounded.
   *  @param  bounds Their bound, the largest's size, and the reach T.
   *  @param  degree Their degree.
   *  @return        The roots; nothing where it is not sure.
   */
  template <std::size_t CountT>
  auto Isolated(const Stretch<CountT> &start, std::array<double, 3> bounds, std::size_t degree)
      -> std::optional<Roots<CountT>>
  {
    Isolation<CountT> isolation = {.Bounds = bounds, .Degree = degree, .Found = {}};
    if (!Isolating(start, isolation))
    {
      return std::nullopt;
    }
    return isolation.Found;
  }

  /**
   *  A root polished: a safeguarded Newton's method in doubles within its stretch, from where the
   *  chord of its ends crosses, to about 2^-30 of the stretch, a grid place there. A step that
   *  lands within the stretch on the place it left is where every step after would land, the same
   *  value at the same place moving the ends to it again, so the polish stops there: the same
   *  place, sooner.
   *
   *  @tparam CountT       How many coefficients the equation holds room for.
   *  @tparam StepsT       How many steps it takes: Polishing for a cubic's equation; more for one
   *                       of a higher degree, whose stretches are wider for as many halvings.
   *  @param  root         The root's stretch.
   *  @param  coefficients f's coefficients, rounded.
   *  @param  degree       The degree.
   *  @return              A place near the root, on the grid of 2^-53.
   */
  template <std::size_t CountT, int StepsT = Polishing>
  auto Guessed(const Root &root, const std::array<double, CountT> &coefficients, std::size_t degree)
      -> double
  {
    std::array<double, CountT> slopes = {};
    for (std::size_t index = 0; index < degree; ++index)
    {
      slopes.at(index) =
          static_cast<double>(degree) * (coefficients.at(index + 1) - coefficients.at(index));
    }
    double       low    = root.Lo;
    double       high   = root.Hi;
    double       atLow  = Valued(coefficients, degree, low);
    const double atHigh = Valued(coefficients, degree, high);
    double       along  = low + ((high - low) * (atLow / (atLow - atHigh)));
    for (int step = 0; step < StepsT; ++step)
    {
      const double value = Valued(coefficients, degree, along);
      if ((value < 0.0) == (atLow < 0.0))
      {
        low   = along;
        atLow = value;
      }
      else
      {
        high = along;
      }
      const double next = along - (value / Steady(Valued(slopes, degree - 1, along)));
      const bool   kept = next >= low && next <= high;
      if (kept && Geoxy::Exactly(next, along))
      {
        break;
      }
      along = kept ? next : low + ((high - low) / 2.0);
    }
    return std::round(along * 0x1p53) * 0x1p-53;
  }
} // namespace Geoxy::Swift

#endif /* GEOXY_SWIFT_ROOTS_HPP */
