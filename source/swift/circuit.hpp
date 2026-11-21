/***************************************************************************************************
 *  @file      source/swift/circuit.hpp
 *  @brief     A ring's equation along a piece, its roots, and each root's crossing certified.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_SWIFT_CIRCUIT_HPP
#define GEOXY_SWIFT_CIRCUIT_HPP

#include "scratch.hpp"
#include "sure.hpp"
#include "swift/casteljau.hpp"
#include "swift/lane.hpp"
#include "swift/piece.hpp"
#include "swift/rim.hpp"
#include "swift/roots.hpp"
#include "swift/section.hpp"
#include "swift/terms.hpp"
#include "trig.hpp"
#include "twofold.hpp"

#include "geoxy/point.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <utility>

namespace Geoxy::Swift
{
  /**
   *  The ring's equation on a piece, H(t) = |A (N(t) - o W(t))|^2 - rho^2 W(t)^2 for the piece's
   *  numerator N and weight W, held by its Bernstein coefficients of degree 2 n in doubles: c_k is
   *  the sum over i + j = k of C(n, i) C(n, j) (G_i . G_j - rho^2 w_i w_j), over C(2 n, k), G_i = A
   *  w_i (P_i - o). H is W^2 times the ring's quadratic Q at the piece's point, and a point at
   *  distance s of the carrier has |Q| at the least s (|A (p - o)| + rho) / kappa: where |H| is
   *  above T = reach W_max^2 (U + rho) / kappa, U the most |A (p - o)| may be on the piece, no more
   *  than on its control points, the piece lies farther than the reach from the carrier
   *  (tests/proof/ring.py's reach).
   *
   *  @tparam CountT The piece's room.
   */
  template <std::size_t CountT>
  struct Equation
  {
    Stretch<Ranks<CountT>> Whole   = {};  /**< The coefficients on the whole piece. */
    std::size_t            Degree  = 0;   /**< 2 n; 3 n for a cubic Bézier curve's. */
    double                 First   = 0.0; /**< How far a coefficient may be off.    */
    double                 Largest = 0.0; /**< The largest coefficient's size.      */
    double                 Reach   = 0.0; /**< T.                                   */
    double                 Heavy   = 0.0; /**< W_max^2; W_max^3 for a cubic's.      */
    double                 Turn    = 0.0; /**< The most |H''| may be.               */
  };

  /**
   *  The ring's equation on a piece. Each G_i is worked out in doubles from the moved control
   *  point: two products and a sum, within 4 u of the terms' sizes m_i, and past that by A's rows'
   *  sizes times what the moved point may be off. A product of two such is within |G_i| e_j + |G_j|
   *  e_i + e_i e_j of the exact one; each term, the sum over i + j = k by the binomial weights, and
   *  rho^2 times the weights' products taken off, then divided, are within 3 u, 3 u and 4 u of the
   *  sizes of their terms for a cubic, 16 u in all with room to spare (tests/proof/ring.py's
   *  product, term, summed and finished); a sum of n + 1 terms adds n u, and the weights' products
   *  and sum (n + 2) u, which 16 u holds for every degree n to Broadest (its running and wide).
   *  Every binomial weight and product of two is a whole number below 2^53, exact.
   *
   *  @tparam CountT The piece's room.
   *  @param  ring   The ring.
   *  @param  piece  The piece.
   *  @param  reach  The reach.
   *  @return        The equation.
   */
  template <std::size_t CountT>
  auto Equated(const Ring &ring, const Piece<CountT> &piece, double reach) -> Equation<CountT>
  {
    const std::size_t     degree = piece.Degree;
    const Vector          first  = ring.Rows[0];
    const Vector          other  = ring.Rows[1];
    Poled<Vector, CountT> terms  = {};
    Poled<Vector, CountT> offs   = {};
    double                most   = 0.0;
    double                heavy  = 0.0;
    for (std::size_t index = 0; index <= degree; ++index)
    {
      const Vector point  = piece.Moved.at(index);
      const Vector slack  = piece.Offs.at(index);
      const double weight = piece.Weights.at(index);
      const Vector term   = {.X = Dot(first, point), .Y = Dot(other, point)};
      terms.at(index)     = term;
      offs.at(index) =
          Vector{.X = (std::abs(first.X) * slack.X)
                    + (std::abs(first.Y) * slack.Y)
                    + (4.0 * Unit * (std::abs(first.X * point.X) + std::abs(first.Y * point.Y))),
                 .Y = (std::abs(other.X) * slack.X)
                    + (std::abs(other.Y) * slack.Y)
                    + (4.0 * Unit * (std::abs(other.X * point.X) + std::abs(other.Y * point.Y)))}
          * Grow;
      most  = std::max(most, Trig::Hypot(std::abs(term.X) + offs.at(index).X,
                                         std::abs(term.Y) + offs.at(index).Y)
                                 / (weight - piece.Heft));
      heavy = std::max(heavy, weight + piece.Heft);
    }
    const double     level    = ring.Level.Value.High;
    const double     drift    = std::abs(ring.Level.Value.Low) + ring.Level.Bound;
    Equation<CountT> equation = {.Degree = 2 * degree};
    for (std::size_t power = 0; power <= equation.Degree; ++power)
    {
      double sum    = 0.0;
      double size   = 0.0;
      double off    = 0.0;
      double weight = 0.0;
      for (std::size_t one = power > degree ? power - degree : 0; one <= std::min(power, degree);
           ++one)
      {
        const std::size_t two      = power - one;
        const double      times    = Choose(degree, one) * Choose(degree, two);
        const Vector      lead     = terms.at(one);
        const Vector      next     = terms.at(two);
        const Vector      leadOff  = offs.at(one);
        const Vector      nextOff  = offs.at(two);
        sum                       += times * Dot(lead, next);
        size   += times * (std::abs(lead.X * next.X) + std::abs(lead.Y * next.Y));
        off    += times
                * ((std::abs(lead.X) * nextOff.X)
                   + (std::abs(next.X) * leadOff.X)
                   + (leadOff.X * nextOff.X)
                   + (std::abs(lead.Y) * nextOff.Y)
                   + (std::abs(next.Y) * leadOff.Y)
                   + (leadOff.Y * nextOff.Y));
        weight += times * piece.Weights.at(one) * piece.Weights.at(two);
        /**
         *  Weights within the heft of the true ones: their product within |w_i| h + |w_j| h + h^2
         *  of the true product (tests/proof/ring.py's product), times rho^2.
         */
        off += times
             * std::abs(level)
             * (((piece.Weights.at(one) + piece.Weights.at(two)) * piece.Heft)
                + (piece.Heft * piece.Heft));
      }
      const double below         = Choose(equation.Degree, power);
      const double value         = (sum - (level * weight)) / below;
      equation.Whole.C.at(power) = value;
      equation.First             = std::max(
          equation.First,
          (((off + (16.0 * Unit * (size + (std::abs(level) * weight))) + (drift * weight)) / below)
           * Grow)
              + Sure::Tiny);
      equation.Largest = std::max(equation.Largest, std::abs(value));
    }
    const Terms<CountT> &values = equation.Whole.C;
    for (std::size_t power = 0; power + 2 <= equation.Degree; ++power)
    {
      equation.Turn =
          std::max(equation.Turn,
                   std::abs(values.at(power + 2) - (2.0 * values.at(power + 1)) + values.at(power))
                       + (8.0 * equation.First)
                       + (8.0 * Unit * equation.Largest));
    }
    const auto order  = static_cast<double>(equation.Degree);
    equation.Turn    *= order * (order - 1.0) * (1.0 + 0x1p-48);
    equation.Heavy    = heavy * heavy;
    equation.Reach    = reach * equation.Heavy * (most + ring.High) / ring.Least * Grow;
    return equation;
  }

  /**
   *  A ring and a piece being met, and what is worked out once about them.
   *
   *  @tparam CountT The piece's room.
   */
  template <std::size_t CountT>
  struct Scene
  {
    const Ring                    *Carrier   = nullptr; /**< The ring whose equation.  */
    const Piece<CountT>           *Part      = nullptr; /**< The piece it is met with. */
    Geoxy::Swift::Equation<CountT> Equation  = {};      /**< Its equation on it.       */
    Geoxy::Tolerance               Tolerance = {};      /**< The tolerance.            */
  };

  /**
   *  A crossing found: its point, the double nearest the true one, and its place on the piece's
   *  curve where the piece is a Bézier curve's, the double nearest the true one.
   */
  struct Found
  {
    Point  At    = {};    /**< The point.                                   */
    double Along = 0.0;   /**< The place on the piece's curve, where owned. */
    bool   Own   = false; /**< Whether the place is the curve's own.        */
  };

  /**
   *  H from a point on the unit circle's plane, the numerator's image A (N - o W), and the weight
   *  W: |A (N - o W)|^2 - rho^2 W^2, in double-double.
   *
   *  @param  ring   The ring.
   *  @param  seen   The image's coordinates.
   *  @param  weight The weight.
   *  @return        H.
   */
  auto Squaring(const Ring &ring, const std::array<Number, 2> &seen, const Number &weight)
      -> Number;

  /**
   *  H at a grid place of a rational span, in double-double: its true homogeneous point, as
   *  Homogeneously gives it, less o times its weight, through A.
   *
   *  @param  ring    The ring.
   *  @param  section The rational span.
   *  @param  along   The place, on the grid.
   *  @return         H there, and the point's numerator's coordinates and its weight.
   */
  auto Raising(const Ring &ring, const Broad &section, double along)
      -> std::pair<Number, std::array<Number, 3>>;

  /**
   *  H at a grid place of a quarter, in double-double: its direction at t through A, k W + B c for
   *  c the direction, the numerator of its point less o times the weight.
   *
   *  @param  ring  The ring.
   *  @param  piece The quarter.
   *  @param  along The place, on the grid.
   *  @return       H there, and no point: a quarter's is worked out from its own map.
   */
  auto Circuited(const Ring &ring, const Piece<Thin> &piece, double along)
      -> std::pair<Number, std::array<Number, 3>>;

  /**
   *  H at a grid place of a polynomial piece, in double-double: its true point from its own
   *  numbers, as Exactly gives it, less o, through A.
   *
   *  @tparam CountT  The piece's room.
   *  @param  ring    The ring.
   *  @param  section The piece.
   *  @param  along   The place, on the grid.
   *  @return         H there, and the point's coordinates and its weight, 1.
   */
  template <std::size_t CountT>
  auto Plained(const Ring &ring, const Sectional<CountT> &section, double along)
      -> std::pair<Number, std::array<Number, 3>>
  {
    const Number                weight = Sure::Exact(1.0);
    const std::array<Number, 2> exact  = Exactly(section, along);
    const Number                gapX   = exact[0] - Hubbed(ring, 0);
    const Number                gapY   = exact[1] - Hubbed(ring, 1);
    const auto                  row    = [&gapX, &gapY](Vector way) -> Number
    {
      return (Sure::Exact(way.X) * gapX) + (Sure::Exact(way.Y) * gapY);
    };
    return {Squaring(ring, {row(ring.Rows[0]), row(ring.Rows[1])}, weight),
            {exact[0], exact[1], weight}};
  }

  /**
   *  H at a grid place, in double-double: a polynomial piece's as Plained gives it, a rational
   *  span's as Raising gives it, a quarter's as Circuited gives it.
   *
   *  @tparam CountT The piece's room.
   *  @param  scene  The ring and the piece.
   *  @param  along  The place, on the grid.
   *  @return        H there, and a spline's or a Bézier curve's point, its numerator's coordinates
   *                 and its weight, unused for a quarter.
   */
  template <std::size_t CountT>
  auto Levelled(const Scene<CountT> &scene, double along)
      -> std::pair<Number, std::array<Number, 3>>
  {
    const Ring          &ring  = *scene.Carrier;
    const Piece<CountT> &piece = *scene.Part;
    if constexpr (CountT == Ample)
    {
      if (piece.Source->Rated)
      {
        return Raising(ring, *piece.Source, along);
      }
    }
    else
    {
      if (piece.Source == nullptr)
      {
        return Circuited(ring, piece, along);
      }
    }
    return Plained(ring, *piece.Source, along);
  }

  /**
   *  A rational span's point at a polished place, by Taylor from the grid place on each of its
   *  homogeneous coordinates, w x, w y and w, each a polynomial as the kernel of a line and a
   *  Bézier curve takes one: N(t*) = N(t0) + N'(t0) D + N''(xi) D^2 / 2, -step within off of D, N'
   *  worked out in doubles from the span's control points, within twice n times the slack and a few
   *  roundings of the most it may be; the point N / W. Where that does not settle, as about a heavy
   *  weight, the span's true point at the place itself, a number within off of t0 - step, by de
   *  Boor's algorithm.
   *
   *  @param  section The span.
   *  @param  point   Its homogeneous point at the grid place.
   *  @param  step    The grid place, the step from it and how far the root may be from the step.
   *  @return         The point's coordinates.
   */
  auto Taylored(const Broad &section, const std::array<Number, 3> &point,
                std::array<double, 3> step) -> std::array<Number, 2>;

  /**
   *  The point a polished place of a quarter stands for, from its own map at the place, a number
   *  within off of t0 - step: o' + M c / W.
   *
   *  @param  piece The quarter.
   *  @param  step  The grid place, the step from it and how far the root may be from the step.
   *  @return       The point's coordinates.
   */
  auto Veered(const Piece<Thin> &piece, std::array<double, 3> step) -> std::array<Number, 2>;

  /**
   *  The point a polished place stands for. A polynomial piece's by Taylor from the grid place, as
   *  the kernel of a line and a Bézier curve works it out: B(t*) = B(t0) + B'(t0) D + B''(xi) D^2 /
   *  2, -step within off of D; a rational span's as Taylored works it out; a quarter's as Veered
   *  does.
   *
   *  @tparam CountT The piece's room.
   *  @param  piece  The piece.
   *  @param  point  The polynomial piece's point at the grid place.
   *  @param  step   The grid place, the step from it and how far the root may be from the step.
   *  @return        The point's coordinates.
   */
  template <std::size_t CountT>
  auto Moved(const Piece<CountT> &piece, const std::array<Number, 3> &point,
             std::array<double, 3> step) -> std::array<Number, 2>
  {
    const auto [along, ahead, off] = step;
    if constexpr (CountT == Ample)
    {
      if (piece.Source->Rated)
      {
        return Taylored(*piece.Source, point, step);
      }
    }
    else
    {
      if (piece.Source == nullptr)
      {
        return Veered(piece, step);
      }
    }
    const Sectional<CountT>    &section = *piece.Source;
    const auto                  order   = static_cast<double>(section.Degree);
    const std::array<double, 2> paced   = Paced(section);
    Poled<double, CountT>       wayX    = {};
    Poled<double, CountT>       wayY    = {};
    for (std::size_t index = 0; index < section.Degree; ++index)
    {
      wayX.at(index) = order * (section.Poles.at(index + 1).X - section.Poles.at(index).X);
      wayY.at(index) = order * (section.Poles.at(index + 1).Y - section.Poles.at(index).Y);
    }
    const double reach = std::abs(ahead) + off;
    const double tail  = ((paced[0] * off)
                          + (paced[0] * 0x1p-48 * std::abs(ahead))
                          + (2.0 * order * section.Slack * std::abs(ahead))
                          + (paced[1] * reach * reach * 0.5))
                       * Grow;
    const auto   moved = [&](const Number &axis, const Poled<double, CountT> &way) -> Number
    {
      const Twofold::Two lean = Twofold::TwoProduct(Valued(way, section.Degree - 1, along), -ahead);
      return axis + Number{.Value = {.High = lean.High, .Low = lean.Low}, .Bound = tail};
    };
    return {moved(point[0], wayX), moved(point[1], wayY)};
  }

  /**
   *  A root's place, polished in doubles within its stretch to a grid place, and the equation's
   *  slope there in doubles, on as many coefficients as a room holds: the same doubles as on the
   *  whole room, a narrow room's copies short.
   *
   *  @tparam RoomT    The room, the degree and one at the least.
   *  @tparam CountT   The piece's room.
   *  @param  equation The equation.
   *  @param  root     The root's stretch.
   *  @return          The place, and the slope there to divide by.
   */
  template <std::size_t RoomT, std::size_t CountT>
  auto Narrowly(const Equation<CountT> &equation, const Root &root) -> std::array<double, 2>
  {
    const std::size_t         degree = equation.Degree;
    const auto                order  = static_cast<double>(degree);
    const auto               &whole  = equation.Whole.C;
    std::array<double, RoomT> values = {};
    std::array<double, RoomT> slopes = {};
    std::ranges::copy(std::span{whole}.first(degree + 1), values.begin());
    for (std::size_t index = 0; index < degree; ++index)
    {
      slopes.at(index) = order * (whole.at(index + 1) - whole.at(index));
    }
    const double along = Guessed<RoomT, 2 * Polishing>(root, values, degree);
    return {along, Steady(Valued(slopes, degree - 1, along))};
  }

  /**
   *  A root's place polished and the slope there, on a cubic's room where the equation's degree is
   *  no more than a cubic piece's, as Narrowly works them out.
   *
   *  @tparam CountT   The piece's room.
   *  @param  equation The equation.
   *  @param  root     The root's stretch.
   *  @return          The place, and the slope there to divide by.
   */
  template <std::size_t CountT>
  auto Burnished(const Equation<CountT> &equation, const Root &root) -> std::array<double, 2>
  {
    if constexpr (CountT == Ample)
    {
      return equation.Degree <= Narrow ? Narrowly<Narrow + 1>(equation, root)
                                       : Narrowly<Ranks<Ample>>(equation, root);
    }
    else
    {
      return Narrowly<Ranks<CountT>>(equation, root);
    }
  }

  /**
   *  How far H' in doubles from the rounded coefficients may be off: by 2 m E and (6 m^2 + 4 m) u
   *  of their largest at the most, grown (tests/proof/strand.py's slack, every degree to sixteen),
   *  within 2 m (E + 32 u largest) up to the tenth degree.
   *
   *  @param  degree  The equation's degree, m.
   *  @param  first   How far a coefficient may be off, E.
   *  @param  largest The largest coefficient's size.
   *  @return         The slack, past the least normal double.
   */
  auto Sloping(std::size_t degree, double first, double largest) -> double;

  /**
   *  A crossing's place on the piece's curve, as Spread gives it: a quarter's is no place of its
   *  curve, which the thorough way places by the point's angle.
   *
   *  @tparam CountT The piece's room.
   *  @param  piece  The piece.
   *  @param  place  The place on the piece, with its bound.
   *  @return        The curve's parameter, with its bound; nothing for a quarter.
   */
  template <std::size_t CountT>
  auto Lodged(const Piece<CountT> &piece, const Number &place) -> Number
  {
    if constexpr (CountT == Thin)
    {
      if (piece.Source == nullptr)
      {
        return Sure::Exact(0.0);
      }
    }
    return Spread(*piece.Source, place);
  }

  /**
   *  Whether a crossing of a quarter lies off its sweep, or near an end of it: its direction at t,
   *  (1 - t^2, 2 t) turned, its place within the bound and a step of the double-double's leading
   *  part, and the margin about it.
   *
   *  @param  piece The quarter.
   *  @param  place The crossing's place on it, with its bound.
   *  @param  apart How far the true root may lie from the place.
   *  @return       In, out, or near the sweep's end.
   */
  auto Swerving(const Piece<Thin> &piece, const Number &place, double apart) -> Side;

  /**
   *  One of the equation's roots: steep, polished, its point worked out; then kept where it lies on
   *  the piece's arc and on the ring's, each decided clear of where the curves may come within the
   *  reach about it, and its point, and a Bézier curve's place, settled.
   *
   *  The sine between the curves at the crossing is |H'| / (W^2 |grad Q| |p'|), |grad Q| = 2 |A^T A
   *  (p - o)| no more than 2 |A| rho there (tests/proof/ring.py's sine). The mean value theorem:
   *  the root t* lies within R = (|H(t0)| + E) / g of t0, g the least |H'| on the stretch, and the
   *  step t1 = t0 - H0 / H'(t0) is within [E + 2u |H0| + F (slack + |H''| R) / g] / |H'(t0)| of it,
   *  F the most |H(t0)| may be. Where the piece comes within the reach of the carrier, |H| is no
   *  more than T, so within T / g of a root on its stretch, and its point within speed T / g of
   *  that root's: a point of the ring that near lies within near = reach + speed T / g of the
   *  crossing, and within near / sigma_min of its direction on the unit circle; a quarter's place
   *  within T / g, its angle within 2 T / g. A direction that near has the crossing's side of each
   *  end of a sweep where the crossing's cross product with the end is clear of that margin
   *  (tests/proof/ring.py's margin).
   *
   *  @tparam CountT The piece's room.
   *  @param  scene  The ring and the piece.
   *  @param  root   The root's stretch.
   *  @param  found  Where the crossing goes, where it is one.
   *  @return        Whether it was decided.
   */
  template <std::size_t CountT>
  auto Crossed(const Scene<CountT> &scene, const Root &root, List<Found> &found) -> bool
  {
    const Ring             &ring     = *scene.Carrier;
    const Piece<CountT>    &piece    = *scene.Part;
    const Equation<CountT> &equation = scene.Equation;
    const double            reach    = Reach(scene.Tolerance);
    const double            sine =
        root.Slope / (equation.Heavy * 2.0 * ring.Most * ring.High * piece.Speed) * Shrink;
    if (!Steep(sine, scene.Tolerance))
    {
      return false;
    }
    const auto [along, slope] = Burnished(equation, root);
    const auto [value, point] = Levelled(scene, along);
    const double slack        = Sloping(equation.Degree, equation.First, equation.Largest);
    const double most         = (Sure::Size(value) + value.Bound) * Grow;
    const double far          = most / root.Slope * Grow;
    const double step         = (value.Value.High + value.Value.Low) / slope;
    const double off          = ((((value.Bound + (2.0 * Unit * Sure::Size(value)))
                                   + (most * (slack + (equation.Turn * far)) / root.Slope))
                                  / std::abs(slope))
                                 + (2.0 * Unit * std::abs(step)))
                              * Grow;
    const Number place = {.Value = Twofold::Add(Twofold::Lifted(along), Twofold::Lifted(-step)),
                          .Bound = off};
    const std::array<Number, 2> crossing = Moved(piece, point, {along, step, off});
    const double                apart    = equation.Reach / root.Slope * Grow;
    if constexpr (CountT == Thin)
    {
      const Side side = piece.Arc != nullptr ? Swerving(piece, place, apart) : Side::In;
      if (side != Side::In)
      {
        return side == Side::Out;
      }
    }
    const double near   = (reach + (piece.Speed * apart)) * Grow;
    const double moved  = Trig::Hypot(crossing[0].Bound + std::abs(crossing[0].Value.Low),
                                      crossing[1].Bound + std::abs(crossing[1].Value.Low));
    auto         seen   = Seen(ring, {.X = crossing[0].Value.High, .Y = crossing[1].Value.High});
    seen.first          = seen.first * ring.Facing;
    seen.second        += ring.Most * moved * Grow;
    const Side arc = ring.Chord != nullptr ? Sidled(*ring.Chord, crossing, near + (moved * Grow))
                                           : Swept(ring, seen, near / (ring.Least * ring.Low));
    if (arc != Side::In)
    {
      return arc == Side::Out;
    }
    /**
     *  A Bézier curve's place is settled with the point; a quarter's is not a place of its curve,
     *  which the thorough way places by the point's angle.
     */
    const std::optional<std::array<double, 3>> rounded =
        Sure::Rounded(std::array{crossing[0], crossing[1], Lodged(piece, place)});
    if (rounded)
    {
      found.push_back({.At    = {.X = (*rounded)[0], .Y = (*rounded)[1]},
                       .Along = (*rounded)[2],
                       .Own   = piece.Source != nullptr});
    }
    return rounded.has_value();
  }

  /**
   *  An equation's roots isolated on as many coefficients as a room holds, the same roots as on the
   *  whole room: a narrow room keeps the halvings' copies short.
   *
   *  @tparam RoomT    The room, the degree and one at the least.
   *  @tparam CountT   The piece's room.
   *  @param  equation The equation.
   *  @return          The roots; nothing where it is not sure.
   */
  template <std::size_t RoomT, std::size_t CountT>
  auto Rooted(const Equation<CountT> &equation) -> std::optional<Roots<Ranks<CountT>>>
  {
    Stretch<RoomT> whole = {};
    std::ranges::copy(std::span{equation.Whole.C}.first(equation.Degree + 1), whole.C.begin());
    const std::optional<Roots<RoomT>> roots =
        Isolated(whole, {equation.First, equation.Largest, equation.Reach}, equation.Degree);
    if (!roots)
    {
      return std::nullopt;
    }
    Roots<Ranks<CountT>> out = {.At = {}, .Count = roots->Count};
    std::ranges::copy(std::span{roots->At}.first(roots->Count), out.At.begin());
    return out;
  }

  /**
   *  An equation's roots isolated, on a cubic's room where its degree is no more than a cubic
   *  piece's, as Rooted isolates them.
   *
   *  @tparam CountT   The piece's room.
   *  @param  equation The equation.
   *  @return          The roots; nothing where it is not sure.
   */
  template <std::size_t CountT>
  auto Isolates(const Equation<CountT> &equation) -> std::optional<Roots<Ranks<CountT>>>
  {
    if constexpr (CountT == Ample)
    {
      return equation.Degree <= Narrow ? Rooted<Narrow + 1>(equation)
                                       : Rooted<Ranks<Ample>>(equation);
    }
    else
    {
      return Isolated(equation.Whole, {equation.First, equation.Largest, equation.Reach},
                      equation.Degree);
    }
  }

  /**
   *  A ring and a piece met: the equation's roots isolated, each crossing decided.
   *
   *  @tparam CountT The piece's room.
   *  @param  scene  The ring, the piece and the tolerance, the equation to be worked out.
   *  @param  found  Where the crossings go.
   *  @return        Whether every one was decided.
   */
  template <std::size_t CountT>
  auto Encountered(Scene<CountT> scene, List<Found> &found) -> bool
  {
    scene.Equation = Equated(*scene.Carrier, *scene.Part, Reach(scene.Tolerance));
    const std::optional<Roots<Ranks<CountT>>> roots = Isolates(scene.Equation);
    return roots
        && std::ranges::all_of(std::span{roots->At}.first(roots->Count),
                               [&scene, &found](const Root &root) -> bool
                               { return Crossed(scene, root, found); });
  }
} // namespace Geoxy::Swift

#endif /* GEOXY_SWIFT_CIRCUIT_HPP */
