/***************************************************************************************************
 *  @file      source/near/plumb.cpp
 *  @brief     The place of a spline nearest a point: its spans arched, sounded and pruned, and the
 *             feet on what is left weighed.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "near/plumb.hpp"

#include "basis.hpp"
#include "eased.hpp"
#include "hull.hpp"
#include "kernel/piece.hpp"
#include "real.hpp"
#include "scratch.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/at.hpp"
#include "geoxy/box.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <optional>
#include <span>
#include <vector>

namespace Geoxy::Near
{
  namespace
  {
    /**
     *  The farthest from zero a coordinate is weighed by branch and bound: so near that no distance
     *  between two points runs past the doubles.
     */
    constexpr double Widest = 0x1p500;

    /**
     *  How many times a span's piece is halved at most to bound how near it comes: to an eighth of
     *  it, all but straight.
     */
    constexpr int Finest = 3;

    /**
     *  How many control points a span weighed by branch and bound has at most: a spline of a higher
     *  degree has every span solved.
     */
    constexpr std::size_t MostPoles = 4;

    /**
     *  Control points of a piece weighed by branch and bound, the first degree and one.
     */
    using Arm = std::array<Point, MostPoles>;

    /**
     *  A span of a polynomial spline as a Bézier piece: its control points, each the span's blossom
     *  at its ends by de Boor's rounds in doubles, the rounds at the start shared by every control
     *  point that takes them first.
     *
     *  Each round is the one Section takes, cell for cell, so for a spline of no weights, whose
     *  every weight Section sums to exactly 1, the points are Section's to the bit, and the first
     *  is At's at the span's start. Every fraction of a round lies from 0 to 1, every cell is a
     *  convex combination of the span's poles, and a round adds at most ten steps u of a double at
     *  their largest coordinate M to what the cells were off by, which it carries as it was
     *  (tests/proof/section.py): each control point is within 11 n u M of the true one, and so is
     *  every point At works out by de Boor's rounds within the span, or a few steps more where all
     *  poles share a weight that its division takes out again.
     */
    struct Arch
    {
      Arm         Poles  = {};  /**< The control points.                      */
      std::size_t Degree = 0;   /**< The degree, below MostPoles.             */
      std::size_t Knot   = 0;   /**< The index of the knot that starts it.    */
      double      Most   = 0.0; /**< The largest coordinate of its poles.     */
      double      Gap    = 0.0; /**< How near the point its padded box comes. */
    };

    /**
     *  A span of a polynomial spline as a Bézier piece: the rounds at the start shared, each
     *  round's two fractions worked out once for every control point that takes it, all rounds
     *  taken in step, each cell as Section works it out.
     *
     *  @param  spline The spline, valid and polynomial, of a degree below MostPoles.
     *  @param  span   The index of the knot that starts the span, which is not empty.
     *  @return        The piece, the gap left for its box to be weighed.
     */
    auto Arched(const Spline &spline, std::size_t span) -> Arch
    {
      const auto                    degree = static_cast<std::size_t>(spline.Degree);
      const std::span<const double> knots  = spline.Knots;
      const double                  start  = knots[span];
      const double                  end    = knots[span + 1];
      Arch arch = {.Poles = {}, .Degree = degree, .Knot = span, .Most = 0.0, .Gap = 0.0};
      /**
       *  The first chain takes every round at the span's start; the one split from it before round
       *  r takes that round and the rest at its end, and gives control point n + 1 - r.
       */
      std::array<Arm, MostPoles> chains = {};
      for (std::size_t index = 0; index <= degree; ++index)
      {
        const Point pole    = spline.Poles[span - degree + index];
        chains[0].at(index) = pole;
        arch.Most           = std::max({arch.Most, std::abs(pole.X), std::abs(pole.Y)});
      }
      const auto mix = [](Arm &cells, std::size_t index, double fraction) -> void
      {
        const double rest = 1.0 - fraction;
        const Point  lhs  = cells.at(index - 1);
        const Point  rhs  = cells.at(index);
        cells.at(index)   = {.X = (lhs.X * rest) + (rhs.X * fraction),
                             .Y = (lhs.Y * rest) + (rhs.Y * fraction)};
      };
      for (std::size_t round = 1; round <= degree; ++round)
      {
        chains.at(round) = chains[0];
        for (std::size_t index = degree; index >= round; --index)
        {
          const double from  = knots[span - degree + index];
          const double until = knots[span + 1 + index - round];
          mix(chains[0], index, (start - from) / (until - from));
          const double toward = (end - from) / (until - from);
          for (std::size_t chain = 1; chain <= round; ++chain)
          {
            mix(chains.at(chain), index, toward);
          }
        }
      }
      arch.Poles[0] = chains[0].at(degree);
      for (std::size_t round = 1; round <= degree; ++round)
      {
        arch.Poles.at(degree + 1 - round) = chains.at(round).at(degree);
      }
      return arch;
    }

    /**
     *  The gap from a point to the box of some control points, padded.
     *
     *  @param  poles The control points.
     *  @param  point The point.
     *  @param  pad   The pad.
     *  @return       The gap.
     */
    auto PaddedGap(std::span<const Point> poles, Point point, double pad) -> double
    {
      Box box;
      for (const Point pole : poles)
      {
        box = Grown(box, pole);
      }
      return Outlying({.Min = {.X = box.Min.X - pad, .Y = box.Min.Y - pad},
                       .Max = {.X = box.Max.X + pad, .Y = box.Max.Y + pad}},
                      point);
    }

    /**
     *  A point asked about a piece, and how far a piece is looked into.
     */
    struct Probing
    {
      Point       At     = {};  /**< The point.                                           */
      double      Least  = 0.0; /**< The length past which a gap is not looked into.      */
      std::size_t Degree = 0;   /**< The piece's degree.                                  */
      double      Unit   = 0.0; /**< n u M, M the largest coordinate of the span's poles. */
    };

    /**
     *  How near a point a part of a span comes at the least, and every point At gives of it: the
     *  gap from the point to the box of its control points, padded by 32 n u M, more than twice
     *  their slack and At's, and 8 n u M more a halving, where that gap is more than the length or
     *  the part has been halved Finest times; else the nearer of its halves'. Each halving, by de
     *  Casteljau's averages, moves a control point by n u M at most.
     *
     *  @param  poles The part's control points.
     *  @param  probe The point, the length and the piece.
     *  @param  depth How many halvings made the part.
     *  @return       The gap, at the least.
     */
    auto Spaced(const Arm &poles, const Probing &probe, int depth) -> double
    {
      const std::size_t degree = probe.Degree;
      const double      pad = ((32.0 + (8.0 * depth)) * probe.Unit * (1.0 + 0x1p-40)) + 0x1p-1000;
      const double      gap = PaddedGap(std::span{poles}.first(degree + 1), probe.At, pad);
      if (!(gap <= probe.Least * (1.0 + 0x1p-40)) || depth == Finest)
      {
        return gap;
      }
      Arm work         = poles;
      Arm before       = {};
      Arm after        = {};
      before[0]        = work[0];
      after.at(degree) = work.at(degree);
      for (std::size_t round = 1; round <= degree; ++round)
      {
        for (std::size_t index = 0; index + round <= degree; ++index)
        {
          work.at(index) = {.X = (work.at(index).X + work.at(index + 1).X) * 0.5,
                            .Y = (work.at(index).Y + work.at(index + 1).Y) * 0.5};
        }
        before.at(round)         = work[0];
        after.at(degree - round) = work.at(degree - round);
      }
      const double first = Spaced(before, probe, depth + 1);
      return first > probe.Least * (1.0 + 0x1p-40)
               ? std::min(first, Spaced(after, probe, depth + 1))
               : first;
    }

    /**
     *  A polynomial in Bernstein form of a degree under twice MostPoles: its coefficients, the
     *  first degree and one used.
     */
    using Lean = std::array<double, 2 * MostPoles>;

    /**
     *  The binomial coefficients of every degree under twice MostPoles, each row as the kernel's
     *  table works it out, by the multiplicative formula: whole numbers, exact.
     */
    constexpr std::array<Lean, 2 * MostPoles> Choices = [] -> std::array<Lean, 2 * MostPoles>
    {
      std::array<Lean, 2 * MostPoles> table = {};
      for (std::size_t degree = 0; degree < 2 * MostPoles; ++degree)
      {
        Lean &row = table.at(degree);
        row.at(0) = 1.0;
        for (std::size_t index = 1; index <= degree; ++index)
        {
          row.at(index) = row.at(index - 1)
                        * static_cast<double>(degree - index + 1)
                        / static_cast<double>(index);
        }
      }
      return table;
    }();

    /**
     *  A coordinate of a piece less the point's, and its derivative, as polynomials.
     */
    struct Strand
    {
      Lean Value = {}; /**< The coordinate less the point's. */
      Lean Slope = {}; /**< Its derivative.                  */
    };

    /**
     *  The product of a polynomial and its derivative, and of their magnitudes, as the kernel's
     *  product of two series works each out: every pair of coefficients times their binomials, in
     *  that order, summed, then over the product's binomial.
     *
     *  @param  leg    The polynomial and its derivative.
     *  @param  degree The polynomial's degree, 1 at least.
     *  @return        The two products, the second of the magnitudes.
     */
    auto Multiplied(const Strand &leg, std::size_t degree) -> std::array<Lean, 2>
    {
      const Lean         &left    = Choices.at(degree);
      const Lean         &right   = Choices.at(degree - 1);
      const Lean         &both    = Choices.at((2 * degree) - 1);
      std::array<Lean, 2> product = {};
      for (std::size_t first = 0; first <= degree; ++first)
      {
        for (std::size_t second = 0; second < degree; ++second)
        {
          const double scale             = left.at(first) * right.at(second);
          product[0].at(first + second) += scale * leg.Value.at(first) * leg.Slope.at(second);
          product[1].at(first + second) +=
              scale * std::abs(leg.Value.at(first)) * std::abs(leg.Slope.at(second));
        }
      }
      for (std::size_t index = 0; index < 2 * degree; ++index)
      {
        product[0].at(index) /= both.at(index);
        product[1].at(index) /= both.at(index);
      }
      return product;
    }

    /**
     *  The feet of a polynomial piece of weights 1, as the kernel's Perpendiculars finds them, step
     *  for step: U U' + V V' for U and V its coordinates less the point's, its bound through their
     *  magnitudes, nothing where every coefficient is within its rounding of 0, else the roots the
     *  kernel finds of it; held in a few doubles, not in series of the kernel's every degree.
     *
     *  @param  arch  The piece.
     *  @param  point The point.
     *  @return       The feet's parameters, from 0 to 1, in order.
     */
    auto Plumbed(const Arch &arch, Point point) -> List<double>
    {
      const std::size_t     degree = arch.Degree;
      const auto            scale  = static_cast<double>(degree);
      std::array<Strand, 2> legs   = {};
      for (std::size_t index = 0; index <= degree; ++index)
      {
        legs[0].Value.at(index) = arch.Poles.at(index).X - point.X;
        legs[1].Value.at(index) = arch.Poles.at(index).Y - point.Y;
      }
      for (Strand &leg : legs)
      {
        for (std::size_t index = 0; index < degree; ++index)
        {
          leg.Slope.at(index) = scale * (leg.Value.at(index + 1) - leg.Value.at(index));
        }
      }
      const std::array<Lean, 2> across  = Multiplied(legs[0], degree);
      const std::array<Lean, 2> rising  = Multiplied(legs[1], degree);
      Lean                      turning = across[0];
      Lean                      bound   = across[1];
      bool                      flat    = true;
      for (std::size_t index = 0; index < 2 * degree; ++index)
      {
        turning.at(index) += rising[0].at(index);
        bound.at(index)   += rising[1].at(index);
        flat               = flat
                          && std::abs(turning.at(index))
                                 <= 64.0 * std::numeric_limits<double>::epsilon() * bound.at(index);
      }
      return flat ? List<double>{}
                  : Kernel::Roots(std::span{turning}.first(2 * degree), Kernel::Budget);
    }

    /**
     *  The place nearest a point so far, as Least chooses among places: the nearest, of those as
     *  near the first along the curve, whatever the order they come in, none of them at no number
     *  of distance.
     */
    struct Weighing
    {
      Point  From  = {};                                      /**< The point.          */
      Place  Best  = {};                                      /**< The nearest so far. */
      double Least = std::numeric_limits<double>::infinity(); /**< How near it is.     */

      /**
       *  A place weighed.
       *
       *  @param place The place.
       */
      void Weigh(const Place &place)
      {
        const double away = Distance(From, place.At);
        if (away < Least || (Exactly(away, Least) && place.Along < Best.Along))
        {
          Best  = place;
          Least = away;
        }
      }

      /**
       *  Whether a gap, at the least, lies farther than the nearest place so far, as every place it
       *  bounds does.
       *
       *  @param  gap The gap.
       *  @return     Whether it does.
       */
      [[nodiscard]] auto Beyond(double gap) const -> bool { return gap > Least * (1.0 + 0x1p-40); }
    };

    /**
     *  A polynomial spline's spans as Bézier pieces, each with how near the point its padded box
     *  comes.
     *
     *  @param  spline The spline, valid and polynomial, of a degree below MostPoles.
     *  @param  point  The point.
     *  @return        The pieces, in order along the spline.
     */
    auto Arches(const Spline &spline, Point point) -> List<Arch>
    {
      const auto                    degree = static_cast<std::size_t>(spline.Degree);
      const std::span<const double> knots  = spline.Knots;
      List<Arch>                    arches;
      for (std::size_t span = degree; span < spline.Poles.size(); ++span)
      {
        if (knots[span] < knots[span + 1])
        {
          Arch         arch = Arched(spline, span);
          const double pad =
              (32.0 * static_cast<double>(degree) * 0x1p-53 * arch.Most * (1.0 + 0x1p-40))
              + 0x1p-1000;
          arch.Gap = PaddedGap(std::span{arch.Poles}.first(degree + 1), point, pad);
          arches.push_back(arch);
        }
      }
      return arches;
    }

    /**
     *  A span's feet, as the spline's parameters, each kept once as Placed keeps them: a spline of
     *  no weights by its piece's lean feet, any other by the kernel's on its section.
     *
     *  @param  spline The spline, valid and polynomial.
     *  @param  arch   The span as a Bézier piece.
     *  @param  point  The point.
     *  @return        The parameters; nothing where one lies within a few steps of a double of the
     *                 span's end, where Placed may keep it once with one from the next span, or a
     *                 seam's.
     */
    auto Rooted(const Spline &spline, const Arch &arch, Point point) -> std::optional<List<double>>
    {
      const Range range = {.Start = spline.Knots[arch.Knot], .End = spline.Knots[arch.Knot + 1]};
      const List<double> roots =
          spline.Weights.empty()
              ? Plumbed(arch, point)
              : Kernel::Perpendiculars(Kernel::Section(spline, arch.Knot), point, Kernel::Budget)
                    .Roots;
      List<double> alongs;
      alongs.reserve(roots.size());
      for (const double root : roots)
      {
        const double along = Kernel::Across(range, root);
        const double close =
            256.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(along));
        if (along <= range.Start + close || along >= range.End - close)
        {
          return std::nullopt;
        }
        alongs.push_back(along);
      }
      Gathered(alongs);
      return alongs;
    }

    /**
     *  A span of a rational spline weighed by branch and bound: the knot that starts it, its piece,
     *  and how near the point the padded box of its control points comes.
     */
    struct Lot
    {
      std::size_t   Knot  = 0;   /**< The index of the knot that starts it.  */
      Kernel::Piece Piece = {};  /**< Its piece, as Kernel::Section cuts it. */
      double        Gap   = 0.0; /**< How near the point its box comes.      */
    };

    /**
     *  What a rational span's pieces are weighed with: the point, the length past which a gap is
     *  not looked into, and a step of a double at the span's size, (w_max / w_min) u M.
     */
    struct Sounding
    {
      Point  At    = {};  /**< The point.                                            */
      double Least = 0.0; /**< The length past which a gap is not looked into.       */
      double Unit  = 0.0; /**< (w_max / w_min) u M, M the span's largest coordinate. */
    };

    /**
     *  How near a point a part of a rational span comes at the least, every point At gives of it
     *  too: the gap to the box of its piece's control points, each its homogeneous point over its
     *  weight, padded; where that is no more than the length and the part has been halved fewer
     *  than Finest times, the nearer of its halves'.
     *
     *  Its weights are above nothing, so the part lies within the hull of those points. De Boor's
     *  mixes, de Casteljau's halvings and the quotients put the piece's points and At's within some
     *  (4 d + 2) steps of a double of the span's largest coordinate M of the curve's, times w_max /
     *  w_min, a halving adding its own few, and the fractions' rounding moves At's point by 12.5 d
     *  (w_max / w_min) steps of the poles' spread at most (tests/proof/pads.py, as Straying in
     *  extents.cpp takes it): 128 (d + 1) such steps, and 16 more a halving, cover them all several
     *  times over.
     *
     *  @param  piece The part.
     *  @param  sound The point, the length and the step.
     *  @param  depth How many halvings made the part.
     *  @return       The gap, at the least.
     */
    auto Sounded(const Kernel::Piece &piece, const Sounding &sound, int depth) -> double
    {
      Box box;
      for (const Homogeneous &pole : Kernel::Poles(piece))
      {
        box = Grown(box, Projected(pole));
      }
      const double pad = (((128.0 * static_cast<double>(piece.Degree + 1)) + (16.0 * depth))
                          * sound.Unit
                          * (1.0 + 0x1p-40))
                       + 0x1p-1000;
      const double gap = Outlying({.Min = {.X = box.Min.X - pad, .Y = box.Min.Y - pad},
                                   .Max = {.X = box.Max.X + pad, .Y = box.Max.Y + pad}},
                                  sound.At);
      if (!(gap <= sound.Least * (1.0 + 0x1p-40)) || depth == Finest)
      {
        return gap;
      }
      const std::array<Kernel::Piece, 2> halves = Kernel::Split(piece, 0.5);
      const double                       first  = Sounded(halves[0], sound, depth + 1);
      return first > sound.Least * (1.0 + 0x1p-40)
               ? std::min(first, Sounded(halves[1], sound, depth + 1))
               : first;
    }

    /**
     *  A step of a double at a rational span's size: (w_max / w_min) u M, M the largest coordinate
     *  of its poles and w its weights.
     *
     *  @param  spline The spline, valid and rational.
     *  @param  span   The index of the knot that starts the span.
     *  @return        The step.
     */
    auto Stepping(const Spline &spline, std::size_t span) -> double
    {
      const auto degree = static_cast<std::size_t>(spline.Degree);
      double     most   = 0.0;
      double     light  = std::numeric_limits<double>::infinity();
      double     heavy  = 0.0;
      for (std::size_t index = span - degree; index <= span; ++index)
      {
        const Point pole = spline.Poles[index];
        most             = std::max({most, std::abs(pole.X), std::abs(pole.Y)});
        light            = std::min(light, spline.Weights[index]);
        heavy            = std::max(heavy, spline.Weights[index]);
      }
      return (heavy / light) * 0x1p-53 * most;
    }

    /**
     *  Whether a rational spline's spans are even, each piece's weights within Gentle times the
     *  lighter of its end weights: else Paved cuts them, and Listed finds the feet on its paving.
     *
     *  @param  lots The spans, their pieces cut.
     *  @return      Whether they are.
     */
    auto Unpaved(std::span<const Lot> lots) -> bool
    {
      return std::ranges::all_of(
          lots,
          [](const Lot &lot) -> bool
          {
            const std::span<const Homogeneous> poles = Kernel::Poles(lot.Piece);
            const double                       light = std::min(poles.front().W, poles.back().W);
            return std::ranges::all_of(poles, [light](const Homogeneous &pole) -> bool
                                       { return pole.W <= Gentle * light; });
          });
    }
  } // namespace

  auto Pruned(const Spline &spline, Point point) -> std::optional<Place>
  {
    const auto degree = static_cast<std::size_t>(spline.Degree);
    const auto within = [](Point corner) -> bool
    {
      return std::abs(corner.X) <= Widest && std::abs(corner.Y) <= Widest;
    };
    if (degree >= MostPoles || !within(point) || !std::ranges::all_of(spline.Poles, within))
    {
      return std::nullopt;
    }
    const bool                    plain    = spline.Weights.empty();
    const std::span<const double> knots    = spline.Knots;
    List<Arch>                    arches   = Arches(spline, point);
    Weighing                      weighing = {.From = point};
    for (const Arch &arch : arches)
    {
      const double start = knots[arch.Knot];
      weighing.Weigh({.At = plain ? arch.Poles[0] : At(spline, start), .Along = start});
    }
    const Arch  &last = arches.back();
    const double end  = knots[last.Knot + 1];
    weighing.Weigh({.At = plain ? last.Poles.at(degree) : At(spline, end), .Along = end});
    std::ranges::sort(arches, {}, &Arch::Gap);
    for (const Arch &arch : arches)
    {
      if (weighing.Beyond(arch.Gap))
      {
        break;
      }
      const Probing probe = {.At     = point,
                             .Least  = weighing.Least,
                             .Degree = degree,
                             .Unit   = static_cast<double>(degree) * 0x1p-53 * arch.Most};
      if (weighing.Beyond(Spaced(arch.Poles, probe, 0)))
      {
        continue;
      }
      const std::optional<List<double>> alongs = Rooted(spline, arch, point);
      if (!alongs)
      {
        return std::nullopt;
      }
      for (const double along : *alongs)
      {
        weighing.Weigh({.At = At(spline, along), .Along = along});
      }
    }
    return weighing.Best;
  }

  auto Rationally(const Spline &spline, Point point) -> std::optional<Place>
  {
    const auto degree = static_cast<std::size_t>(spline.Degree);
    const auto within = [](Point corner) -> bool
    {
      return std::abs(corner.X) <= Widest && std::abs(corner.Y) <= Widest;
    };
    const auto [light, heavy] = std::ranges::minmax(spline.Weights);
    if (degree >= MostPoles
        || !within(point)
        || !std::ranges::all_of(spline.Poles, within)
        || !(heavy <= 0x1p40 * light))
    {
      return std::nullopt;
    }
    const std::span<const double> knots = spline.Knots;
    List<Lot>                     lots;
    for (std::size_t span = degree; span < spline.Poles.size(); ++span)
    {
      if (knots[span] < knots[span + 1])
      {
        lots.push_back({.Knot = span, .Piece = Kernel::Section(spline, span), .Gap = 0.0});
      }
    }
    if (!Unpaved(lots))
    {
      return std::nullopt;
    }
    Weighing weighing = {.From = point};
    for (Lot &lot : lots)
    {
      weighing.Weigh({.At = At(spline, knots[lot.Knot]), .Along = knots[lot.Knot]});
      lot.Gap = Sounded(lot.Piece, {.At = point, .Least = 0.0, .Unit = Stepping(spline, lot.Knot)},
                        Finest);
    }
    const double end = knots[lots.back().Knot + 1];
    weighing.Weigh({.At = At(spline, end), .Along = end});
    std::ranges::sort(lots, {}, &Lot::Gap);
    for (const Lot &lot : lots)
    {
      if (weighing.Beyond(lot.Gap))
      {
        break;
      }
      const Sounding sound = {
          .At = point, .Least = weighing.Least, .Unit = Stepping(spline, lot.Knot)};
      if (weighing.Beyond(Sounded(lot.Piece, sound, 0)))
      {
        continue;
      }
      const Range        range = {.Start = knots[lot.Knot], .End = knots[lot.Knot + 1]};
      const List<double> roots =
          Kernel::Perpendiculars({.Shape = lot.Piece, .Speed = Kernel::Derived(spline, lot.Knot)},
                                 point, Kernel::Budget)
              .Roots;
      List<double> alongs;
      for (const double root : roots)
      {
        const double along = Kernel::Across(range, root);
        const double close =
            256.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(along));
        if (along <= range.Start + close || along >= range.End - close)
        {
          return std::nullopt;
        }
        alongs.push_back(along);
      }
      Gathered(alongs);
      for (const double along : alongs)
      {
        weighing.Weigh({.At = At(spline, along), .Along = along});
      }
    }
    return weighing.Best;
  }

  auto CloserBy(Vector sizes, double near, double far) -> double
  {
    const double middle = (near + far) / 2.0;
    return (2.0 / ((1.0 / sizes.X) + (1.0 / sizes.Y)))
         * (4.0 / (Pi * Pi))
         * ((far - near) / 2.0)
         * std::min(middle, Pi - middle);
  }
} // namespace Geoxy::Near
