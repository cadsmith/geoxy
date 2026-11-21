/***************************************************************************************************
 *  @file      source/swift/pieced.cpp
 *  @brief     Polynomial splines met piece by piece by the swift kernels, their crossings told
 *             apart from one another as the thorough way tells them.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "swift/pieced.hpp"

#include "real.hpp"
#include "scratch.hpp"
#include "sure.hpp"
#include "swift/lane.hpp"
#include "swift/section.hpp"
#include "swift/twin.hpp"

#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <variant>

namespace Geoxy::Swift
{
  namespace
  {
    /**
     *  How far apart, in reaches, two pieces' boxes may lie and still be met: twice as far as the
     *  twin kernel keeps what does not cross, so that pieces not met lie farther apart than every
     *  length its bands are complete to.
     */
    constexpr double ClearReaches = 128.0;

    /**
     *  The most a spline's derivative by its own parameter may be on a span, as Derivative works it
     *  out, and twice that: de Boor's on its slopes, each a pole's step times its degree over the
     *  knots it spans, a convex combination of those, whose coordinates are no larger than theirs.
     *
     *  @param  spline The spline.
     *  @param  knot   The index of the knot that starts the span.
     *  @return        The pace.
     */
    auto Treaded(const Spline &spline, std::size_t knot) -> double
    {
      const auto                    count  = static_cast<std::size_t>(spline.Degree);
      const auto                    degree = static_cast<double>(count);
      const std::span<const double> knots  = spline.Knots;
      double                        most   = 0.0;
      for (std::size_t index = knot - count; index < knot; ++index)
      {
        const Vector step  = spline.Poles[index + 1] - spline.Poles[index];
        const double width = knots[index + count + 1] - knots[index + 1];
        most =
            std::max({most, degree * std::abs(step.X) / width, degree * std::abs(step.Y) / width});
      }
      return 4.0 * most;
    }

    /**
     *  The most a curve's derivative by its own parameter may be on a section, and twice that: a
     *  Bézier curve's de Casteljau's on n times its control steps, a convex combination of those; a
     *  spline's span's as Treaded has it.
     *
     *  @param  section The section.
     *  @return         The pace.
     */
    auto Pace(const Section &section) -> double
    {
      if (section.Source != nullptr)
      {
        return Treaded(*section.Source, section.Knot);
      }
      const auto degree = static_cast<double>(section.Degree);
      double     most   = 0.0;
      for (std::size_t index = 0; index < section.Degree; ++index)
      {
        const Vector step = section.Poles.at(index + 1) - section.Poles.at(index);
        most              = std::max({most, degree * std::abs(step.X), degree * std::abs(step.Y)});
      }
      return 4.0 * most;
    }

    /**
     *  A broad section's pace, a spline's span's: a rational one's 2 n w_max / w_min times its
     *  hull's diagonal over its width, twice; a polynomial one's as Treaded has it.
     *
     *  @param  section The section, a span of a spline.
     *  @return         The pace.
     */
    auto Gaited(const Broad &section) -> double
    {
      return section.Rated ? 4.0 * Fastest(section) / (section.Span.End - section.Span.Start)
                           : Treaded(*section.Source, section.Knot);
    }

    /**
     *  The largest magnitude of a section's parameter: that of either end of its span.
     *
     *  @tparam SectionT The kind of section.
     *  @param  section  The section.
     *  @return          The magnitude.
     */
    template <typename SectionT>
    auto Farthest(const SectionT &section) -> double
    {
      return std::max(std::abs(section.Span.Start), std::abs(section.Span.End));
    }

    /**
     *  Whether the thorough way holds a crossing's candidates to the tolerance alone: where either
     *  curve is a spline it holds them to twice what a step of each place moves its curve where
     *  that is more, and here it is less, each place within its section's span, its step no more
     *  than the step at twice the span's farthest end.
     *
     *  @param  sections  The sections the crossing is on.
     *  @param  tolerance The tolerance.
     *  @return           Whether it does.
     */
    auto Strict(std::array<const Section *, 2> sections, const Tolerance &tolerance) -> bool
    {
      const double coarse = 2.0
                          * ((Pace(*sections[0]) * Above(2.0 * Farthest(*sections[0])))
                             + (Pace(*sections[1]) * Above(2.0 * Farthest(*sections[1]))));
      return 2.0 * coarse <= tolerance.Length;
    }

    /**
     *  A curve as the twin kernel's pieces: a polynomial spline's spans, a Bézier curve whole.
     *
     *  @param  curve The curve.
     *  @return       Its pieces; nothing for any other curve, or a spline the kernels do not take.
     */
    auto Stranded(const Curve &curve) -> std::optional<List<Section>>
    {
      if (const auto *spline = std::get_if<Spline>(&curve))
      {
        return Sectioned(*spline);
      }
      if (const auto *quadratic = std::get_if<Bezier::Quadratic>(&curve))
      {
        return List<Section>{Sectioned(*quadratic)};
      }
      if (const auto *cubic = std::get_if<Bezier::Cubic>(&curve))
      {
        return List<Section>{Sectioned(*cubic)};
      }
      return std::nullopt;
    }

    /**
     *  A crossing the twin kernel proved, its bands as places of each curve.
     */
    struct Proof
    {
      Twain                Found = {}; /**< The crossing, as proved. */
      std::array<Range, 2> Near  = {}; /**< Each curve's places where the curves may come
                                            within the tolerance about it.                 */
      std::array<Range, 2> Twice = {}; /**< Likewise within twice the tolerance. */
    };

    /**
     *  A band of a section's places as places of its curve, rounded outward: the span's start and
     *  its width times each end, each off by a few steps of a double at the largest of them.
     *
     *  @param  section The section.
     *  @param  band    The band, from 0 to 1.
     *  @return         The band on the curve.
     */
    auto Widened(const Section &section, Range band) -> Range
    {
      const double start = section.Span.Start;
      const double width = section.Span.End - start;
      const double low   = start + (width * band.Start);
      const double high  = start + (width * band.End);
      const double pad = 4.0 * Above(std::max({Farthest(section), std::abs(low), std::abs(high)}));
      return {.Start = low - pad, .End = high + pad};
    }

    /**
     *  How far apart two stretches are.
     *
     *  @param  one   One.
     *  @param  other The other.
     *  @return       The gap between them: above zero where they do not meet.
     */
    auto Between(Range one, Range other) -> double
    {
      return std::max(other.Start - one.End, one.Start - other.End);
    }

    /**
     *  The gap between two boxes, along x or along y, whichever is wider: nothing where they meet.
     *
     *  @param  one   One box.
     *  @param  other The other.
     *  @return       The gap.
     */
    auto Spaced(const Box &one, const Box &other) -> double
    {
      return std::max({one.Min.X - other.Max.X, other.Min.X - one.Max.X, one.Min.Y - other.Max.Y,
                       other.Min.Y - one.Max.Y, 0.0});
    }

    /**
     *  Which curve the thorough way looks along between two hits, and that curve's domain.
     */
    struct Looking
    {
      std::size_t Which  = 0;  /**< The curve, 0 the first. */
      Range       Domain = {}; /**< Its domain.             */
    };

    /**
     *  Whether the thorough way, looking halfway between two crossings along a curve, finds it
     *  farther than twice the tolerance from the other, and so keeps them two hits: the place it
     *  looks at lies between the middles of their near bands, as each hit's place lies in its band,
     *  or half the curve's domain from there, the short way round a curve that closes; and the
     *  curves come within twice the tolerance only in the crossings' bands of it, as the twin
     *  kernel's bands show for every pair of pieces met, and pieces not met lie farther apart.
     *
     *  @param  proofs The crossings.
     *  @param  pair   The two, by index.
     *  @param  side   The curve looked along.
     *  @return        Whether it does.
     */
    auto Parted(std::span<const Proof> proofs, std::array<std::size_t, 2> pair, const Looking &side)
        -> bool
    {
      const std::size_t which = side.Which;
      const Range       one   = proofs[pair[0]].Near.at(which);
      const Range       other = proofs[pair[1]].Near.at(which);
      const double      whole = side.Domain.End - side.Domain.Start;
      const double      far =
          std::max({std::abs(one.Start), std::abs(one.End), std::abs(other.Start),
                    std::abs(other.End), std::abs(side.Domain.Start), std::abs(side.Domain.End)});
      const double pad  = 8.0 * Above(2.0 * far);
      const Range  half = {.Start = ((one.Start + other.Start) / 2.0) - pad,
                           .End   = ((one.End + other.End) / 2.0) + pad};
      for (const double shift : {-whole / 2.0, 0.0, whole / 2.0})
      {
        const Range moved = {.Start = half.Start + shift - pad, .End = half.End + shift + pad};
        for (const Proof &proof : proofs)
        {
          if (!(Between(moved, proof.Twice.at(which)) > 0.0))
          {
            return false;
          }
        }
      }
      return true;
    }

    /**
     *  Whether the thorough way keeps every crossing a hit of its own: every two farther apart than
     *  the tolerance, each hit's point being within its band's run of its crossing, and told apart
     *  by the place halfway between them along one curve or the other.
     *
     *  @param  proofs  The crossings.
     *  @param  domains The curves' domains.
     *  @return         Whether it does.
     */
    auto Told(std::span<const Proof> proofs, const std::array<Range, 2> &domains) -> bool
    {
      for (std::size_t one = 0; one < proofs.size(); ++one)
      {
        for (std::size_t other = one + 1; other < proofs.size(); ++other)
        {
          const Twain &first  = proofs[one].Found;
          const Twain &second = proofs[other].Found;
          const double far =
              (std::max(first.Within, second.Within) + first.Runs + second.Runs) * Sure::Grow;
          const bool parted = Parted(proofs, {one, other}, {.Which = 0, .Domain = domains[0]})
                           || Parted(proofs, {one, other}, {.Which = 1, .Domain = domains[1]});
          if (!(Distance(first.At, second.At) > far) || !parted)
          {
            return false;
          }
        }
      }
      return true;
    }
  } // namespace

  auto Pieced(const Curves &curves, const Tolerance &tolerance, Intersections &into)
      -> std::optional<bool>
  {
    const Curve &first  = curves[0];
    const Curve &second = curves[1];
    if (!std::holds_alternative<Spline>(first) && !std::holds_alternative<Spline>(second))
    {
      return std::nullopt;
    }
    if (Laned(first) || Laned(second))
    {
      return std::nullopt;
    }
    const std::optional<List<Section>> ones   = Stranded(first);
    const std::optional<List<Section>> others = ones ? Stranded(second) : std::nullopt;
    if (!others)
    {
      return false;
    }
    const double clear = ClearReaches * Reach(tolerance);
    List<Proof>  proofs;
    for (const Section &one : *ones)
    {
      for (const Section &other : *others)
      {
        if (Spaced(one.Hull, other.Hull) > clear)
        {
          continue;
        }
        const std::optional<Twains> twains = Twined({&one, &other}, tolerance);
        if (!twains || (!twains->empty() && !Strict({&one, &other}, tolerance)))
        {
          return false;
        }
        for (const Twain &twain : *twains)
        {
          proofs.push_back(
              {.Found = twain,
               .Near  = {Widened(one, twain.Near[0]), Widened(other, twain.Near[1])},
               .Twice = {Widened(one, twain.Twice[0]), Widened(other, twain.Twice[1])}});
        }
      }
    }
    if (!Told(proofs, {Domain(first), Domain(second)}))
    {
      return false;
    }
    for (const Proof &proof : proofs)
    {
      const Twain &found = proof.Found;
      into.Hits.push_back({.At     = found.At,
                           .First  = {.At = At(first, found.Along[0]), .Along = found.Along[0]},
                           .Second = {.At = At(second, found.Along[1]), .Along = found.Along[1]},
                           .Kind   = Contact::Cross});
    }
    std::ranges::sort(into.Hits, {}, [](const Hit &hit) -> double { return hit.First.Along; });
    return true;
  }

  auto Fine(const Broad &section, double stride, const Tolerance &tolerance) -> bool
  {
    const double coarse = 2.0 * (stride + (Gaited(section) * Above(2.0 * Farthest(section))));
    return 2.0 * coarse <= tolerance.Length;
  }

  auto Strung(const Lane &lane, const Spline &spline, const Tolerance &tolerance)
      -> std::optional<Spots>
  {
    const std::optional<List<Section>> sections = Sectioned(spline);
    if (!sections)
    {
      return Skeined(lane, spline, tolerance);
    }
    return Strung(lane, std::span<const Section>{*sections}, tolerance);
  }

  auto Strung(const Lane &lane, std::span<const Section> sections, const Tolerance &tolerance)
      -> std::optional<Spots>
  {
    /**
     *  A line's box is its ends'; a ray's or an endless line's runs past every double, and every
     *  piece is met.
     */
    const bool   bounded = lane.Count == 2;
    const Box    reach   = Grown(Grown(Box{}, lane.Ends[0]), lane.Ends[1]);
    const double clear   = ClearReaches * Reach(tolerance);
    Spots        spots;
    for (const Section &section : sections)
    {
      if (bounded && Spaced(section.Hull, reach) > clear)
      {
        continue;
      }
      const std::optional<Spots> found = Bernstein(lane, section, tolerance);
      if (!found)
      {
        return std::nullopt;
      }
      for (const Spot &spot : *found)
      {
        /**
         *  The thorough way holds the crossing's candidates to the tolerance alone where twice what
         *  a step of each place moves its curve is less: a line's place within its band of the one
         *  found, so no farther from nothing than twice it and one.
         */
        const double line   = 2.0 * lane.Speed * Above((2.0 * std::abs(spot.Along[0])) + 1.0);
        const double coarse = 2.0 * (line + (Pace(section) * Above(2.0 * Farthest(section))));
        if (!(2.0 * coarse <= tolerance.Length))
        {
          return std::nullopt;
        }
        spots.push_back(spot);
      }
    }
    return spots;
  }
} // namespace Geoxy::Swift
