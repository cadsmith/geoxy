/***************************************************************************************************
 *  @file      source/intersect.cpp
 *  @brief     Where two curves meet: directly where their kinds and places leave nothing to decide,
 *             the thorough way otherwise; curves too large, too small or too stiff to meet as they
 *             are met through others that stand for them.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/intersect.hpp"

#include "crossing.hpp"
#include "crossings/self.hpp"
#include "eased.hpp"
#include "met/linear.hpp"
#include "met/round.hpp"
#include "met/slight.hpp"
#include "piecewise.hpp"
#include "sane.hpp"
#include "swift/met.hpp"
#include "thorough.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/line.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/xline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  Whether either of two curves is a ray or an endless line.
     *
     *  @param  curves The two curves.
     *  @return        Whether one is.
     */
    auto Unending(const Curves &curves) -> bool
    {
      const auto endless = [](const Curve &curve) -> bool
      {
        return std::holds_alternative<Ray>(curve) || std::holds_alternative<Xline>(curve);
      };
      return endless(curves[0]) || endless(curves[1]);
    }

    /**
     *  Whether two curves are a pair no swift kernel meets: each a line, a ray, an endless line or
     *  an arc, one of them an arc. Swift::Met finds nothing for them, and so far into Directed no
     *  step has touched what they meet in, which came empty: Straightly takes no arc, Piecewise no
     *  curve but a polyline, and Concentrically writes only where it decides.
     *
     *  @param  curves The two curves.
     *  @return        Whether they are.
     */
    auto Swiftless(const Curves &curves) -> bool
    {
      const auto plain = [](const Curve &curve) -> bool
      {
        return std::holds_alternative<Line>(curve)
            || std::holds_alternative<Ray>(curve)
            || std::holds_alternative<Xline>(curve)
            || std::holds_alternative<Arc>(curve);
      };
      return plain(curves[0])
          && plain(curves[1])
          && (std::holds_alternative<Arc>(curves[0].get())
              || std::holds_alternative<Arc>(curves[1].get()));
    }

    /**
     *  Whether either of two curves is a spline.
     *
     *  @param  curves The two curves.
     *  @return        Whether one is.
     */
    auto Splinish(const Curves &curves) -> bool
    {
      return std::holds_alternative<Spline>(curves[0].get())
          || std::holds_alternative<Spline>(curves[1].get());
    }

    /**
     *  Hits put in the thorough way's order along a whole circle first: a hit at its seam, found
     *  there at its first pass and polished to a hair short of its last, first, as Seamwise puts
     *  them.
     *
     *  @param first     The first curve.
     *  @param tolerance The tolerance.
     *  @param into      The meeting, its hits in order along the first curve.
     */
    void Seamed(const Curve &first, const Tolerance &tolerance, Intersections &into)
    {
      const auto *circle = std::get_if<Arc>(&first);
      if (circle == nullptr || std::abs(circle->Sweep.InRadians()) < 2.0 * Pi)
      {
        return;
      }
      Met::Seamwise(*circle, tolerance, into);
    }

    /**
     *  The crossings a direct way finds, roughly, of a straight curve and an arc, either first, or
     *  of two arcs not of one circle: each crossing's places as the curves come, the first's first,
     *  built where it is kept, so that nothing of it is copied.
     *
     *  @param  curves    The two curves.
     *  @param  tolerance The tolerance, no finer than their grain.
     *  @return           The crossings, and the ends at them; nothing for any other pair, or where
     *                    the thorough way must decide.
     */
    auto Roughly(const Curves &curves, const Tolerance &tolerance) -> std::optional<Met::Nearby>
    {
      const Curve &first  = curves[0];
      const Curve &second = curves[1];
      const auto  *round  = std::get_if<Arc>(&first);
      const auto  *across = std::get_if<Arc>(&second);
      if (round != nullptr && across != nullptr)
      {
        return Met::Circularly(*round, *across, tolerance);
      }
      if (round == nullptr && across == nullptr)
      {
        return std::nullopt;
      }
      const std::optional<Met::Straight> straight =
          Met::Straightened(round != nullptr ? second : first);
      if (!straight)
      {
        return std::nullopt;
      }
      if (across != nullptr)
      {
        return Met::Roundly(*straight, *across, tolerance);
      }
      std::optional<Met::Nearby> near = Met::Roundly(*straight, *round, tolerance);
      if (near)
      {
        for (Approach &local : std::span{near->Items}.first(near->Count))
        {
          std::swap(local.Along[0], local.Along[1]);
        }
      }
      return near;
    }

  } // namespace

  /**
   *  Crossings found approximately made hits by the exact crossing the thorough way ends in, in
   *  order along the first curve; nothing if any is not found so.
   *
   *  @param  curves The two curves.
   *  @param  near   Each crossing, roughly.
   *  @param  into   Where the hits go, empty on entry, and left so where one is not found.
   *  @return        Whether every exact crossing was found.
   */
  auto Exactly(const std::array<const Curve *, 2> &curves, std::span<const Approach> near,
               Intersections &into) -> bool
  {
    /**
     *  No crossing to make exact: the carriers' crossings are not worked out for nothing.
     */
    if (near.empty())
    {
      return true;
    }
    const Crossings::Roots roots = Crossings::Quickly(curves);
    for (const Approach &approach : near)
    {
      const std::optional<Crossings::Crossed> exact = Crossings::Crossing(curves, approach, roots);
      if (!exact)
      {
        into.Hits.clear();
        return false;
      }
      into.Hits.push_back({.At     = exact->At,
                           .First  = {.At = exact->On[0], .Along = exact->Along[0]},
                           .Second = {.At = exact->On[1], .Along = exact->Along[1]},
                           .Kind   = approach.Touch ? Contact::Touch : Contact::Cross});
    }
    std::ranges::sort(into.Hits, {}, [](const Hit &hit) -> double { return hit.First.Along; });
    return true;
  }

  namespace
  {
    /**
     *  Two curves met directly where their kinds and places leave nothing to decide, a polyline's
     *  legwork as kept where it was made ready; nothing otherwise.
     *
     *  @param  footing   The two curves, and the legwork kept for each.
     *  @param  tolerance The tolerance, no finer than the curves' grain.
     *  @param  into      Where what they meet in goes, empty on entry, and left so where the
     *                    thorough way must decide.
     *  @return           Whether it was decided directly.
     */
    auto Directed(const Stance &footing, const Tolerance &tolerance, Intersections &into) -> bool
    {
      const Curves curves = {*footing.Shapes[0], *footing.Shapes[1]};
      const Curve &first  = curves[0];
      const Curve &second = curves[1];
      /**
       *  A slight arc that stands off its chord by more than the tolerance is no direct way's: its
       *  circle's center is too far for them. The thorough way takes it as a piece.
       */
      if (Met::Remotely(first, second))
      {
        return false;
      }
      if (Met::Straightly(curves, tolerance, into) || Piecewise(footing, tolerance, into))
      {
        return true;
      }
      const auto *round  = std::get_if<Arc>(&first);
      const auto *across = std::get_if<Arc>(&second);
      if (round != nullptr
          && across != nullptr
          && Met::Concentrically({round, across}, tolerance, into))
      {
        return true;
      }
      if (!Swiftless(curves) && Swift::Met(curves, tolerance, into))
      {
        return true;
      }
      const std::optional<Met::Nearby> near = Roughly(curves, tolerance);
      const bool                       ruled =
          near.transform(
                  [&first, &second, &tolerance, &into](const Met::Nearby &places) -> bool
                  {
                    if (!Exactly({&first, &second}, places.Held(), into))
                    {
                      return false;
                    }
                    for (Hit &hit : into.Hits)
                    {
                      if (hit.Kind == Contact::Cross && places.Ending(hit.At, tolerance.Length))
                      {
                        hit.Kind = Contact::End;
                      }
                    }
                    Seamed(first, tolerance, into);
                    return true;
                  })
              .value_or(false);
      return ruled || Banded(curves, tolerance, into);
    }
  } // namespace

  auto Direct(const Curves &curves, const Tolerance &tolerance, Intersections &into) -> bool
  {
    return Directed({.Shapes = {&curves[0].get(), &curves[1].get()}}, tolerance, into);
  }

  auto Loosened(const Curves &curves, const Tolerance &tolerance, Intersections &into) -> bool
  {
    const auto eased = [](const Curve &curve) -> std::optional<Easing>
    {
      const auto *spline = std::get_if<Spline>(&curve);
      return spline != nullptr ? Eased(*spline, Mild) : std::nullopt;
    };
    const std::array<std::optional<Easing>, 2> ways = {eased(curves[0]), eased(curves[1])};
    if (!ways[0] && !ways[1])
    {
      return false;
    }
    const auto beyond = [](const std::optional<Easing> &way) -> bool
    {
      return way && way->Pieces.empty();
    };
    if (beyond(ways[0]) || beyond(ways[1]))
    {
      /**
       *  Weights farther apart than doubles hold beside each other leave no stand-in to meet
       *  through and no place to tell: nothing is found, and the meeting says a hit may be missing.
       */
      into.Spent = true;
      return true;
    }
    std::array<Curve, 2>           stood = {Curve{Line{}}, Curve{Line{}}};
    std::array<Crossings::Read, 2> reads = {
        Crossings::Read{.Shape = &curves[0].get(), .Stand = &curves[0].get(), .Way = nullptr},
        Crossings::Read{.Shape = &curves[1].get(), .Stand = &curves[1].get(), .Way = nullptr}};
    for (std::size_t which = 0; which < ways.size(); ++which)
    {
      if (const std::optional<Easing> &way = ways.at(which); way.has_value())
      {
        stood.at(which)       = Curve{way->Shape};
        reads.at(which).Stand = &stood.at(which);
        reads.at(which).Way   = std::addressof(*way);
      }
    }
    /**
     *  The stand-ins lie where the splines do, which were found within reach of each other; and the
     *  direct ways are for lines and arcs.
     */
    Thorough(Curves{*reads[0].Stand, *reads[1].Stand}, tolerance, into);
    const auto back = [&reads](std::size_t which, double along) -> double
    {
      const Easing *way = reads.at(which).Way;
      return way != nullptr ? Homed(*way, along) : along;
    };
    for (Hit &hit : into.Hits)
    {
      /**
       *  A crossing was polished on the stand-ins, whose poles are the splines' rounded once: it is
       *  polished again from the splines' own numbers, at the places the stand-ins tell.
       */
      const std::optional<Crossings::Crossed> exact =
          hit.Kind == Contact::Cross
              ? Crossings::Retold(reads, {hit.First.Along, hit.Second.Along}, hit.At)
              : std::nullopt;
      hit.At           = exact ? exact->At : hit.At;
      hit.First.Along  = back(0, hit.First.Along);
      hit.Second.Along = back(1, hit.Second.Along);
    }
    for (Overlap &overlap : into.Overlaps)
    {
      overlap.First  = {.Start = back(0, overlap.First.Start), .End = back(0, overlap.First.End)};
      overlap.Second = {.Start = back(1, overlap.Second.Start), .End = back(1, overlap.Second.End)};
    }
    return true;
  }

  void Entered(const Entering &entering, const Tolerance &sane, Intersections &into)
  {
    const Curves curves = {*entering.Pair.Shapes[0], *entering.Pair.Shapes[1]};
    const double length = entering.Length;
    if ((length > Coarse || length < Fine || Unending(curves))
        && (entering.Kept == nullptr || !Shrinkless(*entering.Kept))
        && Shrunk(curves, sane, into))
    {
      return;
    }
    if (Splinish(curves) && Loosened(curves, sane, into))
    {
      return;
    }
    if (!Directed(entering.Pair, {.Length = length, .Radians = sane.Radians}, into))
    {
      Thorough(curves, sane, into);
    }
  }

  void Deeply(const Curves &curves, const Tolerance &sane, Intersections &into)
  {
    into.Hits.clear();
    into.Overlaps.clear();
    into.Spent = false;
    /**
     *  The curves were within reach of each other with their arcs; with some made chords they are
     *  met within the length their numbers give, whatever their boxes now tell.
     */
    const double length = Reachable(curves, sane).value_or(Usable(curves, sane).Length);
    Entered({.Pair = {.Shapes = {&curves[0].get(), &curves[1].get()}}, .Length = length}, sane,
            into);
  }

  void Intersect(const Curves &curves, const Tolerance &tolerance, Intersections &into)
  {
    const Tolerance sane = Sane(tolerance);
    into.Hits.clear();
    into.Overlaps.clear();
    into.Spent                         = false;
    const std::optional<double> length = Reachable(curves, sane);
    const auto                  met = [](const void *, const Curves &stood, const Tolerance &held,
                                         Intersections &found) -> void
    {
      Intersect(stood, held, found);
    };
    if (!length)
    {
      return;
    }
    /**
     *  Only a polyline has arcs told by a bulge, which may be remote.
     */
    if (!std::holds_alternative<Polyline>(curves[0].get())
        && !std::holds_alternative<Polyline>(curves[1].get()))
    {
      Entered({.Pair = {.Shapes = {&curves[0].get(), &curves[1].get()}}, .Length = *length}, sane,
              into);
      return;
    }
    if (!Met::Slightly(curves, {.Sane = sane, .Length = *length, .Intersect = met, .Deep = Deeply},
                       into))
    {
      const Met::Vetting vetted;
      Entered({.Pair = {.Shapes = {&curves[0].get(), &curves[1].get()}}, .Length = *length}, sane,
              into);
    }
  }

  auto Intersect(const Curve &first, const Curve &second, const Tolerance &tolerance)
      -> Intersections
  {
    Intersections meeting;
    Intersect({first, second}, tolerance, meeting);
    return meeting;
  }

  auto SelfIntersect(const Curve &curve, const Tolerance &tolerance) -> Intersections
  {
    return Crossings::Self(curve, tolerance);
  }
} // namespace Geoxy
