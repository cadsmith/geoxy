/***************************************************************************************************
 *  @file      source/chamfer.cpp
 *  @brief     CHAMFER of two curves and of a polyline's corners: the corner FILLET finds, each
 *             curve cut back from it, between lines to the doubles nearest the exact cuts, and the
 *             line joining the cuts.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/chamfer.hpp"

#include "corner/choose.hpp"
#include "corner/circle.hpp"
#include "corner/guide.hpp"
#include "corner/rework.hpp"
#include "corner/smooth.hpp"
#include "finite.hpp"
#include "sound.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/extents.hpp"
#include "geoxy/fillet.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <expected>
#include <optional>
#include <utility>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  Of the places where a curve meets what cuts it, the one on the side kept of the corner,
     *  nearest it.
     *
     *  @param  guide  The curve's guide.
     *  @param  hits   Where it was met, the curve's places second.
     *  @param  corner The corner's parameter, and the side kept.
     *  @return        The place; none where none lies on that side.
     */
    auto Retained(const Corner::Guide &guide, const std::vector<Hit> &hits,
                  const std::pair<double, int> &corner) -> std::optional<Place>
    {
      std::vector<std::pair<double, Place>> kept;
      for (const Hit &hit : hits)
      {
        const double along = Corner::Near(guide, Corner::Cornered(guide, hit.Second), corner.first);
        const double past  = (along - corner.first) * corner.second;
        if (past > 0.0)
        {
          kept.emplace_back(past, Place{.At = hit.At, .Along = along});
        }
      }
      if (kept.empty())
      {
        return std::nullopt;
      }
      return std::ranges::min_element(kept, {}, [](const std::pair<double, Place> &item) -> double
                                      { return item.first; })
          ->second;
    }

    /**
     *  A corner on a curve, and the side of it kept.
     */
    struct Leaning
    {
      Place At   = {}; /**< The corner, its parameter on the curve. */
      int   Keep = 1;  /**< The side kept.                          */
    };

    /**
     *  Where a curve is the distance from the corner, on the side kept, measured straight.
     *
     *  @param  guide    The curve's guide.
     *  @param  bearing  The corner and the side kept.
     *  @param  distance The distance, not below nothing.
     *  @return          The place: on a curve where Intersect rounds its crossing with the circle
     *                   of the distance about the corner; the corner for a distance of nothing;
     *                   none where the circle does not meet the side kept.
     */
    auto Measuring(const Corner::Guide &guide, const Leaning &bearing, double distance)
        -> std::optional<Place>
    {
      if (!(distance > 0.0))
      {
        return bearing.At;
      }
      const Curve         circle = Geoxy::Circle(bearing.At.At, distance);
      const Intersections meeting =
          Intersect(circle, Corner::Carried(guide), Tolerance::Of(Extents(circle)));
      return Retained(guide, meeting.Hits, {bearing.At.Along, bearing.Keep});
    }

    /**
     *  A corner of two curves, and the side of it each keeps.
     */
    struct Nick
    {
      Corner::Found      At    = {}; /**< The corner, its parameters those of the sides kept. */
      std::array<int, 2> Keeps = {}; /**< The side kept of each.                              */
    };

    /**
     *  Where two curves are cut back from their corner.
     *
     *  @param  guides The curves' guides.
     *  @param  notch  The corner and the sides kept.
     *  @param  bevel  The distances, or the distance and the angle.
     *  @return        The two cuts with their parameters; TooLarge where one is not on its part
     *                 kept.
     */
    auto Cutbacks(const Corner::Guides &guides, const Nick &notch, const Bevel &bevel)
        -> std::expected<std::array<Place, 2>, Misfit>
    {
      const Corner::Found      &corner = notch.At;
      const std::array<int, 2> &keeps  = notch.Keeps;
      const Corner::Guide      &first  = *guides[0];
      const Corner::Guide      &second = *guides[1];
      if (first.Kind == Corner::Lies::Straight && second.Kind == Corner::Lies::Straight)
      {
        Corner::Cutback cutback = {.Distances = {bevel.First, bevel.Second}, .Keeps = keeps};
        if (bevel.Slant)
        {
          const Vector direction = Direction(*bevel.Slant);
          cutback.Slant          = std::array{direction.Y, direction.X};
        }
        const std::optional<std::array<Point, 2>> cuts =
            Corner::Bevelled(first.Holder.Line, second.Holder.Line, cutback);
        if (!cuts)
        {
          return std::unexpected(Misfit::TooLarge);
        }
        return std::array{Place{.At = (*cuts)[0], .Along = Corner::Along(first, (*cuts)[0])},
                          Place{.At = (*cuts)[1], .Along = Corner::Along(second, (*cuts)[1])}};
      }
      const Place                start = {.At = corner.Touch.Touches[0], .Along = corner.Alongs[0]};
      const std::optional<Place> one =
          Measuring(first, {.At = start, .Keep = keeps[0]}, bevel.First);
      if (!one)
      {
        return std::unexpected(Misfit::TooLarge);
      }
      const Place          there = {.At = corner.Touch.Touches[1], .Along = corner.Alongs[1]};
      std::optional<Place> other;
      if (bevel.Slant && bevel.First > 0.0)
      {
        const Vector toward =
            Unit(Corner::Ahead(first, one->Along, one->At) * static_cast<double>(-keeps[0]))
                .value_or(Vector{});
        const Vector onto =
            Unit(Corner::Ahead(second, there.Along, there.At) * static_cast<double>(keeps[1]))
                .value_or(Vector{});
        const Vector        turn  = Direction(*bevel.Slant);
        const double        side  = Cross(toward, onto) < 0.0 ? -1.0 : 1.0;
        const Vector        way   = (toward * turn.X) + (Perpendicular(toward) * (side * turn.Y));
        const Curve         ray   = Ray{.Start = one->At, .Way = way};
        const Curve         along = Corner::Carried(second);
        const Intersections met   = Intersect(ray, along, Tolerance::Of(Extents(along)));
        other                     = Retained(second, met.Hits, {there.Along, keeps[1]});
      }
      else
      {
        other =
            Measuring(second, {.At = there, .Keep = keeps[1]}, bevel.Slant ? 0.0 : bevel.Second);
      }
      if (!other)
      {
        return std::unexpected(Misfit::TooLarge);
      }
      return std::array{*one, *other};
    }

    /**
     *  Whether a bevel is one CHAMFER takes.
     *
     *  @param  bevel The bevel.
     *  @return       NotFinite for a number that is not finite, Degenerate for a distance below
     *                nothing or an angle outside (0, 180) degrees; nothing where it is sound.
     */
    auto Misfitting(const Bevel &bevel) -> std::optional<Misfit>
    {
      const double slant = bevel.Slant ? bevel.Slant->InRadians() : 1.0;
      if (!std::isfinite(bevel.First) || !std::isfinite(bevel.Second) || !std::isfinite(slant))
      {
        return Misfit::NotFinite;
      }
      if (bevel.First < 0.0 || bevel.Second < 0.0 || !(slant > 0.0) || !(slant < Pi))
      {
        return Misfit::Degenerate;
      }
      return std::nullopt;
    }
  } // namespace

  auto Chamfer(const Curve &first, const Curve &second, const Chamfering &chamfering)
      -> std::expected<Bevelled, Misfit>
  {
    const Bevel                &bevel = chamfering.Cut;
    const std::array<Point, 2> &picks = chamfering.Picks;
    if (!Finite(picks[0]) || !Finite(picks[1]) || !Sound(first) || !Sound(second))
    {
      return std::unexpected(Misfit::NotFinite);
    }
    if (const std::optional<Misfit> unsound = Misfitting(bevel))
    {
      return std::unexpected(*unsound);
    }
    const std::expected<Corner::Guide, Misfit> one = Corner::Guided(first, picks[0]);
    const std::expected<Corner::Guide, Misfit> two = Corner::Guided(second, picks[1]);
    if (!one || !two)
    {
      return std::unexpected(Misfit::Degenerate);
    }
    const Corner::Guides guides = {&*one, &*two};
    if (one->Kind == Corner::Lies::Straight
        && two->Kind == Corner::Lies::Straight
        && Corner::Parallel(one->Holder.Line, two->Holder.Line))
    {
      return std::unexpected(Misfit::Parallel);
    }
    const std::optional<Corner::Found> corner =
        Corner::Closest(guides, Corner::Circles(guides, 0.0), picks);
    if (!corner)
    {
      return std::unexpected(Misfit::Apart);
    }
    Corner::Found      placed = Corner::Placed(guides, *corner, 0.0);
    std::array<int, 2> keeps  = {};
    for (std::size_t index = 0; index < 2; ++index)
    {
      const Corner::Side side = Corner::Sided(*guides.at(index), placed.Alongs.at(index));
      placed.Alongs.at(index) = side.Corner;
      keeps.at(index)         = side.Keep != 0 ? side.Keep : 1;
    }
    const std::expected<std::array<Place, 2>, Misfit> cuts =
        Cutbacks(guides, {.At = placed, .Keeps = keeps}, bevel);
    if (!cuts)
    {
      return std::unexpected(cuts.error());
    }
    const std::array<Place, 2>  &cut      = *cuts;
    std::expected<Curve, Misfit> cutFirst = Corner::Cut(*one, cut[0], keeps[0]);
    std::expected<Curve, Misfit> cutOther = Corner::Cut(*two, cut[1], keeps[1]);
    if (!cutFirst || !cutOther)
    {
      return std::unexpected(Misfit::TooLarge);
    }
    Bevelled bevelled = {.First = std::move(*cutFirst), .Second = std::move(*cutOther)};
    if (bevel.First > 0.0 || (!bevel.Slant && bevel.Second > 0.0))
    {
      bevelled.Edge = Line{.Start = cut[0].At, .End = cut[1].At};
    }
    return bevelled;
  }

  auto Chamfer(const Polyline &polyline, const Bevel &bevel, const Tolerance &tolerance)
      -> std::expected<Reshaped, Misfit>
  {
    if (!Sound(Curve{polyline}))
    {
      return std::unexpected(Misfit::NotFinite);
    }
    if (const std::optional<Misfit> unsound = Misfitting(bevel))
    {
      return std::unexpected(*unsound);
    }
    if (!(bevel.First > 0.0) && (bevel.Slant || !(bevel.Second > 0.0)))
    {
      return std::unexpected(Misfit::Degenerate);
    }
    const auto joiner = [&bevel](const Corner::Guides &guides,
                                 Point corner) -> std::expected<Corner::Joint, Misfit>
    {
      if (guides[0]->Kind == Corner::Lies::Straight
          && guides[1]->Kind == Corner::Lies::Straight
          && Corner::Parallel(guides[0]->Holder.Line, guides[1]->Holder.Line))
      {
        return std::unexpected(Misfit::Parallel);
      }
      const Corner::Found notch = {.Touch  = {.Center = corner, .Touches = {corner, corner}},
                                   .Alongs = {1.0, 0.0},
                                   .Exact  = false};
      return Cutbacks(guides, {.At = notch, .Keeps = {-1, 1}}, bevel)
          .transform(
              [](const std::array<Place, 2> &cuts) -> Corner::Joint
              {
                return {.Ends   = {cuts[0].At, cuts[1].At},
                        .Alongs = {cuts[0].Along, cuts[1].Along},
                        .Bulge  = 0.0};
              });
    };
    return Corner::Reworked(polyline, joiner, tolerance);
  }
} // namespace Geoxy
