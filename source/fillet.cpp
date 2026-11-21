/***************************************************************************************************
 *  @file      source/fillet.cpp
 *  @brief     FILLET of two curves: the corner nearest the picks, the circle of the radius in the
 *             angle between the parts kept, settled to the nearest doubles, and each curve cut to
 *             where it touches it.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/fillet.hpp"

#include "corner/choose.hpp"
#include "corner/circle.hpp"
#include "corner/forms.hpp"
#include "corner/guide.hpp"
#include "corner/rework.hpp"
#include "finite.hpp"
#include "sound.hpp"
#include "trig.hpp"
#include "twofold.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <expected>
#include <optional>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  Two curves cut to a choice's touches, and the arc between them.
     *
     *  @param  guides The curves' guides.
     *  @param  choice The circle chosen and the sides kept.
     *  @param  radius The radius; 0 for the corner, with no arc.
     *  @return        The fillet; TooLarge where a touch lies past the far end of a part kept.
     */
    auto Assembled(const Corner::Guides &guides, const Corner::Choice &choice, double radius)
        -> std::expected<Filleted, Misfit>
    {
      const Corner::Found          placed = Corner::Placed(guides, choice.Found, radius);
      const Corner::Touch         &touch  = placed.Touch;
      const std::array<double, 2> &along  = placed.Alongs;
      std::expected<Curve, Misfit> first =
          Corner::Cut(*guides[0], {.At = touch.Touches[0], .Along = along[0]}, choice.Keeps[0]);
      std::expected<Curve, Misfit> second =
          Corner::Cut(*guides[1], {.At = touch.Touches[1], .Along = along[1]}, choice.Keeps[1]);
      if (!first || !second)
      {
        return std::unexpected(Misfit::TooLarge);
      }
      Filleted filleted = {.First   = std::move(*first),
                           .Second  = std::move(*second),
                           .Round   = std::nullopt,
                           .Touches = touch.Touches};
      if (radius > 0.0)
      {
        filleted.Round =
            Arc{.Center = touch.Center,
                .Radius = radius,
                .Start  = Heading(touch.Touches[0] - touch.Center),
                .Sweep  = Radians(Corner::Swept(touch.Center, touch.Touches, choice.Ahead))};
      }
      return filleted;
    }

    /**
     *  Two parallel lines joined by the half circle between them, at the first's end nearest its
     *  pick, the second cut or run on to it: AutoCAD's rule.
     *
     *  @param  guides The lines' guides.
     *  @return        The fillet; Parallel for an endless line first, or lines on one line.
     */
    auto Semicircled(const Corner::Guides &guides) -> std::expected<Filleted, Misfit>
    {
      const Corner::Guide   &first  = *guides[0];
      const Corner::Guide   &second = *guides[1];
      const Twofold::Two     gapX   = Twofold::TwoSum(second.Origin.X, -first.Origin.X);
      const Twofold::Two     gapY   = Twofold::TwoSum(second.Origin.Y, -first.Origin.Y);
      const Corner::Straight across = {.Anchor = first.Origin,
                                       .High   = {.X = gapX.High, .Y = gapY.High},
                                       .Low    = {.X = gapX.Low, .Y = gapY.Low}};
      if (std::holds_alternative<Xline>(first.Shape.get())
          || Corner::Parallel(first.Holder.Line, across))
      {
        return std::unexpected(Misfit::Parallel);
      }
      std::optional<Point> far;
      if (const auto *line = std::get_if<Line>(&first.Shape.get()))
      {
        far = line->End;
      }
      if (const auto *polyline = std::get_if<Polyline>(&first.Shape.get()))
      {
        far = std::get<Line>(Nth(*polyline, first.Segment)).End;
      }
      const Point        pick   = first.Origin + (first.Way * first.Pick);
      const bool         ending = far && Distance(*far, pick) < Distance(first.Origin, pick);
      const Point        from   = ending ? *far : first.Origin;
      const Vector       arrive = ending ? first.Way : -first.Way;
      const Corner::Span span   = Corner::Bridged(second.Holder.Line, from);
      const Point        foot   = span.Touch.Touches[1];
      std::expected<Curve, Misfit> cut;
      if (const auto *line = std::get_if<Line>(&second.Shape.get()))
      {
        cut = Distance(line->Start, foot) <= Distance(line->End, foot)
                ? Line{.Start = foot, .End = line->End}
                : Line{.Start = line->Start, .End = foot};
      }
      else
      {
        cut = Corner::Cut(second, {.At = foot, .Along = Corner::Along(second, foot)},
                          -Corner::Sign(Dot(arrive, second.Way)));
      }
      if (!cut)
      {
        return std::unexpected(cut.error());
      }
      const bool ahead = Cross(from - span.Touch.Center, arrive) > 0.0;
      return Filleted{.First   = first.Shape.get(),
                      .Second  = std::move(*cut),
                      .Round   = Arc{.Center = span.Touch.Center,
                                     .Radius = span.Radius,
                                     .Start  = Heading(from - span.Touch.Center),
                                     .Sweep  = Radians(ahead ? Pi : -Pi)},
                      .Touches = span.Touch.Touches};
    }
  } // namespace

  auto Fillet(const Curve &first, const Curve &second, const Filleting &filleting)
      -> std::expected<Filleted, Misfit>
  {
    const double                radius = filleting.Radius;
    const std::array<Point, 2> &picks  = filleting.Picks;
    if (!std::isfinite(radius)
        || !Finite(picks[0])
        || !Finite(picks[1])
        || !Sound(first)
        || !Sound(second))
    {
      return std::unexpected(Misfit::NotFinite);
    }
    if (radius < 0.0)
    {
      return std::unexpected(Misfit::Degenerate);
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
      return Semicircled(guides);
    }
    const std::vector<Corner::Found>   corners = Corner::Circles(guides, 0.0);
    const std::optional<Corner::Found> corner  = Corner::Closest(guides, corners, picks);
    if (!(radius > 0.0))
    {
      if (!corner)
      {
        return std::unexpected(Misfit::Apart);
      }
      Corner::Choice choice = {.Found = *corner};
      for (std::size_t index = 0; index < 2; ++index)
      {
        const Corner::Side side       = Corner::Sided(*guides.at(index), corner->Alongs.at(index));
        choice.Found.Alongs.at(index) = side.Corner;
        choice.Keeps.at(index)        = side.Keep != 0 ? side.Keep : 1;
      }
      return Assembled(guides, choice, 0.0);
    }
    const std::vector<Corner::Found>    circles = Corner::Circles(guides, radius);
    const std::optional<Corner::Choice> choice =
        corner ? Corner::Nooked(guides, {.At = &*corner, .All = corners}, circles)
               : Corner::Bridging(guides, picks, circles);
    if (!choice)
    {
      return std::unexpected(corner ? Misfit::TooLarge : Misfit::Apart);
    }
    return Assembled(guides, *choice, radius);
  }

  auto Fillet(const Polyline &polyline, double radius, const Tolerance &tolerance)
      -> std::expected<Reshaped, Misfit>
  {
    if (!std::isfinite(radius) || !Sound(Curve{polyline}))
    {
      return std::unexpected(Misfit::NotFinite);
    }
    if (!(radius > 0.0))
    {
      return std::unexpected(Misfit::Degenerate);
    }
    const auto joiner = [radius](const Corner::Guides &guides,
                                 Point corner) -> std::expected<Corner::Joint, Misfit>
    {
      if (guides[0]->Kind == Corner::Lies::Straight
          && guides[1]->Kind == Corner::Lies::Straight
          && Corner::Parallel(guides[0]->Holder.Line, guides[1]->Holder.Line))
      {
        return std::unexpected(Misfit::Parallel);
      }
      const Corner::Found notch = {.Touch  = {.Center = corner, .Touches = {corner, corner}},
                                   .Alongs = {1.0, 0.0}};
      const std::optional<Corner::Choice> choice = Corner::Nooked(
          guides, {.At = &notch, .All = {&notch, 1}}, Corner::Circles(guides, radius));
      if (!choice)
      {
        return std::unexpected(Misfit::TooLarge);
      }
      const Corner::Found placed = Corner::Placed(guides, choice->Found, radius);
      const double sweep = Corner::Swept(placed.Touch.Center, placed.Touch.Touches, choice->Ahead);
      return Corner::Joint{
          .Ends = placed.Touch.Touches, .Alongs = placed.Alongs, .Bulge = Trig::Tan(sweep / 4.0)};
    };
    return Corner::Reworked(polyline, joiner, tolerance);
  }
} // namespace Geoxy
