/***************************************************************************************************
 *  @file      source/corner/guide.cpp
 *  @brief     Curves as a fillet reads them, and cut or run on to a point: lines by their ends
 *             moved, arcs, ellipses, splines and Bézier curves by their parts, polylines by their
 *             vertices.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "corner/guide.hpp"

#include "corner/circle.hpp"
#include "real.hpp"
#include "trig.hpp"
#include "twofold.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/fillet.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/nearest.hpp"
#include "geoxy/part.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <expected>
#include <functional>
#include <optional>
#include <variant>

namespace Geoxy::Corner
{
  namespace
  {
    /**
     *  How far past an end of its domain an arc's or a segment's parameter may come, rounding
     *  alone, and still be at that end: a touch at an arc's end has a parameter worked out from its
     *  heading, within a few steps of a double of the end's.
     */
    constexpr double Slack = 0x1p-40;

    /**
     *  A guide of a line from a point along a way held in two doubles.
     *
     *  @param  guide  The guide so far: its curve and its pick.
     *  @param  origin The line's point at 0.
     *  @param  ahead  Where its parameter is 1: the way is ahead less origin, held exactly.
     *  @return        The guide; Degenerate for a way of no length.
     */
    auto Lined(Guide guide, Point origin, Point ahead) -> std::expected<Guide, Misfit>
    {
      const Twofold::Two x = Twofold::TwoSum(ahead.X, -origin.X);
      const Twofold::Two y = Twofold::TwoSum(ahead.Y, -origin.Y);
      guide.Kind           = Lies::Straight;
      guide.Origin         = origin;
      guide.Way            = {.X = x.High, .Y = y.High};
      guide.Holder         = {.Linear = true,
                              .Line   = {.Anchor = origin,
                                         .High   = {.X = x.High, .Y = y.High},
                                         .Low    = {.X = x.Low, .Y = y.Low}}};
      if (Zero(x.High) && Zero(y.High))
      {
        return std::unexpected(Misfit::Degenerate);
      }
      return guide;
    }

    /**
     *  A guide of a line from a point along a way that is a double.
     *
     *  @param  guide  The guide so far.
     *  @param  origin The line's point at 0.
     *  @param  way    Its way.
     *  @return        The guide; Degenerate for a way of no length.
     */
    auto Lined(Guide guide, Point origin, Vector way) -> std::expected<Guide, Misfit>
    {
      guide.Kind   = Lies::Straight;
      guide.Origin = origin;
      guide.Way    = way;
      guide.Holder = {.Linear = true, .Line = {.Anchor = origin, .High = way, .Low = {}}};
      if (Zero(way.X) && Zero(way.Y))
      {
        return std::unexpected(Misfit::Degenerate);
      }
      return guide;
    }

    /**
     *  A guide of an arc.
     *
     *  @param  guide The guide so far.
     *  @param  arc   The arc.
     *  @return       The guide; Degenerate for an arc of no radius or no sweep.
     */
    auto Bowed(Guide guide, const Arc &arc) -> std::expected<Guide, Misfit>
    {
      const double sweep = std::abs(arc.Sweep.InRadians());
      if (!(arc.Radius > 0.0) || !(sweep > 0.0))
      {
        return std::unexpected(Misfit::Degenerate);
      }
      guide.Kind   = Lies::Round;
      guide.Bow    = arc;
      guide.Holder = {.Linear = false, .Circle = {.Center = arc.Center, .Radius = arc.Radius}};
      guide.Period = 2.0 * Pi / sweep;
      guide.Whole  = !(sweep < 2.0 * Pi);
      return guide;
    }

    /**
     *  The parameter of a point's foot on a line.
     *
     *  @param  guide The line's guide.
     *  @param  point The point.
     *  @return       Its foot's parameter.
     */
    auto Footed(const Guide &guide, Point point) -> double
    {
      return Dot(point - guide.Origin, guide.Way) / Dot(guide.Way, guide.Way);
    }

    /**
     *  A guide of a polyline, by its segment picked.
     *
     *  @param  guide    The guide so far.
     *  @param  polyline The polyline.
     *  @param  pick     Where it was picked.
     *  @return          The guide of the segment, its pick on it; Degenerate for no segment, or a
     *                   straight one of no length.
     */
    auto Picked(Guide guide, const Polyline &polyline, Point pick) -> std::expected<Guide, Misfit>
    {
      const std::size_t count = Count(polyline);
      if (count == 0)
      {
        return std::unexpected(Misfit::Degenerate);
      }
      const double      along = Nearest(polyline, pick).Along;
      const std::size_t index =
          std::min(static_cast<std::size_t>(std::max(std::floor(along), 0.0)), count - 1);
      guide.Segment                      = index;
      guide.Whole                        = polyline.Closed;
      const Segment                piece = Nth(polyline, index);
      std::expected<Guide, Misfit> found =
          std::holds_alternative<Line>(piece)
              ? Lined(guide, std::get<Line>(piece).Start, std::get<Line>(piece).End)
              : Bowed(guide, std::get<Arc>(piece));
      if (found)
      {
        found->Whole = polyline.Closed;
        found->Pick  = along - static_cast<double>(index);
      }
      return found;
    }

    /**
     *  Whether a point on a line lies past one of its points, the way it runs: told on the axis it
     *  runs along the most, where both are doubles and nothing rounds.
     *
     *  @param  point The point.
     *  @param  end   The line's point.
     *  @param  way   The way.
     *  @return       Whether it lies strictly past it.
     */
    auto Past(Point point, Point end, Vector way) -> bool
    {
      if (std::abs(way.X) >= std::abs(way.Y))
      {
        return way.X > 0.0 ? point.X > end.X : point.X < end.X;
      }
      return way.Y > 0.0 ? point.Y > end.Y : point.Y < end.Y;
    }

    /**
     *  The bulge of a part of a bulged segment.
     *
     *  @param  bulge    The segment's bulge.
     *  @param  fraction The part of its sweep kept, past 1 run on.
     *  @return          The part's bulge; none where it would sweep a whole turn or more.
     */
    auto Bulged(double bulge, double fraction) -> std::optional<double>
    {
      const double quarter = Trig::Atan(bulge) * fraction;
      if (!(std::abs(quarter) < Pi / 2.0))
      {
        return std::nullopt;
      }
      return Trig::Tan(quarter);
    }

    /**
     *  A polyline cut or run on along its picked segment to a point.
     *
     *  @param  guide Its guide, of an open polyline.
     *  @param  ends  The point, and its parameter on the segment.
     *  @param  keep  Which side is kept.
     *  @return       The polyline; TooLarge past the far end of the segment, or where a bulged
     *                segment would run on past a turn.
     */
    auto Docked(const Guide &guide, const Place &ends, int keep) -> std::expected<Curve, Misfit>
    {
      const auto       &polyline = std::get<Polyline>(guide.Shape.get());
      const std::size_t index    = guide.Segment;
      const double      bulge    = polyline.Vertices[index].Bulge;
      const double      along    = ends.Along;
      Polyline          cut;
      if (keep > 0)
      {
        const std::optional<double> rest = Bulged(bulge, 1.0 - along);
        if (along > 1.0 + Slack || !rest)
        {
          return std::unexpected(Misfit::TooLarge);
        }
        cut.Vertices.push_back({.Position = ends.At, .Bulge = *rest});
        cut.Vertices.insert(cut.Vertices.end(),
                            polyline.Vertices.begin() + static_cast<std::ptrdiff_t>(index + 1),
                            polyline.Vertices.end());
        return cut;
      }
      const std::optional<double> kept = Bulged(bulge, along);
      if (along < -Slack || !kept)
      {
        return std::unexpected(Misfit::TooLarge);
      }
      cut.Vertices.assign(polyline.Vertices.begin(),
                          polyline.Vertices.begin() + static_cast<std::ptrdiff_t>(index + 1));
      cut.Vertices.back().Bulge = *kept;
      cut.Vertices.push_back({.Position = ends.At, .Bulge = 0.0});
      return cut;
    }

    /**
     *  The stretch kept of a curve whose parameter repeats every turn, from a point to the far end
     *  of its domain or from its start to the point.
     *
     *  @param  guide The guide.
     *  @param  along The point's parameter, the one the part kept runs from or to.
     *  @param  keep  Which side is kept.
     *  @return       The stretch, an end within Slack of the domain's taken as at it; none past the
     *                far end, or past a turn.
     */
    auto Arcing(const Guide &guide, double along, int keep) -> std::optional<Range>
    {
      if (keep > 0)
      {
        if (along > 1.0 + Slack || !(along > 1.0 - guide.Period))
        {
          return std::nullopt;
        }
        return Range{.Start = std::min(along, 1.0), .End = 1.0};
      }
      if (along < -Slack || !(along < guide.Period))
      {
        return std::nullopt;
      }
      return Range{.Start = 0.0, .End = std::max(along, 0.0)};
    }
  } // namespace

  auto Guided(const Curve &curve, Point pick) -> std::expected<Guide, Misfit>
  {
    Guide guide = {.Shape = std::cref(curve)};
    if (const auto *line = std::get_if<Line>(&curve))
    {
      return Lined(guide, line->Start, line->End)
          .transform(
              [pick](Guide found) -> Guide
              {
                found.Pick = Footed(found, pick);
                return found;
              });
    }
    if (const auto *ray = std::get_if<Ray>(&curve))
    {
      return Lined(guide, ray->Start, ray->Way)
          .transform(
              [pick](Guide found) -> Guide
              {
                found.Pick = Footed(found, pick);
                return found;
              });
    }
    if (const auto *xline = std::get_if<Xline>(&curve))
    {
      return Lined(guide, xline->Through, xline->Way)
          .transform(
              [pick](Guide found) -> Guide
              {
                found.Pick = Footed(found, pick);
                return found;
              });
    }
    if (const auto *arc = std::get_if<Arc>(&curve))
    {
      return Bowed(guide, *arc)
          .transform(
              [arc, pick](Guide found) -> Guide
              {
                found.Pick = Nearest(*arc, pick).Along;
                return found;
              });
    }
    if (const auto *polyline = std::get_if<Polyline>(&curve))
    {
      return Picked(guide, *polyline, pick);
    }
    guide.Kind = Lies::Smooth;
    guide.Pick = Nearest(curve, pick).Along;
    if (const auto *ellipse = std::get_if<Ellipse>(&curve))
    {
      const double  sweep = std::abs(ellipse->Sweep.InRadians());
      const Matrix &map   = ellipse->Map;
      if (!(sweep > 0.0) || Zero(Cross(map.X, map.Y)))
      {
        return std::unexpected(Misfit::Degenerate);
      }
      guide.Period = 2.0 * Pi / sweep;
      guide.Whole  = !(sweep < 2.0 * Pi);
      return guide;
    }
    if (const auto *spline = std::get_if<Spline>(&curve))
    {
      if (!Valid(*spline))
      {
        return std::unexpected(Misfit::Degenerate);
      }
      const Range domain = Domain(*spline);
      guide.Whole        = At(*spline, domain.Start) == At(*spline, domain.End);
    }
    return guide;
  }

  auto Along(const Guide &guide, Point point) -> double
  {
    if (guide.Kind == Lies::Straight)
    {
      return Footed(guide, point);
    }
    const double turned = Between(Direction(guide.Bow.Start), point - guide.Bow.Center).InRadians();
    const double along  = turned / guide.Bow.Sweep.InRadians();
    return along < 0.0 ? along + guide.Period : along;
  }

  auto Within(const Guide &guide, double along) -> bool
  {
    if (guide.Whole)
    {
      return true;
    }
    if (guide.Period > 0.0)
    {
      return along <= 1.0 + Slack || along >= guide.Period - Slack;
    }
    const Range domain = std::holds_alternative<Polyline>(guide.Shape.get())
                           ? Range{.Start = 0.0, .End = 1.0}
                           : Domain(guide.Shape.get());
    return along >= domain.Start - Slack && along <= domain.End + Slack;
  }

  auto Near(const Guide &guide, double along, double target) -> double
  {
    if (!(guide.Period > 0.0))
    {
      return along;
    }
    return along + (std::round((target - along) / guide.Period) * guide.Period);
  }

  auto Ahead(const Guide &guide, double along, Point point) -> Vector
  {
    if (guide.Kind == Lies::Straight)
    {
      return guide.Way;
    }
    if (guide.Kind == Lies::Round)
    {
      return Perpendicular(point - guide.Bow.Center)
           * (guide.Bow.Sweep.InRadians() > 0.0 ? 1.0 : -1.0);
    }
    return Derivative(guide.Shape.get(), along);
  }

  auto Cut(const Guide &guide, const Place &place, int keep) -> std::expected<Curve, Misfit>
  {
    const Point  point = place.At;
    const double along = place.Along;
    const Curve &curve = guide.Shape.get();
    if (guide.Whole)
    {
      return curve;
    }
    if (const auto *line = std::get_if<Line>(&curve))
    {
      if (keep > 0)
      {
        if (Past(point, line->End, guide.Way))
        {
          return std::unexpected(Misfit::TooLarge);
        }
        return Line{.Start = point, .End = line->End};
      }
      if (Past(line->Start, point, guide.Way))
      {
        return std::unexpected(Misfit::TooLarge);
      }
      return Line{.Start = line->Start, .End = point};
    }
    if (const auto *ray = std::get_if<Ray>(&curve))
    {
      if (keep > 0)
      {
        return Ray{.Start = point, .Way = ray->Way};
      }
      if (Past(ray->Start, point, guide.Way))
      {
        return std::unexpected(Misfit::TooLarge);
      }
      return Line{.Start = ray->Start, .End = point};
    }
    if (const auto *xline = std::get_if<Xline>(&curve))
    {
      return Ray{.Start = point, .Way = keep > 0 ? xline->Way : -xline->Way};
    }
    if (std::holds_alternative<Polyline>(curve))
    {
      return Docked(guide, place, keep);
    }
    if (guide.Period > 0.0)
    {
      const std::optional<Range> stretch = Arcing(guide, along, keep);
      if (!stretch)
      {
        return std::unexpected(Misfit::TooLarge);
      }
      if (const auto *arc = std::get_if<Arc>(&curve))
      {
        return Part(*arc, *stretch);
      }
      return Part(std::get<Ellipse>(curve), *stretch);
    }
    const Range domain = Domain(curve);
    const Range kept   = keep > 0 ? Range{.Start = along, .End = domain.End}
                                  : Range{.Start = domain.Start, .End = along};
    return Part(curve, kept).value_or(curve);
  }
} // namespace Geoxy::Corner
