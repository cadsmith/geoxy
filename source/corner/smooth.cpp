/***************************************************************************************************
 *  @file      source/corner/smooth.cpp
 *  @brief     Circles of a radius touching curves that are neither lines nor circles: against a
 *             line or a circle, the roots of one offset's distance from the other's, bracketed and
 *             halved; against another such curve, where their offsets' chords cross, polished by
 *             Newton's method; and the corners where such curves meet.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "corner/smooth.hpp"

#include "corner/circle.hpp"
#include "corner/forms.hpp"
#include "corner/guide.hpp"
#include "real.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curvature.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/derivative.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/extents.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <optional>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy::Corner
{
  namespace
  {
    /**
     *  How many chords lay a piece of a smooth curve's offset: a Bézier curve, a span of a spline,
     *  a quarter of an ellipse.
     */
    constexpr std::size_t Chords = 24;

    /**
     *  How many times a root is halved, or a Newton step taken, before it is given up on.
     */
    constexpr int Tries = 200;

    /**
     *  A point of a curve's offset, and the curve's parameter there.
     */
    struct Sample
    {
      double Along = 0.0; /**< The parameter.      */
      Point  At    = {};  /**< The offset's point. */
    };

    /**
     *  A smooth curve's point, its way and its unit normal at a parameter.
     */
    struct Frame
    {
      Point  At     = {};    /**< The point.                                   */
      Vector Way    = {};    /**< The derivative.                              */
      Vector Normal = {};    /**< The unit normal, a quarter turn to its left. */
      bool   Usable = false; /**< Whether the derivative is not nothing.       */
    };

    /**
     *  A smooth curve's frame at a parameter.
     *
     *  @param  curve The curve.
     *  @param  along The parameter.
     *  @return       Its frame.
     */
    auto Framed(const Curve &curve, double along) -> Frame
    {
      const Point  spot = At(curve, along);
      const Vector way  = Derivative(curve, along);
      const double size = Length(way);
      if (!(size > 0.0))
      {
        return {.At = spot, .Way = way};
      }
      return {.At = spot, .Way = way, .Normal = Perpendicular(way) / size, .Usable = true};
    }

    /**
     *  The parameters a smooth curve's offset is laid at: each piece's even steps, an ellipse's
     *  over a whole turn.
     *
     *  @param  guide The curve's guide.
     *  @return       The parameters, in order.
     */
    auto Steps(const Guide &guide) -> std::vector<double>
    {
      std::vector<double> cuts;
      if (guide.Period > 0.0)
      {
        cuts = {0.0, guide.Period};
      }
      else if (const auto *spline = std::get_if<Spline>(&guide.Shape.get()))
      {
        const Range domain = Domain(*spline);
        for (const double knot : spline->Knots)
        {
          cuts.push_back(std::clamp(knot, domain.Start, domain.End));
        }
        Ordered(cuts);
      }
      else
      {
        cuts = {0.0, 1.0};
      }
      const std::size_t   pieces = guide.Period > 0.0 ? 4 : 1;
      std::vector<double> alongs;
      for (std::size_t index = 0; index + 1 < cuts.size(); ++index)
      {
        const double      start = cuts[index];
        const double      span  = cuts[index + 1] - start;
        const std::size_t count = Chords * pieces;
        for (std::size_t step = 0; step < count; ++step)
        {
          alongs.push_back(start + (span * static_cast<double>(step) / static_cast<double>(count)));
        }
      }
      alongs.push_back(cuts.back());
      return alongs;
    }

    /**
     *  The circle of the radius about a point of a smooth curve's offset, as a found circle, the
     *  other curve's touch worked out from it.
     *
     *  @param  frames The smooth curve's frame where it is touched, and its parameter.
     *  @param  other  The other's guide, straight or round.
     *  @param  ask    The radius, the side of the smooth curve and of the other.
     *  @return        The circle, the smooth curve's touch first.
     */
    auto Seating(const std::pair<Frame, double> &frames, const Guide &other, const Ask &ask)
        -> Found
    {
      const Frame &frame  = frames.first;
      const Point  center = frame.At + (frame.Normal * (ask.Radius * ask.First));
      Point        touch  = {};
      if (other.Kind == Lies::Straight)
      {
        const Vector normal = Perpendicular(other.Way) / Length(other.Way);
        touch               = center - (normal * (ask.Radius * ask.Second));
      }
      else
      {
        const double reach = other.Bow.Radius + (ask.Radius * ask.Second);
        touch = other.Bow.Center + ((center - other.Bow.Center) * (other.Bow.Radius / reach));
      }
      return {.Touch  = {.Center = center, .Touches = {frame.At, touch}},
              .Alongs = {frames.second, Along(other, touch)},
              .Ask    = ask,
              .Exact  = false};
    }

    /**
     *  How far a point of a smooth curve's offset is from the other curve's offset: from a line's,
     *  signed, or from a circle's, less its radius.
     *
     *  @param  point The point.
     *  @param  other The other's guide, straight or round.
     *  @param  ask   The radius and the side of the other.
     *  @return       The gap: nothing where the circle about the point touches the other.
     */
    auto Gap(Point point, const Guide &other, const Ask &ask) -> double
    {
      if (other.Kind == Lies::Straight)
      {
        return (Cross(other.Way, point - other.Origin) / Length(other.Way))
             - (ask.Radius * ask.Second);
      }
      return Distance(point, other.Bow.Center)
           - std::abs(other.Bow.Radius + (ask.Radius * ask.Second));
    }

    /**
     *  The roots of a smooth curve's offset's gap from a line's or a circle's: each change of sign
     *  between two steps halved down to the parameter's last bit.
     *
     *  @param  smooth The smooth curve's guide.
     *  @param  other  The other's, straight or round.
     *  @param  ask    The radius, the side of the smooth curve and of the other.
     *  @return        The circles found, the smooth curve's touch first.
     */
    auto Rooted(const Guide &smooth, const Guide &other, const Ask &ask) -> std::vector<Found>
    {
      const Curve &curve = smooth.Shape.get();
      const auto   gap   = [&curve, &other, &ask](double along) -> double
      {
        const Frame frame = Framed(curve, along);
        return Gap(frame.At + (frame.Normal * (ask.Radius * ask.First)), other, ask);
      };
      std::vector<Found>        found;
      const std::vector<double> alongs = Steps(smooth);
      double                    last   = gap(alongs.front());
      for (std::size_t index = 1; index < alongs.size(); ++index)
      {
        const double next = gap(alongs[index]);
        if ((last <= 0.0) != (next <= 0.0))
        {
          const bool below  = last <= 0.0;
          double     low    = alongs[index - 1];
          double     high   = alongs[index];
          double     middle = std::midpoint(low, high);
          while (middle > low && middle < high)
          {
            ((gap(middle) <= 0.0) == below ? low : high) = middle;
            middle                                       = std::midpoint(low, high);
          }
          found.push_back(Seating({Framed(curve, low), low}, other, ask));
        }
        last = next;
      }
      return found;
    }

    /**
     *  A smooth curve's offset laid as points.
     *
     *  @param  guide  The curve's guide.
     *  @param  offset The signed distance, positive to the left.
     *  @return        The points, with their parameters.
     */
    auto Laid(const Guide &guide, double offset) -> std::vector<Sample>
    {
      std::vector<Sample> samples;
      for (const double along : Steps(guide))
      {
        const Frame frame = Framed(guide.Shape.get(), along);
        if (frame.Usable)
        {
          samples.push_back({.Along = along, .At = frame.At + (frame.Normal * offset)});
        }
      }
      return samples;
    }

    /**
     *  Where two chords cross, as the fraction along each.
     *
     *  @param  first  One chord's ends.
     *  @param  second The other's.
     *  @return        The fractions, each in [0, 1]; none where they do not cross.
     */
    auto Crossed(const std::array<Point, 2> &first, const std::array<Point, 2> &second)
        -> std::optional<std::array<double, 2>>
    {
      const Vector one    = first[1] - first[0];
      const Vector other  = second[1] - second[0];
      const Vector apart  = second[0] - first[0];
      const double across = Cross(one, other);
      if (Zero(across))
      {
        return std::nullopt;
      }
      const double along = Cross(apart, other) / across;
      const double onto  = Cross(apart, one) / across;
      if (along < 0.0 || along > 1.0 || onto < 0.0 || onto > 1.0)
      {
        return std::nullopt;
      }
      return std::array{along, onto};
    }

    /**
     *  A point of a smooth curve's offset and its derivative.
     */
    struct Moving
    {
      Point  At  = {}; /**< The offset's point.      */
      Vector Way = {}; /**< The offset's derivative. */
    };

    /**
     *  A smooth curve's offset's point and derivative at a parameter: c + s r n, and c' (1 - s r k)
     *  for k the signed curvature.
     *
     *  @param  curve  The curve.
     *  @param  along  The parameter.
     *  @param  offset s r.
     *  @return        Them; where the curve has no way, its point, and no way.
     */
    auto Moved(const Curve &curve, double along, double offset) -> Moving
    {
      const Frame  frame = Framed(curve, along);
      const double bend  = Curvature(curve, along).value_or(0.0);
      return {.At = frame.At + (frame.Normal * offset), .Way = frame.Way * (1.0 - (offset * bend))};
    }

    /**
     *  The circles of the radius touching two smooth curves on given sides: where their offsets'
     *  chords cross, polished.
     *
     *  @param  guides The curves' guides.
     *  @param  ask    The radius and the sides.
     *  @return        The circles.
     */
    auto Meshed(const std::array<const Guide *, 2> &guides, const Ask &ask) -> std::vector<Found>
    {
      const std::array<double, 2> offsets = {ask.Radius * ask.First, ask.Radius * ask.Second};
      const std::vector<Sample>   one     = Laid(*guides[0], offsets[0]);
      const std::vector<Sample>   other   = Laid(*guides[1], offsets[1]);
      std::vector<Found>          found;
      for (std::size_t first = 0; first + 1 < one.size(); ++first)
      {
        for (std::size_t second = 0; second + 1 < other.size(); ++second)
        {
          const Sample                              &a        = one[first];
          const Sample                              &b        = one[first + 1];
          const Sample                              &c        = other[second];
          const Sample                              &d        = other[second + 1];
          const std::optional<std::array<double, 2>> crossing = Crossed({a.At, b.At}, {c.At, d.At});
          if (!crossing)
          {
            continue;
          }
          const std::array<double, 2> guess = {a.Along + ((b.Along - a.Along) * (*crossing)[0]),
                                               c.Along + ((d.Along - c.Along) * (*crossing)[1])};
          for (const std::array<double, 2> &alongs :
               Polished({.Guides = guides, .Distances = offsets}, guess)
                   .transform(
                       [](const std::array<double, 2> &value) -> std::vector<std::array<double, 2>>
                       { return {value}; })
                   .value_or(std::vector<std::array<double, 2>>{}))
          {
            const Frame near = Framed(guides[0]->Shape.get(), alongs[0]);
            const Frame far  = Framed(guides[1]->Shape.get(), alongs[1]);
            found.push_back({.Touch  = {.Center  = near.At + (near.Normal * offsets[0]),
                                        .Touches = {near.At, far.At}},
                             .Alongs = alongs,
                             .Ask    = ask,
                             .Exact  = false});
          }
        }
      }
      return found;
    }

    /**
     *  Where two curves meet, at least one smooth.
     *
     *  @param  first  One's guide.
     *  @param  second The other's.
     *  @return        The corners.
     */
    auto Corners(const Guide &first, const Guide &second) -> std::vector<Found>
    {
      const Curve one   = Carried(first);
      const Curve other = Carried(second);
      Box         box   = {};
      for (const Curve *curve : {&one, &other})
      {
        const Box extents = Extents(*curve);
        if (std::isfinite((extents.Max.X - extents.Min.X) + (extents.Max.Y - extents.Min.Y)))
        {
          box = Joined(box, extents);
        }
      }
      const Intersections meeting = Intersect(one, other, Tolerance::Of(box));
      std::vector<Found>  found;
      found.reserve(meeting.Hits.size());
      for (const Hit &hit : meeting.Hits)
      {
        found.push_back({.Touch  = {.Center = hit.At, .Touches = {hit.At, hit.At}},
                         .Alongs = {Cornered(first, hit.First), Cornered(second, hit.Second)},
                         .Ask    = {},
                         .Exact  = false});
      }
      return found;
    }
  } // namespace

  /**
   *  A curve as Intersect meets it for a corner: a line's or a circle's carrier whole, an ellipse
   *  whole, a spline or a Bézier curve as it is.
   *
   *  @param  guide The curve's guide.
   *  @return       The curve to meet.
   */
  auto Carried(const Guide &guide) -> Curve
  {
    if (guide.Kind == Lies::Straight)
    {
      return Xline{.Through = guide.Origin, .Way = guide.Way};
    }
    if (guide.Kind == Lies::Round)
    {
      return Geoxy::Circle(guide.Bow.Center, guide.Bow.Radius);
    }
    if (const auto *ellipse = std::get_if<Ellipse>(&guide.Shape.get()))
    {
      Ellipse whole = *ellipse;
      whole.Sweep   = ellipse->Sweep.InRadians() > 0.0 ? Turn : -Turn;
      return whole;
    }
    return guide.Shape.get();
  }

  /**
   *  A corner's parameter on a curve, from where Intersect found it on what it met.
   *
   *  @param  guide The curve's guide.
   *  @param  place Where on what was met.
   *  @return       The parameter on the curve: an ellipse's scaled from its whole turn.
   */
  auto Cornered(const Guide &guide, const Place &place) -> double
  {
    if (guide.Kind != Lies::Smooth)
    {
      return Along(guide, place.At);
    }
    return guide.Period > 0.0 ? place.Along * guide.Period : place.Along;
  }

  auto Polished(const Offsets &offsets, std::array<double, 2> guess)
      -> std::optional<std::array<double, 2>>
  {
    const Point  from  = At(offsets.Guides[0]->Shape.get(), guess[0]);
    const double scale = std::max({1.0, std::abs(from.X), std::abs(from.Y)});
    for (int step = 0; step < Tries; ++step)
    {
      const Moving one    = Moved(offsets.Guides[0]->Shape.get(), guess[0], offsets.Distances[0]);
      const Moving other  = Moved(offsets.Guides[1]->Shape.get(), guess[1], offsets.Distances[1]);
      const Vector miss   = one.At - other.At;
      const double turned = Cross(one.Way, other.Way);
      const double across =
          std::copysign(std::max(std::abs(turned), std::numeric_limits<double>::min()), turned);
      const double first  = Cross(miss, other.Way) / across;
      const double second = Cross(miss, one.Way) / across;
      const double size   = std::max({1.0, std::abs(guess[0]), std::abs(guess[1])});
      guess               = {guess[0] - first, guess[1] - second};
      if (std::abs(first) + std::abs(second) <= 0x1p-50 * size && Length(miss) <= 0x1p-30 * scale)
      {
        return guess;
      }
    }
    return std::nullopt;
  }

  auto Smoothly(const Guide &first, const Guide &second, double radius) -> std::vector<Found>
  {
    if (!(radius > 0.0))
    {
      return Corners(first, second);
    }
    std::vector<Found> found;
    const std::array   signs = {1, -1};
    for (const int near : signs)
    {
      for (const int far : signs)
      {
        const Ask  ask = {.Radius = radius, .First = near, .Second = far};
        const bool hollowed =
            (first.Kind == Lies::Round && near < 0 && Exactly(first.Bow.Radius, radius))
            || (second.Kind == Lies::Round && far < 0 && Exactly(second.Bow.Radius, radius));
        if (hollowed)
        {
          continue;
        }
        if (first.Kind == Lies::Smooth && second.Kind == Lies::Smooth)
        {
          const std::vector<Found> meshed = Meshed({&first, &second}, ask);
          found.insert(found.end(), meshed.begin(), meshed.end());
        }
        else if (first.Kind == Lies::Smooth)
        {
          const std::vector<Found> rooted = Rooted(first, second, ask);
          found.insert(found.end(), rooted.begin(), rooted.end());
        }
        else
        {
          const Ask swapped = {.Radius = radius, .First = far, .Second = near};
          for (Found circle : Rooted(second, first, swapped))
          {
            std::swap(circle.Touch.Touches[0], circle.Touch.Touches[1]);
            std::swap(circle.Alongs[0], circle.Alongs[1]);
            circle.Ask = ask;
            found.push_back(circle);
          }
        }
      }
    }
    return found;
  }
} // namespace Geoxy::Corner
