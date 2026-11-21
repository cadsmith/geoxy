/***************************************************************************************************
 *  @file      source/near/feet.cpp
 *  @brief     The feet of the perpendiculars from a point to a curve: straight runs and arcs in
 *             closed form, every other curve by the kernel's feet on its pieces, listed in order.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "near/feet.hpp"

#include "basis.hpp"
#include "bow.hpp"
#include "elliptic.hpp"
#include "kernel/piece.hpp"
#include "lean.hpp"
#include "real.hpp"
#include "scratch.hpp"
#include "secular.hpp"
#include "seek.hpp"
#include "straight.hpp"
#include "twofold.hpp"
#include "walk.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <optional>
#include <span>
#include <variant>
#include <vector>

namespace Geoxy::Near
{
  namespace
  {
    /**
     *  Whether ways are each from 2^-450 to 2^450 long in their larger coordinate: their products
     *  and squares stay among the normal doubles, and none is nothing.
     *
     *  @param  ways The ways, exact.
     *  @return      Whether they are.
     */
    auto Temperate(std::initializer_list<Lean::Way> ways) -> bool
    {
      return std::ranges::all_of(ways,
                                 [](const Lean::Way &way) -> bool
                                 {
                                   const double most =
                                       std::max(std::abs(way[0].High), std::abs(way[1].High));
                                   return most >= 0x1p-450 && most < 0x1p450;
                                 });
    }

    /**
     *  A coordinate moved by a double-double amount, rounded once.
     *
     *  @param  start The coordinate.
     *  @param  step  The amount, its low part within a step of a double of its high one.
     *  @return       The double nearest start + step but a hair from halfway.
     */
    auto Nudged(double start, Twofold::Double step) -> double
    {
      const Twofold::Two sum = Twofold::TwoSum(start, step.High);
      return sum.High + (sum.Low + step.Low);
    }

    /**
     *  A double-double brought back to its low part within a step of a double of its high one.
     *
     *  @param  value The number.
     *  @return       It, normalized.
     */
    auto Normed(Twofold::Double value) -> Twofold::Double
    {
      const Twofold::Two out = Twofold::FastTwoSum(value.High, value.Low);
      return {.High = out.High, .Low = out.Low};
    }

    /**
     *  A polyline's foot of a point at a parameter, as Trod works it out: at a vertex, the vertex
     *  as At gives it; on a bulge, the side its circle's point there lies on from the center, as At
     *  tells it, which for a foot is clearly one or the other.
     *
     *  @param  polyline The polyline, with a segment.
     *  @param  along    The foot's parameter.
     *  @param  point    The point.
     *  @return          The foot.
     */
    auto Planted(const Polyline &polyline, double along, Point point) -> Point
    {
      const Step step = Locate(polyline, along);
      if (Exactly(step.Fraction, 0.0) || Exactly(step.Fraction, 1.0))
      {
        return At(polyline, along);
      }
      double side = 1.0;
      if (!Straightness(step))
      {
        const Hoop hoop = Hooped(step.Start, step.End, step.Bulge);
        side = Dot(At(polyline, along) - hoop.Center, point - hoop.Center) < 0.0 ? -1.0 : 1.0;
      }
      return Trod(step, point, side);
    }

    /**
     *  The parameters of an arc's feet: where it passes the way from its center to the point, and
     *  the way back.
     *
     *  @param  arc   The arc.
     *  @param  point The point.
     *  @return       The parameters, unordered; none for the center.
     */
    auto Around(const Arc &arc, Point point) -> List<double>
    {
      List<double> alongs;
      if (point == arc.Center)
      {
        return alongs;
      }
      const double  toward  = Heading(point - arc.Center).InRadians();
      const Ellipse ellipse = Elliptic(arc);
      for (const double angle : {toward, toward + Pi})
      {
        const List<double> passes = Passes(ellipse, angle);
        alongs.insert(alongs.end(), passes.begin(), passes.end());
      }
      return alongs;
    }

    /**
     *  The parameters of a polyline segment's feet.
     *
     *  @param  segment The segment.
     *  @param  point   The point.
     *  @return         The segment's parameters, from 0 to 1, unordered.
     */
    auto Around(const Segment &segment, Point point) -> List<double>
    {
      if (const Line *line = std::get_if<Line>(&segment))
      {
        return Foot(Running(*line), point);
      }
      return Around(std::get<Arc>(segment), point);
    }

    /**
     *  What a polyline's segment has of the feet: its own, and which way it runs from the point at
     *  its ends, for the turns of the distance at its vertices.
     *
     *  @param  segment The segment.
     *  @param  point   The point.
     *  @return         The feet's parameters, and the segment's run away from the point at its
     *                  ends.
     */
    auto Solved(const Segment &segment, Point point) -> Kernel::Solution
    {
      Kernel::Solution solution = Leaned(segment, point, [](Vector away, Vector tangent) -> double
                                         { return Dot(away, tangent); });
      solution.Roots            = Around(segment, point);
      return solution;
    }

    /**
     *  The feet of a piece.
     *
     *  @param  point The point.
     *  @return       What solves for the feet on a piece.
     */
    auto Footed(Point point)
    {
      return [point](const Kernel::Piece &piece) -> Kernel::Solution
      {
        return Kernel::Perpendiculars(piece, point, Kernel::Budget);
      };
    }
  } // namespace

  void Lineup(std::span<Place> places)
  {
    for (std::size_t index = 1; index < places.size(); ++index)
    {
      const Place moving = places[index];
      std::size_t into   = index;
      while (into > 0 && moving.Along < places[into - 1].Along)
      {
        places[into] = places[into - 1];
        --into;
      }
      places[into] = moving;
    }
  }

  auto Dropped(Point center, double radius, Point point) -> Point
  {
    const Lean::Way away = Lean::Apart(point, center);
    if (!Temperate({away}))
    {
      const Vector plain = point - center;
      return center + (plain * (radius / Length(plain)));
    }
    const Lean::Sized     square = Lean::Dot(away, away);
    const Lean::Rooted    length = Lean::Root(square.Value);
    const Twofold::Double ratio =
        Lean::Over({.High = radius, .Low = 0.0}, length.Value, length.Inverse);
    return {.X = Lean::Moved(center.X, away[0], ratio).Value,
            .Y = Lean::Moved(center.Y, away[1], ratio).Value};
  }

  auto Trod(const Step &step, Point point, double side) -> Point
  {
    const auto plain = [&step] -> Point
    {
      if (Straightness(step))
      {
        return Lerp(step.Start, step.End, step.Fraction);
      }
      const Bight bow = Bowed(step.Start, step.End, step.Bulge);
      return Along(bow, bow.Sweep * step.Fraction);
    };
    const Lean::Way chord = Lean::Apart(step.End, step.Start);
    const Lean::Way away  = Lean::Apart(point, step.Start);
    if (!Temperate({chord, away}))
    {
      return plain();
    }
    if (Straightness(step))
    {
      const Twofold::Double square = Lean::Dot(chord, chord).Value;
      const Twofold::Double share =
          Normed(Lean::Over(Lean::Dot(away, chord).Value, square, 1.0 / square.High));
      return {.X = Lean::Moved(step.Start.X, chord[0], share).Value,
              .Y = Lean::Moved(step.Start.Y, chord[1], share).Value};
    }
    const double          bulge   = step.Bulge;
    const Twofold::Double one     = {.High = 1.0, .Low = 0.0};
    const Twofold::Two    product = Twofold::TwoProduct(bulge, bulge);
    const Twofold::Double shrunk  = {.High = product.High, .Low = product.Low};
    const Twofold::Double offset =
        Normed(Lean::Over(Lean::Plus(one, Lean::Minus(shrunk)), {.High = 4.0 * bulge, .Low = 0.0},
                          1.0 / (4.0 * bulge)));
    const std::array<Twofold::Double, 2> turned = {Lean::Minus(chord[1]), chord[0]};
    Lean::Way                            center = {};
    Lean::Way                            moved  = {};
    for (std::size_t axis = 0; axis < 2; ++axis)
    {
      center.at(axis) =
          Lean::Plus({.High = chord.at(axis).High / 2.0, .Low = chord.at(axis).Low / 2.0},
                     Lean::Times(turned.at(axis), offset));
      moved.at(axis) = Lean::Plus(away.at(axis), Lean::Minus(center.at(axis)));
    }
    if (!Temperate({moved}))
    {
      return plain();
    }
    const Twofold::Double reach = Lean::Dot(moved, moved).Value;
    const Twofold::Double share =
        Normed(Lean::Over(Lean::Dot(center, center).Value, reach, 1.0 / reach.High));
    const Twofold::Double scale = Lean::Root(share).Value;
    const Twofold::Double ratio = {.High = side * scale.High, .Low = side * scale.Low};
    return {.X = Nudged(step.Start.X, Lean::Plus(center[0], Lean::Times(moved[0], ratio))),
            .Y = Nudged(step.Start.Y, Lean::Plus(center[1], Lean::Times(moved[1], ratio)))};
  }

  auto Listed(const Line &line, Point point) -> List<Place>
  {
    return Placed(line, Foot(Running(line), point));
  }

  auto Listed(const Ray &ray, Point point) -> List<Place>
  {
    return Placed(ray, Foot(Running(ray), point));
  }

  auto Listed(const Xline &xline, Point point) -> List<Place>
  {
    return Placed(xline, Foot(Running(xline), point));
  }

  auto Listed(const Arc &arc, Point point) -> List<Place>
  {
    List<Place> places;
    if (point == arc.Center)
    {
      return places;
    }
    const double  toward  = Heading(point - arc.Center).InRadians();
    const Ellipse ellipse = Elliptic(arc);
    for (const double side : {1.0, -1.0})
    {
      const Point foot = Dropped(arc.Center, side * arc.Radius, point);
      for (const double along : Passes(ellipse, side > 0.0 ? toward : toward + Pi))
      {
        places.push_back({.At = foot, .Along = along});
      }
    }
    std::ranges::stable_sort(places, {}, &Place::Along);
    const auto close = [](const Place &lhs, const Place &rhs) -> bool
    {
      const double scale = std::max({1.0, std::abs(lhs.Along), std::abs(rhs.Along)});
      return rhs.Along - lhs.Along <= 64.0 * std::numeric_limits<double>::epsilon() * scale;
    };
    const auto repeated = std::ranges::unique(places, close);
    places.erase(repeated.begin(), repeated.end());
    return places;
  }

  auto Listed(const Ellipse &ellipse, Point point) -> List<Place>
  {
    const Vector x = ellipse.Map.X;
    const Vector y = ellipse.Map.Y;
    if (Zero(Dot(x, y)) && Exactly(Dot(x, x), Dot(y, y)))
    {
      /**
       *  A circle, however turned or mirrored: its feet straight out from its center, in the unit
       *  circle's frame, as an arc's, never by the kernel's quarters.
       */
      const std::optional<Matrix> unit = Inverse(ellipse.Map);
      const Point                 seen = unit ? *unit * point : Point{};
      if (!unit || seen == Point{})
      {
        return {};
      }
      const double       toward = Heading(seen - Point{}).InRadians();
      List<double>       alongs = Passes(ellipse, toward);
      const List<double> back   = Passes(ellipse, toward + Pi);
      alongs.insert(alongs.end(), back.begin(), back.end());
      return Placed(ellipse, alongs);
    }
    if (const std::optional<Footholds> holds = Plumbs(ellipse, point))
    {
      List<double> alongs;
      for (const double angle : std::span{holds->Angles}.first(holds->Count))
      {
        const List<double> passes = Passes(ellipse, angle);
        alongs.insert(alongs.end(), passes.begin(), passes.end());
      }
      return Placed(ellipse, alongs);
    }
    return Placed(ellipse, Rounded(ellipse, Footed(point),
                                   Residual(point, [](Vector away, Vector tangent) -> double
                                            { return Dot(away, tangent); })));
  }

  auto Listed(const Polyline &polyline, Point point) -> List<Place>
  {
    List<double> alongs = Segmented(polyline, [point](const Segment &segment) -> Kernel::Solution
                                    { return Solved(segment, point); });
    Gathered(alongs);
    List<Place> places;
    places.reserve(alongs.size());
    for (const double along : alongs)
    {
      places.push_back({.At = Planted(polyline, along, point), .Along = along});
    }
    return places;
  }

  auto Listed(const Spline &spline, Point point) -> List<Place>
  {
    if (Rational(spline))
    {
      const auto feet = [point](const Kernel::Tread &tread) -> Kernel::Solution
      {
        return Kernel::Perpendiculars(tread, point, Kernel::Budget);
      };
      return Sought(spline, feet,
                    [&spline, &feet](std::size_t span) -> Kernel::Solution
                    {
                      return feet({.Shape = Kernel::Section(spline, span),
                                   .Speed = Kernel::Derived(spline, span)});
                    });
    }
    return Placed(spline, Kernel::Solved(spline, Footed(point)));
  }

  auto Listed(const Bezier::Quadratic &quadratic, Point point) -> List<Place>
  {
    return Placed(quadratic, Footed(point)(Kernel::Polynomial(quadratic)).Roots);
  }

  auto Listed(const Bezier::Cubic &cubic, Point point) -> List<Place>
  {
    return Placed(cubic, Footed(point)(Kernel::Polynomial(cubic)).Roots);
  }

  auto Least(Point point, std::span<const Place> places) -> Place
  {
    const Place *first = &places.front();
    for (const Place &place : places)
    {
      if (place.Along < first->Along)
      {
        first = &place;
      }
    }
    Place  best  = *first;
    double least = Distance(point, best.At);
    for (const Place &place : places)
    {
      const double away = Distance(point, place.At);
      if (away < least || (Exactly(away, least) && place.Along < best.Along))
      {
        best  = place;
        least = away;
      }
    }
    return best;
  }

  auto Everywhere(const Polyline &polyline, Point point) -> Place
  {
    List<Place>       places = Listed(polyline, point);
    const std::size_t count  = Count(polyline);
    for (std::size_t index = 0; index <= count; ++index)
    {
      const auto along = static_cast<double>(index);
      places.push_back({.At = At(polyline, along), .Along = along});
    }
    return Least(point, places);
  }
} // namespace Geoxy::Near
