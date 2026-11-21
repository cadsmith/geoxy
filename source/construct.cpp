/***************************************************************************************************
 *  @file      source/construct.cpp
 *  @brief     Arcs, circles and lines built from what they pass through or touch.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/construct.hpp"

#include "real.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/line.hpp"
#include "geoxy/orient.hpp"
#include "geoxy/point.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  A sweep from one direction to another about a center, the way it must run.
     *
     *  @param  start The start, seen from the center.
     *  @param  end   The end, seen from the center.
     *  @param  ahead Whether it runs counterclockwise.
     *  @return       The angle from the start to the end, in (0, 2 pi] counterclockwise and in [-2
     *                pi, 0) clockwise, in radians.
     */
    auto Sweeping(Vector start, Vector end, bool ahead) -> Angle
    {
      const double between = Between(start, end).InRadians();
      if (ahead && !(between > 0.0))
      {
        return Radians(between + (2.0 * Pi));
      }
      if (!ahead && !(between < 0.0))
      {
        return Radians(between - (2.0 * Pi));
      }
      return Radians(between);
    }

    /**
     *  An arc about a center from a start to an end, the way it must run.
     *
     *  @param  center The center.
     *  @param  ends   Where it starts and where it ends, the end giving only its direction.
     *  @param  ahead  Whether it runs counterclockwise.
     *  @return        The arc, its radius the start's distance from the center.
     */
    auto Swung(Point center, std::array<Point, 2> ends, bool ahead) -> Arc
    {
      const Vector start = ends[0] - center;
      return {.Center = center,
              .Radius = Length(start),
              .Start  = Heading(start),
              .Sweep  = Sweeping(start, ends[1] - center, ahead)};
    }

    /**
     *  Whether a point's coordinates are finite.
     *
     *  @param  point The point.
     *  @return       Whether both are: their halves' sum is, which neither overflows nor hides an
     *                infinity or no number.
     */
    auto Representable(Point point) -> bool
    {
      return std::isfinite((point.X * 0.5) + (point.Y * 0.5));
    }

    /**
     *  The three unknowns of a circle tangent to curves: its center's coordinates and its radius.
     */
    using Unknowns = std::array<double, 3>;

    /**
     *  One of Apollonius's equations, linear in the unknowns: Row · (x, y, r) = Value.
     */
    struct Plane
    {
      Unknowns Row   = {};  /**< The factors of x, y and r. */
      double   Value = 0.0; /**< What they sum to.          */
    };

    /**
     *  One of Apollonius's equations for a circle: |c - Center|^2 = (Radius + Side r)^2.
     */
    struct Round
    {
      Point  Center = {};  /**< The circle's center, about the origin of the work. */
      double Radius = 0.0; /**< Its radius, not below nothing.                     */
      double Side   = 1.0; /**< 1 to touch it from without, -1 from within.        */
    };

    /**
     *  A curve a circle is built tangent to: a line, taken whole, or a circle.
     */
    struct Guide
    {
      bool   Straight = true; /**< Whether it is a line.                             */
      Point  Anchor   = {};   /**< A line's start, a circle's center.                */
      Vector Normal   = {};   /**< A line's unit normal, a quarter turn to its left. */
      double Radius   = 0.0;  /**< A circle's radius, not below nothing.             */
    };

    /**
     *  A line or an arc as a guide.
     *
     *  @param  segment The line or the arc.
     *  @return         The guide; none for a line of no length.
     */
    auto Guided(const Segment &segment) -> std::optional<Guide>
    {
      if (const Line *line = std::get_if<Line>(&segment))
      {
        return Unit(line->End - line->Start)
            .transform(
                [line](Vector way) -> Guide
                {
                  return {.Straight = true,
                          .Anchor   = line->Start,
                          .Normal   = Perpendicular(way),
                          .Radius   = 0.0};
                });
      }
      const Arc &arc = std::get<Arc>(segment);
      return Guide{
          .Straight = false, .Anchor = arc.Center, .Normal = {}, .Radius = std::abs(arc.Radius)};
    }

    /**
     *  A quotient, where it is one.
     *
     *  @param  numerator   The numerator.
     *  @param  denominator The denominator.
     *  @return             Their quotient; none over nothing.
     */
    auto Quotient(double numerator, double denominator) -> std::optional<double>
    {
      if (Zero(denominator))
      {
        return std::nullopt;
      }
      return numerator / denominator;
    }

    /**
     *  The determinant of three columns.
     *
     *  @param  first  The first column.
     *  @param  second The second.
     *  @param  third  The third.
     *  @return        first · (second × third).
     */
    auto Volume(const Unknowns &first, const Unknowns &second, const Unknowns &third) -> double
    {
      return (first[0] * ((second[1] * third[2]) - (second[2] * third[1])))
           - (first[1] * ((second[0] * third[2]) - (second[2] * third[0])))
           + (first[2] * ((second[0] * third[1]) - (second[1] * third[0])));
    }

    /**
     *  The solution of three linear equations, by Cramer's rule.
     *
     *  @param  lines The equations.
     *  @return       The unknowns; none where the equations do not fix them.
     */
    auto Cramered(std::span<const Plane, 3> lines) -> std::vector<Unknowns>
    {
      const Unknowns across = {lines[0].Row[0], lines[1].Row[0], lines[2].Row[0]};
      const Unknowns upward = {lines[0].Row[1], lines[1].Row[1], lines[2].Row[1]};
      const Unknowns round  = {lines[0].Row[2], lines[1].Row[2], lines[2].Row[2]};
      const Unknowns values = {lines[0].Value, lines[1].Value, lines[2].Value};
      const double   whole  = Volume(across, upward, round);
      if (Zero(whole))
      {
        return {};
      }
      return {{Volume(values, upward, round) / whole, Volume(across, values, round) / whole,
               Volume(across, upward, values) / whole}};
    }

    /**
     *  The solutions of two linear equations and a circle's: two unknowns follow the third, the one
     *  whose removal leaves the largest determinant, and the circle's equation is a quadratic in
     *  it.
     *
     *  @param  lines The linear equations, two.
     *  @param  base  The circle's.
     *  @return       The unknowns, none, one or two.
     */
    auto Cramered(std::span<const Plane, 2> lines, const Round &base) -> std::vector<Unknowns>
    {
      const std::span<const double> one         = lines[0].Row;
      const std::span<const double> two         = lines[1].Row;
      const auto                    determinant = [one, two](std::size_t free) -> double
      {
        const std::size_t first  = (free + 1) % 3;
        const std::size_t second = (free + 2) % 3;
        return (one[first] * two[second]) - (one[second] * two[first]);
      };
      const std::array<std::size_t, 3> frees = {0, 1, 2};
      const std::size_t                free =
          *std::ranges::max_element(frees, {}, [&determinant](std::size_t index) -> double
                                    { return std::abs(determinant(index)); });
      const double pivot = determinant(free);
      if (Zero(pivot))
      {
        return {};
      }
      /**
       *  Each unknown as start + slope t, t the free one.
       */
      const std::size_t       first  = (free + 1) % 3;
      const std::size_t       second = (free + 2) % 3;
      Unknowns                start  = {};
      Unknowns                slope  = {};
      const std::span<double> starts{start};
      const std::span<double> slopes{slope};
      const double            left  = lines[0].Value;
      const double            right = lines[1].Value;
      slopes[free]                  = 1.0;
      starts[first]                 = ((left * two[second]) - (right * one[second])) / pivot;
      slopes[first]  = ((two[free] * one[second]) - (one[free] * two[second])) / pivot;
      starts[second] = ((one[first] * right) - (two[first] * left)) / pivot;
      slopes[second] = ((two[first] * one[free]) - (one[first] * two[free])) / pivot;

      const double x      = start[0] - base.Center.X;
      const double y      = start[1] - base.Center.Y;
      const double reach  = base.Radius + (base.Side * start[2]);
      const double grows  = base.Side * slope[2];
      const double square = (slope[0] * slope[0]) + (slope[1] * slope[1]) - (grows * grows);
      const double linear = 2.0 * ((x * slope[0]) + (y * slope[1]) - (grows * reach));
      const double level  = (x * x) + (y * y) - (reach * reach);
      /**
       *  The roots of square t^2 + linear t + level, each as the quotient that does not cancel:
       *  half / square and level / half, half the larger in magnitude of the two sums. A root
       *  repeated is one, the first; where square is nothing, the first is none and the second the
       *  equation's one.
       */
      const double discriminant = (linear * linear) - (4.0 * square * level);
      if (!(discriminant >= 0.0))
      {
        return {};
      }
      const double half = -(linear + std::copysign(std::sqrt(discriminant), linear)) / 2.0;
      std::vector<std::optional<double>> roots = {Quotient(half, square)};
      if (discriminant > 0.0)
      {
        roots.push_back(Quotient(level, half));
      }
      std::vector<Unknowns> found;
      for (const std::optional<double> root : roots)
      {
        const double value = root.value_or(0.0);
        if (root)
        {
          found.push_back({start[0] + (slope[0] * value), start[1] + (slope[1] * value),
                           start[2] + (slope[2] * value)});
        }
      }
      return found;
    }

    /**
     *  The circles tangent to guides on given sides, about an origin near them.
     *
     *  @param  guides The guides, about the origin.
     *  @param  sides  Which side of each, 1 or -1.
     *  @param  radius The radius where it is given; none where it is to be found.
     *  @return        The centers and radii, about the origin.
     */
    auto Centers(std::span<const Guide> guides, std::span<const double> sides,
                 std::optional<double> radius) -> std::vector<Unknowns>
    {
      std::vector<Plane>   lines;
      std::optional<Round> base;
      for (std::size_t index = 0; index < guides.size(); ++index)
      {
        const Guide &guide = guides[index];
        if (guide.Straight)
        {
          lines.push_back({.Row   = {guide.Normal.X, guide.Normal.Y, -sides[index]},
                           .Value = Dot(guide.Normal, guide.Anchor - Point{})});
          continue;
        }
        const Round round = {.Center = guide.Anchor, .Radius = guide.Radius, .Side = sides[index]};
        if (!base)
        {
          base = round;
          continue;
        }
        /**
         *  Less the base circle's, a circle's equation is linear: -2 (C - B) · c - 2 (s R - t S) r
         *  = |B|^2 - |C|^2 + R^2 - S^2, for B, S, t the base's center, radius and side.
         */
        const Vector across = round.Center - base->Center;
        lines.push_back(
            {.Row   = {-2.0 * across.X, -2.0 * across.Y,
                       -2.0 * ((round.Side * round.Radius) - (base->Side * base->Radius))},
             .Value = Dot(base->Center - Point{}, base->Center - Point{})
                    - Dot(round.Center - Point{}, round.Center - Point{})
                    + ((round.Radius - base->Radius) * (round.Radius + base->Radius))});
      }
      if (radius)
      {
        lines.push_back({.Row = {0.0, 0.0, 1.0}, .Value = *radius});
      }
      if (base)
      {
        return Cramered(std::span<const Plane, 2>{lines.data(), 2}, *base);
      }
      return Cramered(std::span<const Plane, 3>{lines.data(), 3});
    }

    /**
     *  Where a circle touches a guide, and how it lies against it.
     *
     *  @param  guide  The guide.
     *  @param  circle The circle's center and radius.
     *  @param  side   Which side of the guide, 1 or -1.
     *  @return        The touch; none where the circle and a circle guide share their center.
     */
    auto Touchpoint(const Guide &guide, const Unknowns &circle, double side)
        -> std::optional<Grazing>
    {
      const Point center = {.X = circle[0], .Y = circle[1]};
      if (guide.Straight)
      {
        return Grazing{.At   = center - (guide.Normal * (side * circle[2])),
                       .Kind = side > 0.0 ? Qualifier::Left : Qualifier::Right};
      }
      const std::optional<Vector> away = Unit(center - guide.Anchor);
      if (!away)
      {
        return std::nullopt;
      }
      if (side > 0.0)
      {
        return Grazing{.At = guide.Anchor + (*away * guide.Radius), .Kind = Qualifier::Outside};
      }
      if (circle[2] < guide.Radius)
      {
        return Grazing{.At = guide.Anchor + (*away * guide.Radius), .Kind = Qualifier::Enclosed};
      }
      return Grazing{.At = guide.Anchor - (*away * guide.Radius), .Kind = Qualifier::Enclosing};
    }

    /**
     *  A circle tangent to up to three curves, with where it touches each.
     */
    struct Solution
    {
      Arc                    Circle  = {}; /**< The circle.                       */
      std::array<Grazing, 3> Touches = {}; /**< Where it touches each, the first. */
    };

    /**
     *  Curves as guides.
     *
     *  @param  segments The curves.
     *  @return          Their guides; none where a line has no length.
     */
    auto Guides(std::span<const Segment *const> segments) -> std::optional<std::vector<Guide>>
    {
      std::vector<Guide> guides;
      for (const Segment *segment : segments)
      {
        const std::optional<Guide> guide = Guided(*segment);
        if (!guide)
        {
          return std::nullopt;
        }
        guides.push_back(*guide);
      }
      return guides;
    }

    /**
     *  The sides of the guides a choice puts a circle on.
     *
     *  @param  choice The choice: its bit at an index set for the side -1 of that guide.
     *  @param  count  How many guides.
     *  @return        The sides, 1 or -1 each.
     */
    auto Signs(std::size_t choice, std::size_t count) -> std::vector<double>
    {
      std::vector<double> sides;
      sides.reserve(count);
      for (std::size_t index = 0; index < count; ++index)
      {
        sides.push_back(((choice >> index) & 1U) == 0 ? 1.0 : -1.0);
      }
      return sides;
    }

    /**
     *  Whether a choice of sides repeats another: within a circle of no radius, a point, is without
     *  it.
     *
     *  @param  guides The guides.
     *  @param  sides  The sides.
     *  @return        Whether it does.
     */
    auto Repeated(std::span<const Guide> guides, std::span<const double> sides) -> bool
    {
      for (std::size_t index = 0; index < guides.size(); ++index)
      {
        if (sides[index] < 0.0 && !guides[index].Straight && Zero(guides[index].Radius))
        {
          return true;
        }
      }
      return false;
    }

    /**
     *  A circle found, with where it touches each guide.
     *
     *  @param  guides The guides.
     *  @param  sides  The sides of them it lies on.
     *  @param  found  Its center and radius.
     *  @return        The circle; none where it shares its center with a circle it touches.
     */
    auto Positioned(std::span<const Guide> guides, std::span<const double> sides,
                    const Unknowns &found) -> std::optional<Solution>
    {
      Solution circle = {.Circle  = Geoxy::Circle({.X = found[0], .Y = found[1]}, found[2]),
                         .Touches = {}};
      const std::span<Grazing> touches = circle.Touches;
      for (std::size_t index = 0; index < guides.size(); ++index)
      {
        const std::optional<Grazing> touch = Touchpoint(guides[index], found, sides[index]);
        if (!touch)
        {
          return std::nullopt;
        }
        touches[index] = *touch;
      }
      return circle;
    }

    /**
     *  Every circle tangent to curves, on every choice of sides.
     *
     *  @param  segments The curves, two or three.
     *  @param  radius   The radius where it is given.
     *  @return          The circles, with their touches.
     */
    auto Apollonius(std::span<const Segment *const> segments, std::optional<double> radius)
        -> std::vector<Solution>
    {
      const std::optional<std::vector<Guide>> guides = Guides(segments);
      if (!guides)
      {
        return {};
      }
      /**
       *  The work is about the first guide's anchor, so that its sums keep the digits of the
       *  curves' sizes, not of their places.
       */
      const Vector       origin = guides->front().Anchor - Point{};
      std::vector<Guide> local  = *guides;
      for (Guide &guide : local)
      {
        guide.Anchor = guide.Anchor - origin;
      }
      std::vector<Solution> circles;
      for (std::size_t choice = 0; choice < (std::size_t{1} << guides->size()); ++choice)
      {
        const std::vector<double> sides = Signs(choice, guides->size());
        if (Repeated(*guides, sides))
        {
          continue;
        }
        for (Unknowns found : Centers(local, sides, radius))
        {
          found[0] += origin.X;
          found[1] += origin.Y;
          if (!(found[2] > 0.0) || !Representable({.X = found[0], .Y = found[1]}))
          {
            continue;
          }
          found[2]                             = radius.value_or(found[2]);
          const std::optional<Solution> circle = Positioned(*guides, sides, found);
          if (circle)
          {
            circles.push_back(*circle);
          }
        }
      }
      return circles;
    }

    /**
     *  Circles tangent to curves, with as many touches as there are curves.
     *
     *  @tparam CountT    How many curves.
     *  @param  solutions The circles.
     *  @return           The circles, each with its first touches.
     */
    template <std::size_t CountT>
    auto Narrowed(const std::vector<Solution> &solutions) -> std::vector<Tangential<CountT>>
    {
      std::vector<Tangential<CountT>> circles;
      for (const Solution &solution : solutions)
      {
        Tangential<CountT> circle = {.Circle = solution.Circle, .Touches = {}};
        std::ranges::copy(std::span{solution.Touches}.first(CountT), circle.Touches.begin());
        circles.push_back(circle);
      }
      return circles;
    }
  } // namespace

  auto Arc::Through(Point start, Point through, Point end) -> std::optional<Arc>
  {
    const Orientation turn   = Orient(start, through, end);
    const Vector      first  = through - start;
    const Vector      second = end - start;
    const double      twice  = 2.0 * Cross(first, second);
    if (turn == Orientation::Collinear || Zero(twice))
    {
      return std::nullopt;
    }
    const double near   = Dot(first, first);
    const double far    = Dot(second, second);
    const Point  center = start
                        + Vector{.X = ((second.Y * near) - (first.Y * far)) / twice,
                                 .Y = ((first.X * far) - (second.X * near)) / twice};
    if (!Representable(center))
    {
      return std::nullopt;
    }
    return Swung(center, {start, end}, turn == Orientation::Counterclockwise);
  }

  auto Arc::Around(Point start, Point center, Point end) -> std::optional<Arc>
  {
    if (start == center || end == center)
    {
      return std::nullopt;
    }
    return Swung(center, {start, end}, true);
  }

  auto Arc::Spanning(Point start, Point end, double radius) -> std::optional<Arc>
  {
    const Vector chord = end - start;
    const double span  = Length(chord);
    const double half  = span / 2.0;
    const double reach = std::abs(radius);
    if (!(span > 0.0) || !(reach >= half))
    {
      return std::nullopt;
    }
    const double rise = std::sqrt((reach - half) * (reach + half));
    const Point  center =
        Midpoint(start, end) + (Perpendicular(chord) * ((radius > 0.0 ? rise : -rise) / span));
    Arc arc    = Swung(center, {start, end}, true);
    arc.Radius = reach;
    return arc;
  }

  auto Arc::Leaving(Point start, Point end, Vector way) -> std::optional<Arc>
  {
    const std::optional<Vector> along = Unit(way);
    const Vector                chord = end - start;
    if (!along || Zero(Cross(*along, chord)))
    {
      return std::nullopt;
    }
    const double reach  = Dot(chord, chord) / (2.0 * Cross(*along, chord));
    const Point  center = start + (Perpendicular(*along) * reach);
    if (!Representable(center))
    {
      return std::nullopt;
    }
    return Swung(center, {start, end}, reach > 0.0);
  }

  auto Circle(Point first, Point second) -> Arc
  {
    return Circle(Midpoint(first, second), Distance(first, second) / 2.0);
  }

  auto Circle(Point first, Point second, Point third) -> std::optional<Arc>
  {
    return Arc::Through(first, second, third)
        .transform([](const Arc &arc) -> Arc { return Circle(arc.Center, arc.Radius); });
  }

  auto TangentCircles(const Segment &first, const Segment &second, double radius)
      -> std::vector<Tangential<2>>
  {
    if (!(radius > 0.0))
    {
      return {};
    }
    const std::array segments = {&first, &second};
    return Narrowed<2>(Apollonius(segments, radius));
  }

  auto TangentCircles(const Segment &first, const Segment &second, const Segment &third)
      -> std::vector<Tangential<3>>
  {
    const std::array segments = {&first, &second, &third};
    return Narrowed<3>(Apollonius(segments, std::nullopt));
  }

  auto TangentLines(Point from, const Arc &circle) -> std::vector<Line>
  {
    const Vector away     = from - circle.Center;
    const double distance = Length(away);
    const double radius   = std::abs(circle.Radius);
    if (!(distance >= radius) || !(distance > 0.0))
    {
      return {};
    }
    const double rise = std::sqrt((distance - radius) * (distance + radius));
    if (Zero(rise))
    {
      return {{.Start = from, .End = from}};
    }
    /**
     *  Each touch is radius / distance of the way toward the point and rise / distance across it,
     *  times the radius.
     */
    const Vector way    = away / distance;
    const Vector toward = way * (radius * radius / distance);
    const Vector across = Perpendicular(way) * (radius * rise / distance);
    return {{.Start = from, .End = circle.Center + (toward + across)},
            {.Start = from, .End = circle.Center + (toward - across)}};
  }

  auto TangentLines(const Arc &first, const Arc &second) -> std::vector<Line>
  {
    const Vector across   = second.Center - first.Center;
    const double distance = Length(across);
    if (!(distance > 0.0))
    {
      return {};
    }
    const Vector      way  = across / distance;
    const double      near = std::abs(first.Radius);
    const double      far  = std::abs(second.Radius);
    std::vector<Line> lines;
    for (const double side : {1.0, -1.0})
    {
      /**
       *  The line n · x + k = 0 has the first center at n · c1 + k = r1 and the second at side r2:
       *  n · (c2 - c1) = side r2 - r1, the rest of n turned from the centers' line.
       */
      const double along = ((side * far) - near) / distance;
      if (!(std::abs(along) <= 1.0))
      {
        continue;
      }
      const double turned = std::sqrt((1.0 - along) * (1.0 + along));
      for (const double sense : {1.0, -1.0})
      {
        const Vector normal = (way * along) + (Perpendicular(way) * (sense * turned));
        lines.push_back({.Start = first.Center - (normal * near),
                         .End   = second.Center - (normal * (side * far))});
      }
    }
    return lines;
  }
} // namespace Geoxy
