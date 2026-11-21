/***************************************************************************************************
 *  @file      source/belt/carrier.cpp
 *  @brief     Where two pieces' lines or circles meet, in doubles and in two doubles a number, and
 *             how far a carrier's point is from a shape's piece.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "belt/carrier.hpp"

#include "belt/piece.hpp"
#include "twofold.hpp"

#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <array>
#include <cmath>

namespace Geoxy::Belt
{
  namespace
  {
    using Twofold::Add;
    using Twofold::Difference;
    using Twofold::Divide;
    using Twofold::Double;
    using Twofold::Lifted;
    using Twofold::Multiply;
    using Twofold::Root;

    /**
     *  How near touching, as a share of the radius squared, a line and a circle or two circles are
     *  read again in two doubles.
     */
    constexpr double Near = 1e-6;

    /**
     *  How far from touching, as a share of the radius squared, doubles take a line and a circle or
     *  two circles to touch.
     */
    constexpr double Slack = 1e-13;

    /**
     *  How near touching, as a share of the size of the numbers that place them, a line and a
     *  circle or two circles are taken to touch where two doubles find them to miss.
     */
    constexpr double Tight = 1e-12;

    /**
     *  The same in two doubles, as a share of the terms' product.
     */
    constexpr double Fine = 1e-26;

    /**
     *  The sum of two doubles, exactly.
     *
     *  @param  first  One.
     *  @param  second The other.
     *  @return        Their sum.
     */
    auto Sum(double first, double second) -> Double
    {
      const Twofold::Two sum = Twofold::TwoSum(first, second);
      return {.High = sum.High, .Low = sum.Low};
    }

    /**
     *  The product of two doubles, exactly.
     *
     *  @param  first  One.
     *  @param  second The other.
     *  @return        Their product.
     */
    auto Product(double first, double second) -> Double
    {
      const Twofold::Two product = Twofold::TwoProduct(first, second);
      return {.High = product.High, .Low = product.Low};
    }

    /**
     *  The two points a way either side of a foot, or the foot alone, or none: as the square of the
     *  half chord says.
     *
     *  @param  foot   The foot, midway between the two.
     *  @param  half   The way from the foot to one of them, as long as 1 before it is scaled by the
     *                 half chord.
     *  @param  square The half chord squared, and the slack within which it is 0.
     *  @return        The points.
     */
    auto Chorded(Point foot, Vector half, const std::array<double, 2> &square) -> Hits
    {
      if (square[0] < -square[1])
      {
        return {};
      }
      if (square[0] <= square[1])
      {
        return {.At = {foot}, .Count = 1};
      }
      const Vector reach = half * std::sqrt(square[0]);
      return {.At = {foot - reach, foot + reach}, .Count = 2};
    }

    /**
     *  The points a line and a circle share.
     *
     *  @param  line   The line.
     *  @param  circle The circle.
     *  @return        The points.
     */
    auto Sliced(const Piece &line, const Piece &circle) -> Hits
    {
      const Vector way    = line.End - line.Start;
      const Vector unit   = way / Long(way);
      const Vector off    = circle.Center - line.Start;
      const Point  foot   = line.Start + (unit * Dot(off, unit));
      const double aside  = Cross(unit, off);
      const double square = (circle.Radius - aside) * (circle.Radius + aside);
      const double radius = circle.Radius * circle.Radius;
      Hits         hits   = Chorded(foot, unit, {square, Slack * radius});
      /**
       *  A circle far smaller than the numbers that place it touches within their rounding, however
       *  far that is as a share of its radius.
       */
      const double hair = Tight
                        * (std::abs(line.Start.X)
                           + std::abs(line.Start.Y)
                           + std::abs(circle.Center.X)
                           + std::abs(circle.Center.Y)
                           + circle.Radius);
      hits.Graze        = std::abs(std::abs(aside) - circle.Radius) <= hair;
      hits.Close        = hits.Graze || std::abs(square) <= Near * radius;
      hits.Touch        = foot;
      return hits;
    }

    /**
     *  The points two circles share.
     *
     *  @param  first  One circle.
     *  @param  second The other.
     *  @return        The points; none for two about one center.
     */
    auto Circled(const Piece &first, const Piece &second) -> Hits
    {
      const Vector between = second.Center - first.Center;
      const double apart2  = Dot(between, between);
      if (!(apart2 > 0.0))
      {
        return {};
      }
      const double apart = std::sqrt(apart2);
      const double along =
          (apart2 + ((first.Radius - second.Radius) * (first.Radius + second.Radius)))
          / (2.0 * apart);
      const double square = (first.Radius - along) * (first.Radius + along);
      const double radius = first.Radius * first.Radius;
      const Point  foot   = first.Center + (between * (along / apart));
      Hits         hits =
          Chorded(foot, Vector{.X = between.Y, .Y = -between.X} / apart, {square, Slack * radius});
      /**
       *  Two circles about centers within a thousandth of a radius of each other cross where
       *  doubles place the point poorly, by as many times their rounding as the radius is the
       *  centers' distance: they too are met again in two doubles a number.
       */
      const double hair = Tight
                        * (std::abs(first.Center.X)
                           + std::abs(first.Center.Y)
                           + std::abs(second.Center.X)
                           + std::abs(second.Center.Y)
                           + first.Radius
                           + second.Radius);
      hits.Graze        = std::abs(std::abs(along) - first.Radius) * apart
                       <= hair * (apart + first.Radius + second.Radius);
      hits.Close = hits.Graze || std::abs(square) <= Near * radius || apart2 <= 1e-6 * radius;
      hits.Touch = foot;
      return hits;
    }

    /**
     *  The points a line's carrier and a circle's share, in two doubles.
     *
     *  @param  line   The line.
     *  @param  circle The circle.
     *  @return        The points.
     */
    auto Sliced(const Carrier &line, const Carrier &circle) -> Hits
    {
      const Double wayX    = line.WayX;
      const Double wayY    = line.WayY;
      const Double offX    = Difference(circle.X, line.X);
      const Double offY    = Difference(circle.Y, line.Y);
      const Double length2 = Add(Multiply(wayX, wayX), Multiply(wayY, wayY));
      const Double along   = Add(Multiply(offX, wayX), Multiply(offY, wayY));
      const Double aside   = Difference(Multiply(wayX, offY), Multiply(wayY, offX));
      const Double square  = Difference(Multiply(Multiply(circle.Radius, circle.Radius), length2),
                                        Multiply(aside, aside));
      const double slack   = Fine * circle.Radius.High * circle.Radius.High * length2.High;
      const Double place   = Divide(along, length2);
      const Double footX   = Add(line.X, Multiply(wayX, place));
      const Double footY   = Add(line.Y, Multiply(wayY, place));
      if (square.High < -slack)
      {
        return {};
      }
      if (square.High <= slack)
      {
        return {.At = {Point{.X = footX.High, .Y = footY.High}}, .Count = 1};
      }
      const Double half = Divide(Root(square), length2);
      return {.At    = {Point{.X = Difference(footX, Multiply(wayX, half)).High,
                              .Y = Difference(footY, Multiply(wayY, half)).High},
                        Point{.X = Add(footX, Multiply(wayX, half)).High,
                              .Y = Add(footY, Multiply(wayY, half)).High}},
              .Count = 2};
    }

    /**
     *  The points two circles' carriers share, in two doubles.
     *
     *  @param  first  One circle.
     *  @param  second The other.
     *  @return        The points; none for two about one center.
     */
    auto Circled(const Carrier &first, const Carrier &second) -> Hits
    {
      const Double betweenX = Difference(second.X, first.X);
      const Double betweenY = Difference(second.Y, first.Y);
      const Double apart2   = Add(Multiply(betweenX, betweenX), Multiply(betweenY, betweenY));
      if (!(apart2.High > 0.0))
      {
        return {};
      }
      const Double radius2 = Multiply(first.Radius, first.Radius);
      const Double along   = Divide(
          Difference(Add(apart2, radius2), Multiply(second.Radius, second.Radius)), Lifted(2.0));
      const Double square = Difference(Multiply(radius2, apart2), Multiply(along, along));
      const double slack  = Fine * first.Radius.High * first.Radius.High * apart2.High;
      const Double place  = Divide(along, apart2);
      const Double footX  = Add(first.X, Multiply(betweenX, place));
      const Double footY  = Add(first.Y, Multiply(betweenY, place));
      if (square.High < -slack)
      {
        return {};
      }
      if (square.High <= slack)
      {
        return {.At = {Point{.X = footX.High, .Y = footY.High}}, .Count = 1};
      }
      const Double half = Divide(Root(square), apart2);
      return {.At    = {Point{.X = Difference(footX, Multiply(betweenY, half)).High,
                              .Y = Add(footY, Multiply(betweenX, half)).High},
                        Point{.X = Add(footX, Multiply(betweenY, half)).High,
                              .Y = Difference(footY, Multiply(betweenX, half)).High}},
              .Count = 2};
    }

    /**
     *  The point two lines' carriers share, in two doubles.
     *
     *  @param  first  One line.
     *  @param  second The other.
     *  @return        The point; none for two along one direction.
     */
    auto Lined(const Carrier &first, const Carrier &second) -> Hits
    {
      const Double cross =
          Difference(Multiply(first.WayX, second.WayY), Multiply(first.WayY, second.WayX));
      if (!(std::abs(cross.High) > 0.0))
      {
        return {};
      }
      const Double offX = Difference(second.X, first.X);
      const Double offY = Difference(second.Y, first.Y);
      const Double along =
          Divide(Difference(Multiply(offX, second.WayY), Multiply(offY, second.WayX)), cross);
      return {.At    = {Point{.X = Add(first.X, Multiply(first.WayX, along)).High,
                              .Y = Add(first.Y, Multiply(first.WayY, along)).High}},
              .Count = 1};
    }

    /**
     *  A point in two doubles a coordinate.
     */
    struct Spot
    {
      Double X = {}; /**< Its x. */
      Double Y = {}; /**< Its y. */
    };

    /**
     *  The point of a carrier nearest a point.
     *
     *  @param  carrier The carrier.
     *  @param  point   The point, not a circle's center.
     *  @return         The carrier's point.
     */
    auto Footed(const Carrier &carrier, Point point) -> Spot
    {
      const Double offX = Difference(Lifted(point.X), carrier.X);
      const Double offY = Difference(Lifted(point.Y), carrier.Y);
      if (carrier.Round)
      {
        const Double scale =
            Divide(carrier.Radius, Root(Add(Multiply(offX, offX), Multiply(offY, offY))));
        return {.X = Add(carrier.X, Multiply(offX, scale)),
                .Y = Add(carrier.Y, Multiply(offY, scale))};
      }
      const Double wayX  = carrier.WayX;
      const Double wayY  = carrier.WayY;
      const Double place = Divide(Add(Multiply(offX, wayX), Multiply(offY, wayY)),
                                  Add(Multiply(wayX, wayX), Multiply(wayY, wayY)));
      return {.X = Add(carrier.X, Multiply(wayX, place)),
              .Y = Add(carrier.Y, Multiply(wayY, place))};
    }

    /**
     *  The square of how far a spot is from a point, less a square.
     *
     *  @param  spot   The spot.
     *  @param  point  The point.
     *  @param  square The square taken off.
     *  @return        The difference.
     */
    auto Ended(const Spot &spot, Point point, Double square) -> Double
    {
      const Double offX = Difference(spot.X, Lifted(point.X));
      const Double offY = Difference(spot.Y, Lifted(point.Y));
      return Difference(Add(Multiply(offX, offX), Multiply(offY, offY)), square);
    }
  } // namespace

  auto Met(const Piece &first, const Piece &second) -> Hits
  {
    /**
     *  A nearly straight arc's center lies too far off for doubles to meet its circle from: it is
     *  always met again in two doubles a number.
     */
    if (first.Shallow || second.Shallow)
    {
      return {.Close = true};
    }
    if (first.Round && second.Round)
    {
      return Circled(first, second);
    }
    if (first.Round != second.Round)
    {
      return first.Round ? Sliced(second, first) : Sliced(first, second);
    }
    const Vector way   = first.End - first.Start;
    const Vector other = second.End - second.Start;
    const double cross = Cross(way, other);
    if (!(std::abs(cross) > 1e-14 * std::sqrt(Dot(way, way) * Dot(other, other))))
    {
      return {};
    }
    const double along = Cross(second.Start - first.Start, other) / cross;
    Point        point = first.Start + (way * along);
    /**
     *  A level or an upright line holds one of the point's two numbers itself, exactly.
     */
    for (const Piece *line : {&first, &second})
    {
      if (!(std::abs(line->End.X - line->Start.X) > 0.0))
      {
        point.X = line->Start.X;
      }
      if (!(std::abs(line->End.Y - line->Start.Y) > 0.0))
      {
        point.Y = line->Start.Y;
      }
    }
    /**
     *  Two lines within a thousandth of one direction cross where doubles place the point poorly:
     *  they are met again in two doubles a number.
     */
    return {.At    = {point},
            .Count = 1,
            .Close = std::abs(cross) <= 1e-3 * std::sqrt(Dot(way, way) * Dot(other, other)),
            .Touch = point};
  }

  auto Carried(const Piece &piece, double bulge) -> Carrier
  {
    const Double chordX = Sum(piece.End.X, -piece.Start.X);
    const Double chordY = Sum(piece.End.Y, -piece.Start.Y);
    if (!piece.Round)
    {
      return {.Round = false,
              .X     = Lifted(piece.Start.X),
              .Y     = Lifted(piece.Start.Y),
              .WayX  = chordX,
              .WayY  = chordY};
    }
    const Double factor =
        Divide(Difference(Lifted(1.0), Product(bulge, bulge)), Lifted(4.0 * bulge));
    const Double centerX =
        Difference(Divide(Sum(piece.Start.X, piece.End.X), Lifted(2.0)), Multiply(chordY, factor));
    const Double centerY =
        Add(Divide(Sum(piece.Start.Y, piece.End.Y), Lifted(2.0)), Multiply(chordX, factor));
    const Double armX = Difference(Lifted(piece.Start.X), centerX);
    const Double armY = Difference(Lifted(piece.Start.Y), centerY);
    return {.Round  = true,
            .X      = centerX,
            .Y      = centerY,
            .Radius = Root(Add(Multiply(armX, armX), Multiply(armY, armY)))};
  }

  auto Shifted(const Carrier &carrier, bool left, double by) -> Carrier
  {
    Carrier moved = carrier;
    if (carrier.Round)
    {
      const Double reach = Lifted(std::abs(by));
      moved.Radius =
          (by > 0.0) == left ? Difference(carrier.Radius, reach) : Add(carrier.Radius, reach);
      return moved;
    }
    const Double length =
        Root(Add(Multiply(carrier.WayX, carrier.WayX), Multiply(carrier.WayY, carrier.WayY)));
    const Double factor = Divide(Lifted(by), length);
    moved.X             = Difference(carrier.X, Multiply(carrier.WayY, factor));
    moved.Y             = Add(carrier.Y, Multiply(carrier.WayX, factor));
    return moved;
  }

  auto Met(const Carrier &first, const Carrier &second) -> Hits
  {
    if (!first.Round && !second.Round)
    {
      return Lined(first, second);
    }
    if (first.Round && second.Round)
    {
      return Circled(first, second);
    }
    return first.Round ? Sliced(second, first) : Sliced(first, second);
  }

  auto Clear(const Probe &probe, const Piece &piece, const Carrier &wide) -> Double
  {
    const Spot   spot   = Footed(*probe.On, probe.Near);
    const Double square = Product(probe.Distance, probe.Distance);
    if (!piece.Round)
    {
      const Vector way   = piece.End - piece.Start;
      const double along = Dot(probe.Near - piece.Start, way) / Dot(way, way);
      if (along <= 0.0 || along >= 1.0)
      {
        return Ended(spot, along <= 0.0 ? piece.Start : piece.End, square);
      }
      const Double cross = Difference(Multiply(wide.WayX, Difference(spot.Y, wide.Y)),
                                      Multiply(wide.WayY, Difference(spot.X, wide.X)));
      return Difference(Divide(Multiply(cross, cross),
                               Add(Multiply(wide.WayX, wide.WayX), Multiply(wide.WayY, wide.WayY))),
                        square);
    }
    if (Facing(piece, probe.Near))
    {
      const Double offX = Difference(spot.X, wide.X);
      const Double offY = Difference(spot.Y, wide.Y);
      const Double gap =
          Difference(Root(Add(Multiply(offX, offX), Multiply(offY, offY))), wide.Radius);
      return Difference(Multiply(gap, gap), square);
    }
    const Double start = Ended(spot, piece.Start, square);
    const Double end   = Ended(spot, piece.End, square);
    return end.High < start.High ? end : start;
  }

  auto Sharpened(const Carrier &first, const Carrier &second, Point near) -> Point
  {
    /**
     *  How far off a carrier a point lies, to a scale, and the way that grows fastest.
     */
    struct Miss
    {
      double Value = 0.0; /**< A line: the cross of its way with the way to the point. A circle:
                               the square of the point's distance from its center less its
                               radius squared.                                              */
      Vector Slope = {};  /**< How the value changes with the point. */
    };

    const auto missed = [](const Carrier &carrier, Point point) -> Miss
    {
      const Double offX = Difference(Lifted(point.X), carrier.X);
      const Double offY = Difference(Lifted(point.Y), carrier.Y);
      if (carrier.Round)
      {
        return {.Value = Difference(Add(Multiply(offX, offX), Multiply(offY, offY)),
                                    Multiply(carrier.Radius, carrier.Radius))
                             .High,
                .Slope = {.X = 2.0 * offX.High, .Y = 2.0 * offY.High}};
      }
      return {.Value = Difference(Multiply(carrier.WayX, offY), Multiply(carrier.WayY, offX)).High,
              .Slope = {.X = -carrier.WayY.High, .Y = carrier.WayX.High}};
    };
    /**
     *  One step: the doubles lie within a few of their own steps of the point, and what a second
     *  would add is the square of that.
     */
    const Miss   one   = missed(first, near);
    const Miss   other = missed(second, near);
    const double turn  = Cross(one.Slope, other.Slope);
    /**
     *  Two that run one way there are stepped onto no further.
     */
    if (!(std::abs(turn) > 0.0))
    {
      return near;
    }
    const Vector move = {.X = ((other.Value * one.Slope.Y) - (one.Value * other.Slope.Y)) / turn,
                         .Y = ((one.Value * other.Slope.X) - (other.Value * one.Slope.X)) / turn};
    const double far  = std::abs(move.X) + std::abs(move.Y);
    /**
     *  A step out of all measure with the point, as two that all but touch give, is not taken.
     */
    return far <= 0x1p-20 * (std::abs(near.X) + std::abs(near.Y) + far) ? near + move : near;
  }

  auto Clear(const Probe &probe, Point point) -> Double
  {
    return Ended(Footed(*probe.On, probe.Near), point, Product(probe.Distance, probe.Distance));
  }

  auto Beside(const Probe &probe, const Carrier &edge) -> Double
  {
    const Spot   spot = Footed(*probe.On, probe.Near);
    const Double offX = Difference(spot.X, edge.X);
    const Double offY = Difference(spot.Y, edge.Y);
    if (edge.Round)
    {
      return Difference(Add(Multiply(offX, offX), Multiply(offY, offY)),
                        Multiply(edge.Radius, edge.Radius));
    }
    return Difference(Multiply(edge.WayX, offY), Multiply(edge.WayY, offX));
  }
} // namespace Geoxy::Belt
