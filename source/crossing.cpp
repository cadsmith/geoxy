/***************************************************************************************************
 *  @file      source/crossing.cpp
 *  @brief     Where two lines or conics cross, worked out exactly and rounded once.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "crossing.hpp"

#include "bow.hpp"
#include "crossings/bows.hpp"
#include "elliptic.hpp"
#include "exact/whole.hpp"
#include "expansion.hpp"
#include "met/scale.hpp"
#include "real.hpp"
#include "sure.hpp"
#include "tangency.hpp"
#include "twofold.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/nearest.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <numbers>
#include <optional>
#include <utility>
#include <variant>

namespace Geoxy::Crossings
{
  namespace
  {
    using Expansion::Difference;
    using Expansion::Number;
    using Expansion::Of;
    using Expansion::Product;
    using Expansion::Quotient;
    using Expansion::Sign;
    using Expansion::Sum;
    using Twofold::Add;
    using Twofold::Difference;
    using Twofold::Divide;
    using Twofold::Double;
    using Twofold::Lifted;
    using Twofold::Multiply;
    using Twofold::Negated;
    using Twofold::Root;
    using Twofold::Two;
    using Twofold::TwoProduct;
    using Twofold::TwoSum;

    /**
     *  A pair of exact numbers: a point or a vector.
     */
    using Pair = std::array<Number, 2>;

    /**
     *  A line, exactly: a point on it and its way.
     */
    struct Straight
    {
      Pair Through = {}; /**< A point on it. */
      Pair Way     = {}; /**< Its way.       */
    };

    /**
     *  An ellipse's carrier, exactly: the points P where |adj (P - origin)|^2 = det^2, adj and det
     *  those of the map of the unit circle; and, for a circle, its center and squared radius.
     */
    struct Conic
    {
      Pair                Origin   = {}; /**< The center. */
      std::array<Pair, 2> Adjugate = {}; /**< The map's adjugate, a row a pair; a
                                              circle's is not asked for.          */
      Number Det     = {};               /**< The map's determinant, likewise. */
      bool   Circle  = false;            /**< Whether it is a circle.          */
      Number Squared = {};               /**< A circle's squared radius.       */
    };

    /**
     *  A polyline's bulge, exactly: its two vertices and its bulge b. Its circle is the rational
     *  one through both, its center m + k n, m the chord's middle, n the chord turned a quarter to
     *  the left and k = (1 - b^2) / (4 b): the points X where 4 b (X - P0).(X - P1) - 2 (1 - b^2)
     *  (X - P0).n is zero, with no radius and no division.
     */
    struct Bulge
    {
      Pair   Start = {}; /**< Its first vertex, P0.   */
      Pair   End   = {}; /**< Its second vertex, P1.  */
      Number Bend  = {}; /**< Its bulge, b, not zero. */
    };

    /**
     *  The carrier of a part of a curve: a line, a conic or a bulge's circle.
     */
    using Carrier = std::variant<Straight, Conic, Bulge>;

    /**
     *  The dot product of two exact pairs.
     *
     *  @param  first  One.
     *  @param  second The other.
     *  @return        The product.
     */
    auto Dot(const Pair &first, const Pair &second) -> Number
    {
      return Sum(Product(first[0], second[0]), Product(first[1], second[1]));
    }

    /**
     *  The cross product of two exact pairs.
     *
     *  @param  first  One.
     *  @param  second The other.
     *  @return        The product.
     */
    auto Cross(const Pair &first, const Pair &second) -> Number
    {
      return Difference(Product(first[0], second[1]), Product(first[1], second[0]));
    }

    /**
     *  A point as an exact pair.
     *
     *  @param  point The point.
     *  @return       The pair.
     */
    auto Exact(Point point) -> Pair
    {
      return {Of(point.X), Of(point.Y)};
    }

    /**
     *  A vector as an exact pair.
     *
     *  @param  vector The vector.
     *  @return        The pair.
     */
    auto Exact(Vector vector) -> Pair
    {
      return {Of(vector.X), Of(vector.Y)};
    }

    /**
     *  The difference of two points, exactly.
     *
     *  @param  head The point.
     *  @param  tail The point taken from it.
     *  @return      The difference.
     */
    auto Between(Point head, Point tail) -> Pair
    {
      return {Difference(Of(head.X), Of(tail.X)), Difference(Of(head.Y), Of(tail.Y))};
    }

    /**
     *  The carrier of the ellipse a map of the unit circle draws.
     *
     *  @param  origin Where the unit circle's center goes.
     *  @param  axes   Where its x and y axes go.
     *  @return        The carrier.
     */
    auto Conical(Point origin, const std::array<Vector, 2> &axes) -> Conic
    {
      const Vector &x      = axes[0];
      const Vector &y      = axes[1];
      Conic         conic  = {.Origin   = Exact(origin),
                              .Adjugate = {Pair{Of(y.Y), Of(-y.X)}, Pair{Of(-x.Y), Of(x.X)}},
                              .Det      = Cross(Exact(x), Exact(y)),
                              .Circle   = false,
                              .Squared  = {}};
      const Number  dot    = Dot(Exact(x), Exact(y));
      const Number  across = Dot(Exact(x), Exact(x));
      conic.Circle  = Sign(dot) == 0 && Sign(Difference(across, Dot(Exact(y), Exact(y)))) == 0;
      conic.Squared = across;
      return conic;
    }

    /**
     *  Whether a place names a polyline's segment: the polyline has one, and the place is a finite
     *  number.
     *
     *  @param  polyline The polyline.
     *  @param  along    The place.
     *  @return          Whether it does.
     */
    auto Seatable(const Polyline &polyline, double along) -> bool
    {
      return Count(polyline) > 0 && std::isfinite(along);
    }

    /**
     *  A polyline's segment a place names: the whole part of the place, the first segment before it
     *  and the last past it.
     *
     *  @pre             The place names one, as Seatable tells.
     *  @param  polyline The polyline.
     *  @param  along    The place.
     *  @return          The segment's two vertices.
     */
    auto Seated(const Polyline &polyline, double along) -> std::array<const Vertex *, 2>
    {
      const auto index = static_cast<std::size_t>(
          std::clamp(std::floor(along), 0.0, static_cast<double>(Count(polyline) - 1)));
      return {&polyline.Vertices[index],
              &polyline.Vertices[(index + 1) % polyline.Vertices.size()]};
    }

    /**
     *  Whether a polyline's segment is one Nth makes a line, told without making its arc.
     *
     *  @param  start Its first vertex.
     *  @param  end   Its second.
     *  @return       Whether it is.
     */
    auto Flatly(const Vertex &start, const Vertex &end) -> bool
    {
      return !(std::abs(start.Bulge) >= Geoxy::Straight) || start.Position == end.Position;
    }

    /**
     *  What carries a curve's part at a place, as Carried finds it.
     */
    enum class Bearer : std::uint8_t
    {
      None,  /**< Nothing: a Bézier curve, a spline, a polyline with no segment there. */
      Line,  /**< A line.                                                              */
      Conic, /**< An arc's circle or an ellipse.                                       */
      Bow,   /**< A polyline's bulge's circle.                                         */
    };

    /**
     *  What carries a curve's part at a place, told without working the carrier out.
     *
     *  @param  curve The curve.
     *  @param  along The place.
     *  @return       What Carried finds there.
     */
    auto Borne(const Curve &curve, double along) -> Bearer
    {
      if (std::holds_alternative<Line>(curve)
          || std::holds_alternative<Ray>(curve)
          || std::holds_alternative<Xline>(curve))
      {
        return Bearer::Line;
      }
      if (std::holds_alternative<Arc>(curve) || std::holds_alternative<Ellipse>(curve))
      {
        return Bearer::Conic;
      }
      const auto *polyline = std::get_if<Polyline>(&curve);
      if (polyline == nullptr || !Seatable(*polyline, along))
      {
        return Bearer::None;
      }
      const std::array<const Vertex *, 2> seat = Seated(*polyline, along);
      return Flatly(*seat[0], *seat[1]) ? Bearer::Line : Bearer::Bow;
    }

    /**
     *  Whether Solved finds nothing where two carriers cross, told before they are worked out in
     *  expansions: a part with no carrier, or a bulge's circle met with anything but a line.
     *
     *  @param  one   What carries the first curve's part.
     *  @param  other What carries the second's.
     *  @return       Whether it does.
     */
    auto Unsolvable(Bearer one, Bearer other) -> bool
    {
      return one == Bearer::None
          || other == Bearer::None
          || (one != Bearer::Line
              && other != Bearer::Line
              && (one == Bearer::Bow || other == Bearer::Bow));
    }

    /**
     *  The carrier of a curve's part at a place: a line, a conic or a bulge's circle.
     *
     *  @pre          The curve's part there has one, as Borne tells.
     *  @param  curve The curve.
     *  @param  along The place.
     *  @return       The carrier.
     */
    auto Carried(const Curve &curve, double along) -> Carrier
    {
      if (const auto *line = std::get_if<Line>(&curve))
      {
        return Straight{.Through = Exact(line->Start), .Way = Between(line->End, line->Start)};
      }
      if (const auto *ray = std::get_if<Ray>(&curve))
      {
        return Straight{.Through = Exact(ray->Start), .Way = Exact(ray->Way)};
      }
      if (const auto *xline = std::get_if<Xline>(&curve))
      {
        return Straight{.Through = Exact(xline->Through), .Way = Exact(xline->Way)};
      }
      if (const auto *arc = std::get_if<Arc>(&curve))
      {
        /**
         *  A circle, as Conical would find it: its center and squared radius are all a circle's
         *  carrier is solved by.
         */
        return Conic{.Origin   = Exact(arc->Center),
                     .Adjugate = {},
                     .Det      = {},
                     .Circle   = true,
                     .Squared  = Product(Of(arc->Radius), Of(arc->Radius))};
      }
      if (const auto *ellipse = std::get_if<Ellipse>(&curve))
      {
        return Conical(ellipse->Map.Origin, {ellipse->Map.X, ellipse->Map.Y});
      }
      const std::array<const Vertex *, 2> seat  = Seated(std::get<Polyline>(curve), along);
      const Vertex                       &start = *seat[0];
      const Vertex                       &end   = *seat[1];
      if (Flatly(start, end))
      {
        return Straight{.Through = Exact(start.Position),
                        .Way     = Between(end.Position, start.Position)};
      }
      return Bulge{
          .Start = Exact(start.Position), .End = Exact(end.Position), .Bend = Of(start.Bulge)};
    }

    /**
     *  A curve's part at a place as whole numbers read a straight carrier: a line's or a polyline's
     *  straight segment's two ends, a ray's or an endless line's point and way.
     *
     *  @param  curve The curve.
     *  @param  along The place.
     *  @return       The carrier; nothing where the part there is not straight.
     */
    auto Ruled(const Curve &curve, double along) -> std::optional<Exact::Straightway>
    {
      if (const auto *line = std::get_if<Line>(&curve))
      {
        return Exact::Straightway{.Through = line->Start, .Toward = line->End, .Way = false};
      }
      if (const auto *ray = std::get_if<Ray>(&curve))
      {
        return Exact::Straightway{
            .Through = ray->Start, .Toward = {.X = ray->Way.X, .Y = ray->Way.Y}, .Way = true};
      }
      if (const auto *xline = std::get_if<Xline>(&curve))
      {
        return Exact::Straightway{.Through = xline->Through,
                                  .Toward  = {.X = xline->Way.X, .Y = xline->Way.Y},
                                  .Way     = true};
      }
      if (Borne(curve, along) != Bearer::Line)
      {
        return std::nullopt;
      }
      const std::array<const Vertex *, 2> seat = Seated(std::get<Polyline>(curve), along);
      return Exact::Straightway{
          .Through = seat[0]->Position, .Toward = seat[1]->Position, .Way = false};
    }

    /**
     *  Whether sums of doubles hold every product the exact crossing of two straight carriers works
     *  out and rounds by, so that it is exact. Of their numbers, every one below 2^a, and none
     *  setting a bit below 2^b, a step of a double below the least: each way and gap lies below
     *  2^(a + 1) on 2^b's grid, the cross products below 2^(2a + 3) on 2^2b's, each numerator below
     *  2^(3a + 5) on 2^3b's, which overflows for no a up to 339 and drops no bit for b from -358
     *  up. The coordinate, a quotient of the two, lies below 2^(3a + 5 - 2b), within the doubles
     *  where 3a less 2b is at most 1017; not zero, at or above 2^(3b - 2a - 3), normal where 3b
     *  less 2a is at least -1019. The halfway points Rounded compares it with set no bit below half
     *  a step of a double at its size taken a binade down, 2^(3b - 2a - 57), and their products
     *  with a cross product none below 2^(5b - 2a - 57), which drops none where 5b less 2a is at
     *  least -1016; the last implies the one before. Numbers of like size hold so from 2^-250 to
     *  2^339; numbers that part further, one a sliver beside the rest, may not.
     *
     *  @param  ways The two carriers, as whole numbers read them.
     *  @return      Whether they are held; so too where a number is not finite, or all are zero,
     *               which the sums decide as they always have.
     */
    auto Holdable(const std::array<Exact::Straightway, 2> &ways) -> bool
    {
      int most  = std::numeric_limits<int>::min();
      int least = std::numeric_limits<int>::max();
      for (const Exact::Straightway &way : ways)
      {
        for (const double value : {way.Through.X, way.Through.Y, way.Toward.X, way.Toward.Y})
        {
          if (!std::isfinite(value))
          {
            return true;
          }
          if (!Zero(value))
          {
            const int place = std::ilogb(value) + 1;
            most            = std::max(most, place);
            least           = std::min(least, place);
          }
        }
      }
      if (least == std::numeric_limits<int>::max())
      {
        return true;
      }
      const int high = most;
      const int low  = least - 53;
      return std::min({339 - high, low + 358, 1017 - ((3 * high) - (2 * low)),
                       (5 * low) - (2 * high) + 1016})
          >= 0;
    }

    /**
     *  A point that may hold a square root: each coordinate's quotient, the same root.
     */
    struct Candidate
    {
      Quotient X       = {};    /**< Its x.                                         */
      Quotient Y       = {};    /**< Its y.                                         */
      bool     Root    = false; /**< Whether it has one.                            */
      bool     Tangent = false; /**< Whether the root is of zero: one point, twice. */
    };

    /**
     *  Where two lines cross: the one point, nothing where they are parallel.
     *
     *  @param  first  One line.
     *  @param  second The other.
     *  @return        The point.
     */
    auto Solved(const Straight &first, const Straight &second) -> std::optional<Candidate>
    {
      const Number turn = Cross(first.Way, second.Way);
      if (Sign(turn) == 0)
      {
        return std::nullopt;
      }
      const Pair   gap   = {Difference(second.Through[0], first.Through[0]),
                            Difference(second.Through[1], first.Through[1])};
      const Number along = Cross(gap, second.Way);
      Candidate    out;
      for (std::size_t axis = 0; axis < 2; ++axis)
      {
        Quotient &quotient = axis == 0 ? out.X : out.Y;
        quotient.Top =
            Sum(Product(first.Through.at(axis), turn), Product(first.Way.at(axis), along));
        quotient.Below = turn;
      }
      return out;
    }

    /**
     *  Where a line crosses a conic: the points t of the line where |adj (P + t d - O)|^2 = det^2,
     *  a quadratic a t^2 + 2 b t + c, the roots (-b +- sqrt(b^2 - a c)) / a; for a circle, where |P
     *  + t d - O|^2 = r^2.
     *
     *  @param  line  The line.
     *  @param  conic The conic.
     *  @return       The points, the root's sign to choose; nothing where they do not meet.
     */
    auto Solved(const Straight &line, const Conic &conic) -> std::optional<Candidate>
    {
      const Pair gap = {Difference(line.Through[0], conic.Origin[0]),
                        Difference(line.Through[1], conic.Origin[1])};
      /**
       *  A circle's adjugate is its radius times a turn, its determinant plus or minus its squared
       *  radius: the quadratic is the circle's own, |P + t d - O|^2 = r^2, times r^2, and its roots
       *  are the same, worked in numbers four factors of r shorter.
       */
      const auto coefficients = [&line, &conic, &gap] -> std::array<Number, 3>
      {
        if (conic.Circle)
        {
          return {Dot(line.Way, line.Way), Dot(gap, line.Way),
                  Difference(Dot(gap, gap), conic.Squared)};
        }
        const Pair p = {Dot(conic.Adjugate[0], gap), Dot(conic.Adjugate[1], gap)};
        const Pair q = {Dot(conic.Adjugate[0], line.Way), Dot(conic.Adjugate[1], line.Way)};
        return {Dot(q, q), Dot(p, q), Difference(Dot(p, p), Product(conic.Det, conic.Det))};
      };
      const std::array<Number, 3> terms = coefficients();
      const Number               &a     = terms[0];
      const Number               &b     = terms[1];
      const Number               &c     = terms[2];
      const Number                d     = Difference(Product(b, b), Product(a, c));
      if (Sign(a) == 0 || Sign(d) < 0)
      {
        return std::nullopt;
      }
      Candidate out;
      out.Root    = true;
      out.Tangent = Sign(d) == 0;
      for (std::size_t axis = 0; axis < 2; ++axis)
      {
        Quotient &quotient = axis == 0 ? out.X : out.Y;
        quotient.Top = Difference(Product(line.Through.at(axis), a), Product(line.Way.at(axis), b));
        quotient.Bottom = line.Way.at(axis);
        quotient.Root   = d;
        quotient.Below  = a;
      }
      return out;
    }

    /**
     *  Where a line crosses a bulge's circle: its quadratic times 4 b, so that no division is left,
     *  a t^2 + 2 c t + e with a = 4 b d.d, c = 2 b d.(h + j) - (1 - b^2) d.n and e = 4 b h.j - 2 (1
     *  - b^2) h.n, h = P - P0, j = P - P1 and n the chord turned a quarter to the left; the roots
     *    (-c +- sqrt(c^2 - a e)) / a, the radius never worked out.
     *
     *  @param  line  The line.
     *  @param  bulge The bulge.
     *  @return       The points, the root's sign to choose; nothing where they do not meet.
     */
    auto Solved(const Straight &line, const Bulge &bulge) -> std::optional<Candidate>
    {
      const Pair   from   = {Difference(line.Through[0], bulge.Start[0]),
                             Difference(line.Through[1], bulge.Start[1])};
      const Pair   onto   = {Difference(line.Through[0], bulge.End[0]),
                             Difference(line.Through[1], bulge.End[1])};
      const Pair   both   = {Sum(from[0], onto[0]), Sum(from[1], onto[1])};
      const Pair   across = {Difference(bulge.Start[1], bulge.End[1]),
                             Difference(bulge.End[0], bulge.Start[0])};
      const Number twice  = Sum(bulge.Bend, bulge.Bend);
      const Number four   = Sum(twice, twice);
      const Number rest   = Difference(Of(1.0), Product(bulge.Bend, bulge.Bend));
      const Number a      = Product(four, Dot(line.Way, line.Way));
      const Number b =
          Difference(Product(twice, Dot(line.Way, both)), Product(rest, Dot(line.Way, across)));
      const Number c =
          Difference(Product(four, Dot(from, onto)), Product(Sum(rest, rest), Dot(from, across)));
      const Number d = Difference(Product(b, b), Product(a, c));
      if (Sign(a) == 0 || Sign(d) < 0)
      {
        return std::nullopt;
      }
      Candidate out;
      out.Root    = true;
      out.Tangent = Sign(d) == 0;
      for (std::size_t axis = 0; axis < 2; ++axis)
      {
        Quotient &quotient = axis == 0 ? out.X : out.Y;
        quotient.Top = Difference(Product(line.Through.at(axis), a), Product(line.Way.at(axis), b));
        quotient.Bottom = line.Way.at(axis);
        quotient.Root   = d;
        quotient.Below  = a;
      }
      return out;
    }

    /**
     *  Where two circles cross: on their radical line, c1 + (k d +- sqrt(K) d') / (2 |d|^2), with d
     *  = c2 - c1, d' it turned a quarter, k (the chord below) = |d|^2 + r1^2 - r2^2 and K = 4 |d|^2
     *  r1^2 - k^2.
     *
     *  @param  first  One circle.
     *  @param  second The other.
     *  @return        The points, nothing where they do not meet or are one.
     */
    auto Solved(const Conic &first, const Conic &second) -> std::optional<Candidate>
    {
      if (!first.Circle || !second.Circle)
      {
        return std::nullopt;
      }
      const Pair   gap   = {Difference(second.Origin[0], first.Origin[0]),
                            Difference(second.Origin[1], first.Origin[1])};
      const Number apart = Dot(gap, gap);
      const Number chord = Difference(Sum(apart, first.Squared), second.Squared);
      const Number big =
          Difference(Product(Product(Of(4.0), apart), first.Squared), Product(chord, chord));
      if (Sign(apart) == 0 || Sign(big) < 0)
      {
        return std::nullopt;
      }
      Candidate out;
      out.Root           = true;
      out.Tangent        = Sign(big) == 0;
      const Number twice = Sum(apart, apart);
      out.X              = {.Top    = Sum(Product(first.Origin[0], twice), Product(chord, gap[0])),
                            .Bottom = Difference(Number{}, gap[1]),
                            .Root   = big,
                            .Below  = twice};
      out.Y              = {.Top    = Sum(Product(first.Origin[1], twice), Product(chord, gap[1])),
                            .Bottom = gap[0],
                            .Root   = big,
                            .Below  = twice};
      return out;
    }

    /**
     *  Where two bulges' circles cross, as Crossbows works it out, exactly: P0 + (8 b B A + (A'.q
     *  +- sqrt(D)) A') / (8 b |A|^2), P0 the first's start.
     *
     *  @param  first  One bulge, met with the radical line.
     *  @param  second The other.
     *  @return        The points, the root's sign to choose; nothing where they do not meet or are
     *                 of one center.
     */
    auto Solved(const Bulge &first, const Bulge &second) -> std::optional<Candidate>
    {
      const Bulge &small  = first;
      const Bulge &large  = second;
      const auto   pulled = [](const Bulge &bulge) -> Pair
      {
        const Pair   chord = {Difference(bulge.End[0], bulge.Start[0]),
                              Difference(bulge.End[1], bulge.Start[1])};
        const Number four  = Expansion::Scaled(bulge.Bend, 4.0);
        const Number rest =
            Expansion::Scaled(Difference(Of(1.0), Product(bulge.Bend, bulge.Bend)), 2.0);
        return {Difference(Product(four, chord[0]), Product(rest, chord[1])),
                Sum(Product(four, chord[1]), Product(rest, chord[0]))};
      };
      const Pair   q     = pulled(small);
      const Pair   p     = pulled(large);
      const Pair   gap   = {Difference(large.Start[0], small.Start[0]),
                            Difference(large.Start[1], small.Start[1])};
      const Number both  = Expansion::Scaled(Product(small.Bend, large.Bend), 8.0);
      const Pair   a     = {Sum(Difference(Product(both, gap[0]), Product(large.Bend, q[0])),
                                Product(small.Bend, p[0])),
                            Sum(Difference(Product(both, gap[1]), Product(large.Bend, q[1])),
                                Product(small.Bend, p[1]))};
      const Number lead  = Sum(Product(Expansion::Scaled(both, 0.5), Dot(gap, gap)),
                               Product(small.Bend, Dot(gap, p)));
      const Pair   turn  = {Difference(Number{}, a[1]), a[0]};
      const Number span  = Dot(a, a);
      const Number eight = Expansion::Scaled(small.Bend, 8.0);
      const Number lean  = Difference(Dot(a, q), Product(eight, lead));
      const Number disc  = Difference(Product(span, Dot(q, q)), Product(lean, lean));
      if (Sign(span) == 0 || Sign(disc) < 0)
      {
        return std::nullopt;
      }
      const Number below = Product(eight, span);
      const Number foot  = Product(eight, lead);
      const Number along = Dot(turn, q);
      Candidate    out;
      out.Root    = true;
      out.Tangent = Sign(disc) == 0;
      for (std::size_t axis = 0; axis < 2; ++axis)
      {
        Quotient &quotient = axis == 0 ? out.X : out.Y;
        quotient.Top       = Sum(Product(small.Start.at(axis), below),
                                 Sum(Product(foot, a.at(axis)), Product(along, turn.at(axis))));
        quotient.Bottom    = turn.at(axis);
        quotient.Root      = disc;
        quotient.Below     = below;
      }
      return out;
    }

    /**
     *  Where two carriers cross, the first's points first.
     *
     *  @pre           Solved finds something for them, as Unsolvable tells: a line and any carrier,
     *                 or two conics.
     *  @param  first  One carrier.
     *  @param  second The other.
     *  @return        The points.
     */
    auto Solved(const Carrier &first, const Carrier &second) -> std::optional<Candidate>
    {
      const auto *one   = std::get_if<Straight>(&first);
      const auto *other = std::get_if<Straight>(&second);
      if (one != nullptr && other != nullptr)
      {
        return Solved(*one, *other);
      }
      if (one != nullptr || other != nullptr)
      {
        const Straight &line = one != nullptr ? *one : *other;
        const Carrier  &rest = one != nullptr ? second : first;
        if (const auto *conic = std::get_if<Conic>(&rest))
        {
          return Solved(line, *conic);
        }
        return Solved(line, std::get<Bulge>(rest));
      }
      return Solved(std::get<Conic>(first), std::get<Conic>(second));
    }

    /**
     *  A point's place on a line where the foot is nearest it by far: the foot, as Nearest would
     *  take it, found with no root taken, by squared distances. Where the foot's is below the
     *  nearer end's by more than their roundings could close, its distance as Nearest works it out
     *  is below the ends'. A foot past an end is taken at that end, whose own distance it then
     *  equals, and so is left to Nearest; so is a point at an end, or so near one that its square
     *  falls below the doubles' middle range. Above that range the foot's is held below 2^900.
     *
     *  @param  line  The line.
     *  @param  along The foot's place, as Footing works it out.
     *  @param  point The point.
     *  @return       The place; nothing where Nearest must decide.
     */
    auto Footed(const Line &line, double along, Point point) -> std::optional<Place>
    {
      constexpr double Low    = 0x1p-900;
      constexpr double High   = 0x1p900;
      constexpr double Margin = 1.0 - 0x1p-40;
      const double     within = std::clamp(along, 0.0, 1.0);
      const Point      foot   = Lerp(line.Start, line.End, within);
      const Vector     off    = point - foot;
      const Vector     start  = point - Lerp(line.Start, line.End, 0.0);
      const Vector     end    = point - Lerp(line.Start, line.End, 1.0);
      const double     square = Dot(off, off);
      const double     nearer = std::min(Dot(start, start), Dot(end, end));
      if (Low <= nearer && square < std::min(nearer, High) * Margin)
      {
        return Place{.At = foot, .Along = within};
      }
      return std::nullopt;
    }

    /**
     *  The angle on an ellipse's unit circle of a point, in radians: the heading of adj(M) (p - o),
     *  turned by the sign of det(M), which points where M^-1 (p - o) does, adj(M) and det(M) those
     *  of the ellipse's map M and o its center. The difference p - o is exact in two doubles and
     *  the products in double-double, so each coordinate of adj(M) (p - o) is the double nearest it
     *  but for a 2^-100 of its terms; the heading is the double nearest the angle of those.
     *
     *  @param  ellipse The ellipse.
     *  @param  point   The point.
     *  @return         The angle; nothing for an ellipse flat to within 2^-90 of its map's terms,
     *                  or for a point at its center.
     */
    auto Facing(const Ellipse &ellipse, Point point) -> std::optional<double>
    {
      const Vector x   = ellipse.Map.X;
      const Vector y   = ellipse.Map.Y;
      const Two    one = TwoProduct(x.X, y.Y);
      const Two    two = TwoProduct(x.Y, y.X);
      const Double det =
          Difference({.High = one.High, .Low = one.Low}, {.High = two.High, .Low = two.Low});
      const double terms = std::abs(one.High) + std::abs(two.High);
      if (!(std::abs(det.High) > 0x1p-90 * terms))
      {
        return std::nullopt;
      }
      const Two    wide   = TwoSum(point.X, -ellipse.Map.Origin.X);
      const Two    high   = TwoSum(point.Y, -ellipse.Map.Origin.Y);
      const Double gapX   = {.High = wide.High, .Low = wide.Low};
      const Double gapY   = {.High = high.High, .Low = high.Low};
      const Double across = Difference(Multiply(Lifted(y.Y), gapX), Multiply(Lifted(y.X), gapY));
      const Double upward = Difference(Multiply(Lifted(x.X), gapY), Multiply(Lifted(x.Y), gapX));
      const double sign   = det.High > 0.0 ? 1.0 : -1.0;
      const Vector seen   = {.X = sign * across.High, .Y = sign * upward.High};
      if (Zero(seen.X) && Zero(seen.Y))
      {
        return std::nullopt;
      }
      return Heading(seen).InRadians();
    }

    /**
     *  How far past a line's or a segment's ends a point on its carrier may lie and be its: a few
     *  steps of a double, as the point is the true one rounded.
     */
    constexpr double Slack = 16.0 * std::numeric_limits<double>::epsilon();

    /**
     *  A point's place along a way from a start, as a fraction of the way; no number along a way
     *  shorter than the square root of the least double, whose square is nothing, so that the point
     *  is on no such run and Intersect keeps the hit as it found it.
     *
     *  @param  start Where the way starts.
     *  @param  way   The way.
     *  @param  there The point.
     *  @return       The fraction; NaN along a way of no square.
     */
    auto Inside(Point start, Vector way, Point there) -> double
    {
      const double square = Dot(way, way);
      if (!(square > 0.0))
      {
        return std::numeric_limits<double>::quiet_NaN();
      }
      return Dot(there - start, way) / square;
    }

    /**
     *  A place, where a point lies on its curve.
     *
     *  @param  onIt  Whether it does.
     *  @param  place The place.
     *  @return       The place; nothing where it does not.
     */
    auto Kept(bool onIt, const Place &place) -> std::optional<Place>
    {
      return onIt ? std::optional<Place>{place} : std::nullopt;
    }

    /**
     *  Where on a polyline a point on its carrier lies, taken on the segment the hit's place names,
     *  so that a polyline that passes a point twice keeps the hit's pass: on a straight segment by
     *  its foot; on a bulge by the angle phi from the chord to it at the bulge's start, the chord
     *  to a point a part s of the sweep S along being the whole chord turned by (s - S) / 2, so the
     *  part is 1 + 2 phi / S; worked from the start along the chord, as the bulge's points are,
     *  never from its center, which a small bulge puts far off. A point nearer the start than 2^-20
     *  of its way to the end, as a vertex put on the other curve is, is placed by its foot on the
     *  bulge's tangent there: the way to it from the start is a few roundings of the point, its
     *  angle with the chord noise, which put such a point anywhere along the bulge; along the
     *  tangent its rounding moves the place by no more than its size, and the bulge leaves the
     *  tangent by the cube of the part, nothing so near. The start itself is such a point, its foot
     *  the start.
     *
     *  @param  polyline The polyline.
     *  @param  along    The hit's place on it.
     *  @param  point    The point.
     *  @return          The place; nothing where the point is not on the segment.
     */
    auto Strung(const Polyline &polyline, double along, Point point) -> std::optional<Place>
    {
      const auto index = static_cast<std::size_t>(
          std::clamp(std::floor(along), 0.0, static_cast<double>(Count(polyline) - 1)));
      const Vertex &start = polyline.Vertices[index];
      const Vertex &end   = polyline.Vertices[(index + 1) % polyline.Vertices.size()];
      if (Flatly(start, end))
      {
        /**
         *  The foot's place, where it is nearest by far, as Nearest would take it; else Nearest's.
         */
        const Line                 segment = {.Start = start.Position, .End = end.Position};
        const double               t    = Inside(segment.Start, segment.End - segment.Start, point);
        const std::optional<Place> foot = Footed(segment, t, point);
        const double               placed =
            static_cast<double>(index) + (foot ? foot->Along : Nearest(segment, point).Along);
        const bool within = -Slack <= t && t <= 1.0 + Slack;
        if (!within)
        {
          return std::nullopt;
        }
        /**
         *  The point at the place as At works it out: on this segment, its part along it as Locate
         *  takes it, the place less the segment's index, exactly; where the place rounds onto the
         *  next segment's start, by At.
         */
        const double part = placed - static_cast<double>(index);
        return Place{.At =
                         part < 1.0 ? Lerp(segment.Start, segment.End, part) : At(polyline, placed),
                     .Along = placed};
      }
      const Bight  bow    = Bowed(start.Position, end.Position, start.Bulge);
      const Vector toward = point - start.Position;
      const Vector rest   = end.Position - point;
      const bool   early  = Dot(toward, toward) < Dot(rest, rest);
      /**
       *  On an arc of a quarter turn or less, its bulge no more than sqrt(2) - 1, the part is told
       *  from how far along the chord the point lies, u, which is 1 / 2 + sin((part - 1 / 2) sweep)
       *  / (2 sin(sweep / 2)): the angle at the start from the chord to the point, which tells it
       *  on a longer arc, is rounded by a step of a double over itself, and on a slight arc or near
       *  the end that is as much as the arc is long.
       */
      const bool   slight = std::abs(start.Bulge) <= std::numbers::sqrt2 - 1.0;
      const double within =
          slight ? ((2.0 * Dot(toward, bow.Chord) / Dot(bow.Chord, bow.Chord)) - 1.0) * bow.Half
                 : 0.0;
      const double angled =
          slight ? 0.5 + (std::asin(std::clamp(within, -1.0, 1.0)) / bow.Sweep)
                 : 1.0 + (2.0 * Geoxy::Between(bow.Chord, toward).InRadians() / bow.Sweep);
      /**
       *  On a longer arc, by the angle from the chord at whichever end lies farther from the point,
       *  where the way to the point is told best: half the sweep left to run at the start, half the
       *  sweep run at the end.
       */
      const double part =
          early && !slight ? 2.0 * Geoxy::Between(bow.Chord, rest).InRadians() / bow.Sweep : angled;
      if (!(-Slack <= part) || !(part <= 1.0 + Slack))
      {
        return std::nullopt;
      }
      /**
       *  The point at the place as At works it out: inside the segment, from the bow worked out
       *  here already; at an end, where the place may round onto the next segment, by At.
       */
      const double placed   = static_cast<double>(index) + std::clamp(part, 0.0, 1.0);
      const double fraction = placed - static_cast<double>(index);
      const bool   inside   = fraction > 0.0 && fraction < 1.0;
      return Place{.At    = inside ? Along(bow, bow.Sweep * fraction) : At(polyline, placed),
                   .Along = placed};
    }

    /**
     *  Where on an arc a point on its circle lies, by its heading from the center, the pass the
     *  hit's place is on.
     *
     *  @param  arc   The arc.
     *  @param  along The hit's place on it.
     *  @param  point The point.
     *  @return       The place; nothing at the center or off the sweep.
     */
    auto Abeam(const Arc &arc, double along, Point point) -> std::optional<double>
    {
      if (point == arc.Center)
      {
        return std::nullopt;
      }
      return Abreast(Elliptic(arc), Heading(point - arc.Center).InRadians(), along);
    }

    /**
     *  Where on an arc or an ellipse a point on its carrier lies, by its angle: an arc's by its
     *  heading from the center, an ellipse's by its angle on the unit circle; within a rounding of
     *  the curve, the foot Nearest would find, nearest of all by far, without its search or
     *  weighing the other foot and the ends. A curve swept round more than once passes the point on
     *  every turn, and the hit keeps the turn its place names, as a polyline's keeps its segment.
     *
     *  @param  curve The curve.
     *  @param  along The hit's place on it.
     *  @param  point The point.
     *  @return       The place; nothing for other curves, off the sweep, at the center, or on an
     *                ellipse too flat to be told so, where Nearest decides.
     */
    auto Angled(const Curve &curve, double along, Point point) -> std::optional<Place>
    {
      std::optional<double> swept;
      if (const auto *arc = std::get_if<Arc>(&curve))
      {
        swept = Abeam(*arc, along, point);
      }
      if (const auto *ellipse = std::get_if<Ellipse>(&curve))
      {
        swept = Facing(*ellipse, point)
                    .and_then([ellipse, along](double angle) -> std::optional<double>
                              { return Abreast(*ellipse, angle, along); });
      }
      return swept.transform([&curve](double place) -> Place
                             { return {.At = At(curve, place), .Along = place}; });
    }

    /**
     *  What a bound below adds for the products and quotients that fall near the normal doubles'
     *  end, where a rounding is a step of the least double and no share of its size: a few such
     *  steps, many times over. Double-double arithmetic's bounds are shares of its terms' sizes,
     *  and a share of a term so small rounds to nothing while the term's own rounding does not; so
     *  a bound over such a term is no less than this, and a number it settles is no less than
     *  2^-947, its step past twice the bound.
     */
    constexpr double Underflow = 0x1p-1000;

    /**
     *  What a bound adds for its terms, products and quotients, near the subnormals: Underflow
     *  where any not zero lies below 2^-900, where a rounding below the normal doubles is no share
     *  of the bound's own; nothing otherwise, so that a term worked out exactly, as a zero along an
     *  axis, settles as it always has.
     *
     *  @param  terms The terms.
     *  @return       What the bound adds.
     */
    auto Underflowing(std::initializer_list<double> terms) -> double
    {
      for (const double term : terms)
      {
        if (!Zero(term) && std::abs(term) < 0x1p-900)
        {
          return Underflow;
        }
      }
      return 0.0;
    }

    /**
     *  The double nearest a number held in two doubles, where its error bound leaves no halfway
     *  point between doubles near it, the nearer taken as half the smaller step either side;
     *  nothing otherwise, a number or a bound past every double or not a number included, as no
     *  comparison with them holds.
     *
     *  @param  value The number.
     *  @param  bound How far the number may be from the exact value.
     *  @return       The double nearest the exact value, where certain.
     */
    auto Certain(Double value, double bound) -> std::optional<double>
    {
      const double high = value.High;
      if (Spacing(high) - (2.0 * std::abs(value.Low)) > 2.0 * bound)
      {
        return high;
      }
      return std::nullopt;
    }

    /**
     *  A difference of two doubles, exactly, in two doubles.
     *
     *  @param  head The double.
     *  @param  tail The double taken from it.
     *  @return      The difference.
     */
    auto Apart(double head, double tail) -> Double
    {
      const Two difference = TwoSum(head, -tail);
      return {.High = difference.High, .Low = difference.Low};
    }

    /**
     *  A straight curve as double-double works it: a point on it and its way, exact in two doubles,
     *  a line's from its start to its end.
     */
    struct Lane
    {
      Point                 Start = {}; /**< Where its parameter is 0.   */
      std::array<Double, 2> Way   = {}; /**< How far it goes for each 1. */
    };

    /**
     *  A line, a ray or an endless line as double-double works it; nothing for any other curve.
     *
     *  @param  curve The curve.
     *  @return       It, straight.
     */
    auto Laned(const Curve *curve) -> std::optional<Lane>
    {
      if (const auto *line = std::get_if<Line>(curve))
      {
        return Lane{.Start = line->Start,
                    .Way = {Apart(line->End.X, line->Start.X), Apart(line->End.Y, line->Start.Y)}};
      }
      if (const auto *ray = std::get_if<Ray>(curve))
      {
        return Lane{.Start = ray->Start, .Way = {Lifted(ray->Way.X), Lifted(ray->Way.Y)}};
      }
      if (const auto *xline = std::get_if<Xline>(curve))
      {
        return Lane{.Start = xline->Through, .Way = {Lifted(xline->Way.X), Lifted(xline->Way.Y)}};
      }
      return std::nullopt;
    }

    /**
     *  Where two straight curves' carriers cross, in double-double with a bound on its error: each
     *  way and the gap exact in two doubles, the cross products and their quotient each within
     *  2^-98 of their terms' sizes, 64 times what double-double makes; the point the double nearest
     *  it where the bound leaves no halfway point near, which nearly every crossing is. Nothing
     *  where the bound does not settle it, or the lines run so near parallel that it cannot, for
     *  the exact way to decide.
     *
     *  @param  first  One line.
     *  @param  second The other.
     *  @return        The double nearest the crossing, where certain.
     */
    auto Quick(const Lane &first, const Lane &second) -> std::optional<Point>
    {
      constexpr double             Grain = 0x1p-98;
      const std::array<Double, 2> &way   = first.Way;
      const std::array<Double, 2> &other = second.Way;
      const std::array<Double, 2>  gap   = {Apart(second.Start.X, first.Start.X),
                                            Apart(second.Start.Y, first.Start.Y)};
      const Double turn  = Difference(Multiply(way[0], other[1]), Multiply(way[1], other[0]));
      const Double reach = Difference(Multiply(gap[0], other[1]), Multiply(gap[1], other[0]));
      const double turnBound =
          (Grain * (std::abs(way[0].High * other[1].High) + std::abs(way[1].High * other[0].High)))
          + Underflowing({way[0].High * other[1].High, way[1].High * other[0].High});
      const double reachBound =
          (Grain * (std::abs(gap[0].High * other[1].High) + std::abs(gap[1].High * other[0].High)))
          + Underflowing({gap[0].High * other[1].High, gap[1].High * other[0].High});
      if (!(std::abs(turn.High) > 4.0 * turnBound))
      {
        return std::nullopt;
      }
      const Double along      = Divide(reach, turn);
      const double size       = std::abs(along.High);
      const double alongBound = (2.0
                                 * (reachBound + (2.0 * size * turnBound))
                                 / (std::abs(turn.High) - (2.0 * turnBound)))
                              + (Grain * size)
                              + Underflowing({along.High});
      const auto   axis = [&along, alongBound](double start, Double step) -> std::optional<double>
      {
        const Double value = Add(Lifted(start), Multiply(step, along));
        const double bound = (2.0 * std::abs(step.High) * alongBound)
                           + (Grain * (std::abs(start) + std::abs(step.High * along.High)))
                           + Underflowing({step.High * along.High});
        return Certain(value, bound);
      };
      const std::optional<double> x = axis(first.Start.X, way[0]);
      const std::optional<double> y = axis(first.Start.Y, way[1]);
      if (!x || !y)
      {
        return std::nullopt;
      }
      return Point{.X = *x, .Y = *y};
    }

    /**
     *  A double-double and how far it may be from the value it stands for.
     */
    struct Bounded
    {
      Double Value = {};  /**< The value.         */
      double Bound = 0.0; /**< How far it may be. */
    };

    /**
     *  A crossing's point in double-double, each coordinate with its bound, for the root chosen
     *  later by which point lies nearer.
     */
    struct Rooted
    {
      std::array<Bounded, 2> At = {}; /**< Each coordinate. */
    };

    /**
     *  A crossing's point rounded: the double nearest each coordinate, where its bound settles it.
     *
     *  @param  root The point in double-double, with its bounds.
     *  @return      The point; nothing where either coordinate is not certain.
     */
    auto Pointed(const Rooted &root) -> std::optional<Point>
    {
      const std::optional<double> x = Certain(root.At[0].Value, root.At[0].Bound);
      const std::optional<double> y = Certain(root.At[1].Value, root.At[1].Bound);
      if (!x || !y)
      {
        return std::nullopt;
      }
      return Point{.X = *x, .Y = *y};
    }

    /**
     *  The two points a square root's two signs give, each rounded where its bound settles it, for
     *  the nearer to be chosen by each hit.
     *
     *  @param  roots The two points, the root added first, then taken.
     *  @return       Both, each where certain.
     */
    auto Rounded(const std::array<Rooted, 2> &roots) -> Roots
    {
      return {.At = {Pointed(roots[0]), Pointed(roots[1])}, .Twice = true, .Apart = true};
    }

    /**
     *  Whether a radicand is certainly not zero: its value farther from zero than four times its
     *  bound, either way.
     *
     *  @param  radicand The radicand, with its bound.
     *  @return          Whether it is.
     */
    auto Clear(const Bounded &radicand) -> bool
    {
      return std::abs(radicand.Value.High) > 4.0 * radicand.Bound;
    }

    /**
     *  A square root in double-double, and its bound from its radicand's: nothing where the
     *  radicand is not clearly above zero, as near a tangency.
     *
     *  @param  radicand Under the root, with its bound.
     *  @return          The root, with its bound.
     */
    auto Rooting(const Bounded &radicand) -> std::optional<Bounded>
    {
      constexpr double Grain = 0x1p-98;
      if (!(radicand.Value.High > 4.0 * radicand.Bound))
      {
        return std::nullopt;
      }
      const Double root = Root(radicand.Value);
      return Bounded{.Value = root,
                     .Bound =
                         (2.0 * radicand.Bound / std::sqrt(radicand.Value.High - radicand.Bound))
                         + (Grain * std::abs(root.High))};
    }

    /**
     *  What a line's two crossings with a circle share: the quadratic's a and b, and the root of
     *  its discriminant, in double-double with their bounds.
     */
    struct Chord
    {
      Bounded                A     = {};    /**< u.u, for the line P + t u.                      */
      Bounded                B     = {};    /**< u.(P - C).                                      */
      std::optional<Bounded> Root  = {};    /**< sqrt(b^2 - a c), where clear above zero.        */
      bool                   Apart = false; /**< Whether b^2 - a c is clear of zero, either way. */
    };

    /**
     *  What a line's two crossings with a circle share: the line P + t u, the circle |X - C| = r, t
     *  = (-b +- sqrt(b^2 - a c)) / a for a = u.u, b = u.(P - C), c = |P - C|^2 - r^2.
     *
     *  @param  line The line, ray or endless line.
     *  @param  arc  The arc whose circle.
     *  @return      a, b and the root; no root where the discriminant is not clearly above zero, as
     *               near a tangency, or for a line of no length, whose a and b are zero.
     */
    auto Chorded(const Lane &line, const Arc &arc) -> Chord
    {
      constexpr double             Grain   = 0x1p-98;
      const std::array<Double, 2> &way     = line.Way;
      const std::array<Double, 2>  gap     = {Apart(line.Start.X, arc.Center.X),
                                              Apart(line.Start.Y, arc.Center.Y)};
      const Two                    squared = TwoProduct(arc.Radius, arc.Radius);
      const Double                 radius  = {.High = squared.High, .Low = squared.Low};
      const auto                   dot     = [](const std::array<Double, 2> &one,
                                                const std::array<Double, 2> &other) -> Bounded
      {
        return {.Value = Add(Multiply(one[0], other[0]), Multiply(one[1], other[1])),
                .Bound = (Grain
                          * (std::abs(one[0].High * other[0].High)
                             + std::abs(one[1].High * other[1].High)))
                       + Underflowing({one[0].High * other[0].High, one[1].High * other[1].High})};
      };
      const Bounded a    = dot(way, way);
      const Bounded b    = dot(way, gap);
      const Bounded away = dot(gap, gap);
      const Bounded c    = {.Value = Difference(away.Value, radius),
                            .Bound = away.Bound
                                   + (Grain * (std::abs(away.Value.High) + radius.High))
                                   + Underflowing({radius.High})};
      const Double  lead = Multiply(b.Value, b.Value);
      const Double  rest = Multiply(a.Value, c.Value);
      const Bounded d    = {.Value = Difference(lead, rest),
                            .Bound = 2.0
                                   * ((Grain * (std::abs(lead.High) + std::abs(rest.High)))
                                      + (2.0 * std::abs(b.Value.High) * b.Bound)
                                      + (std::abs(a.Value.High) * c.Bound)
                                      + (std::abs(c.Value.High) * a.Bound)
                                      + Underflowing({lead.High, rest.High}))};
      return {.A = a, .B = b, .Root = Rooting(d), .Apart = Clear(d)};
    }

    /**
     *  One of a line's two crossings with a circle, in double-double with its bounds: the further
     *  along the line where the root is added, the nearer where it is taken.
     *
     *  @param  line  The line, ray or endless line.
     *  @param  chord What its two crossings share, with a root.
     *  @param  sign  Which: 0 for the root added, 1 for the root taken.
     *  @return       The crossing's point.
     */
    auto Sided(const Lane &line, const Chord &chord, std::size_t sign) -> Rooted
    {
      constexpr double            Grain = 0x1p-98;
      const Bounded              &a     = chord.A;
      const Bounded              &b     = chord.B;
      const Bounded               root  = chord.Root.value_or(Bounded{});
      const std::array<double, 2> start = {line.Start.X, line.Start.Y};
      const Double                top =
          sign == 0 ? Add(Negated(b.Value), root.Value) : Difference(Negated(b.Value), root.Value);
      const double topBound =
          b.Bound + root.Bound + (Grain * (std::abs(b.Value.High) + std::abs(root.Value.High)));
      const Double along = Divide(top, a.Value);
      const double size  = std::abs(along.High);
      const double alongBound =
          (2.0 * (topBound + (size * a.Bound)) / (std::abs(a.Value.High) - a.Bound))
          + (Grain * size)
          + Underflowing({along.High});
      Rooted rooted;
      for (std::size_t axis = 0; axis < 2; ++axis)
      {
        const Double step  = line.Way.at(axis);
        rooted.At.at(axis) = {
            .Value = Add(Lifted(start.at(axis)), Multiply(step, along)),
            .Bound = (2.0 * std::abs(step.High) * alongBound)
                   + (Grain * (std::abs(start.at(axis)) + std::abs(step.High * along.High)))
                   + Underflowing({step.High * along.High})};
      }
      return rooted;
    }

    /**
     *  Where a line's carrier crosses a circle, in double-double with error bounds: the line P + t
     *  u, the circle |X - C| = r, t = (-b +- sqrt(b^2 - a c)) / a for a = u.u, b = u.(P - C), c =
     *  |P - C|^2 - r^2; both roots, each coordinate the double nearest it where its bound settles
     *  it, for each hit to take the nearer as the exact way would. Nothing where the root does not
     *  settle, for the exact way to decide.
     *
     *  @param  line The line, ray or endless line.
     *  @param  arc  The arc whose circle.
     *  @return      Both crossings, each where certain.
     */
    auto Quick(const Lane &line, const Arc &arc) -> Roots
    {
      const Chord chord = Chorded(line, arc);
      if (!chord.Root)
      {
        return {.At = {}, .Twice = false, .Apart = chord.Apart};
      }
      return Rounded({Sided(line, chord, 0), Sided(line, chord, 1)});
    }

    /**
     *  Where two circles cross, in double-double with error bounds: on their radical line, c1 + (k
     *  g +- sqrt(4 |g|^2 r1^2 - k^2) g') / (2 |g|^2), g = c2 - c1, g' it turned a quarter, k =
     *  |g|^2 + r1^2 - r2^2; both roots, each coordinate the double nearest it where its bound
     *  settles it. Nothing where the root does not settle.
     *
     *  @param  first  One arc's circle.
     *  @param  second The other's.
     *  @return        Both crossings, each where certain.
     */
    auto Quick(const Arc &first, const Arc &second) -> Roots
    {
      constexpr double            Grain = 0x1p-98;
      const std::array<Double, 2> gap   = {Apart(second.Center.X, first.Center.X),
                                           Apart(second.Center.Y, first.Center.Y)};
      const Two                   one   = TwoProduct(first.Radius, first.Radius);
      const Two                   other = TwoProduct(second.Radius, second.Radius);
      const Double                inner = {.High = one.High, .Low = one.Low};
      const Double                outer = {.High = other.High, .Low = other.Low};
      const Double                apart = Add(Multiply(gap[0], gap[0]), Multiply(gap[1], gap[1]));
      const double  spread = (Grain * std::abs(apart.High))
                           + Underflowing({gap[0].High * gap[0].High, gap[1].High * gap[1].High});
      const Double  chord  = Difference(Add(apart, inner), outer);
      const double  chordBound = spread
                               + (Grain * (std::abs(apart.High) + inner.High + outer.High))
                               + Underflowing({inner.High, outer.High});
      const Double  four       = Multiply(Multiply(Lifted(4.0), apart), inner);
      const Double  squared    = Multiply(chord, chord);
      const Bounded big        = {.Value = Difference(four, squared),
                                  .Bound = 2.0
                                         * ((Grain * (std::abs(four.High) + std::abs(squared.High)))
                                            + (4.0 * inner.High * spread)
                                            + (2.0 * std::abs(chord.High) * chordBound)
                                            + Underflowing({four.High, squared.High}))};
      /**
       *  Circles on one center have |g|^2 zero, and so no root: |g|^2 is clear of zero where there
       *  is one.
       */
      const std::optional<Bounded> root = Rooting(big);
      if (!root)
      {
        return {.At = {}, .Twice = false, .Apart = Clear(big)};
      }
      const Double                twice  = Add(apart, apart);
      const std::array<double, 2> centre = {first.Center.X, first.Center.Y};
      const std::array<Double, 2> across = {Negated(gap[1]), gap[0]};
      std::array<Rooted, 2>       roots  = {};
      for (std::size_t sign = 0; sign < 2; ++sign)
      {
        for (std::size_t axis = 0; axis < 2; ++axis)
        {
          const Double lean = Multiply(across.at(axis), root->Value);
          const Double top  = sign == 0 ? Add(Multiply(chord, gap.at(axis)), lean)
                                        : Difference(Multiply(chord, gap.at(axis)), lean);
          const double topBound =
              (std::abs(gap.at(axis).High) * chordBound)
              + (std::abs(across.at(axis).High) * root->Bound)
              + (Grain * (std::abs(chord.High * gap.at(axis).High) + std::abs(lean.High)))
              + Underflowing({chord.High * gap.at(axis).High, lean.High});
          const Double share = Divide(top, twice);
          const double size  = std::abs(share.High);
          const double shareBound =
              (2.0 * (topBound + (2.0 * size * spread)) / (std::abs(twice.High) - (2.0 * spread)))
              + (Grain * size)
              + Underflowing({share.High});
          roots.at(sign).At.at(axis) = {.Value = Add(Lifted(centre.at(axis)), share),
                                        .Bound = (2.0 * shareBound)
                                               + (Grain * (std::abs(centre.at(axis)) + size))};
        }
      }
      return Rounded(roots);
    }

    /**
     *  The crossing a hit is at, of its carriers' crossings found in double-double: the one, or of
     *  two the nearer the hit, as the exact way chooses, by the rounded points' distances, the
     *  first where they are as near; nothing where any it chooses among is not certain.
     *
     *  @param  roots The crossings.
     *  @param  near  Where the hit is.
     *  @return       The crossing, rounded.
     */
    auto Picked(const Roots &roots, Point near) -> std::optional<Point>
    {
      if (!roots.Twice)
      {
        return roots.At[0];
      }
      return roots.At[0].and_then(
          [&roots, near](Point one) -> std::optional<Point>
          {
            return roots.At[1].transform(
                [one, near](Point other) -> Point
                { return Distance(one, near) <= Distance(other, near) ? one : other; });
          });
    }

    /**
     *  Of the points two carriers cross at, rounded once, the one nearer a hit, the first where
     *  they are as near.
     *
     *  @param  crossed The points.
     *  @param  near    Where the hit is.
     *  @return         The point.
     */
    auto Chosen(const Candidate &crossed, Point near) -> Point
    {
      const auto rounded = [&crossed](int side) -> Point
      {
        return {.X = Expansion::Rounded(crossed.X, side), .Y = Expansion::Rounded(crossed.Y, side)};
      };
      Point point = rounded(1);
      if (crossed.Root)
      {
        const Point other = rounded(-1);
        point             = Distance(point, near) <= Distance(other, near) ? point : other;
      }
      return point;
    }

    /**
     *  Where two curves' carriers cross near a hit, worked out in expansions from their parts'
     *  carriers and rounded once, the value taken 2^lift times smaller first: the root nearer the
     *  hit, the hit taken as many times smaller, the first where they are as near. Each quotient's
     *  divisor is scaled up by as much as 2^900, which keeps every bit of its few sums of products
     *  of numbers below 2, and its numerator down by the rest, which drops none of numbers scaled
     *  up from below the normal doubles: so a crossing of curves scaled up from there rounds to the
     *  double nearest it at their own size.
     *
     *  @param  curves   The two curves, their numbers below 2 where lifted.
     *  @param  approach The hit's place on each, and where it is.
     *  @param  lift     The power of two the crossing is taken times smaller, not below 0.
     *  @return          The point; nothing where the carriers do not cross, or a part has none.
     */
    auto Solving(const std::array<const Curve *, 2> &curves, const Approach &approach, int lift)
        -> std::optional<Point>
    {
      const std::array<double, 2> &places = approach.Along;
      if (Unsolvable(Borne(*curves[0], places[0]), Borne(*curves[1], places[1])))
      {
        return std::nullopt;
      }
      std::optional<Candidate> crossed =
          Solved(Carried(*curves[0], places[0]), Carried(*curves[1], places[1]));
      if (!crossed || (approach.Touch && !crossed->Tangent))
      {
        return std::nullopt;
      }
      if (lift > 0)
      {
        const int    raised  = std::min(lift, 900);
        const double larger  = std::ldexp(1.0, raised);
        const double smaller = std::ldexp(1.0, raised - lift);
        for (Quotient *quotient : {&crossed->X, &crossed->Y})
        {
          quotient->Top    = Expansion::Scaled(quotient->Top, smaller);
          quotient->Bottom = Expansion::Scaled(quotient->Bottom, smaller);
          quotient->Below  = Expansion::Scaled(quotient->Below, larger);
        }
      }
      const Point near = {.X = std::ldexp(approach.At.X, -lift),
                          .Y = std::ldexp(approach.At.Y, -lift)};
      return Chosen(*crossed, near);
    }

    /**
     *  Where two polylines' bulges cross near a hit, from Crossbows: of their circles' two
     *  crossings the one nearer the hit, as the exact way chooses it, where the two lie clearly
     *  apart in distance from it, the rounded points' distances then ordered alike; the double
     *  nearest it where its bounds settle it.
     *
     *  @param  bellies The two bulges, the smaller circle's first, as Narrower tells.
     *  @param  near    Where the hit is.
     *  @return         The point; nothing where it is not settled so, for the exact way.
     */
    auto Bowline(const std::array<Belly, 2> &bellies, Point near) -> std::optional<Point>
    {
      const std::optional<Crossbow> bows = Crossbows(bellies);
      if (!bows)
      {
        return std::nullopt;
      }
      double                slack = 0.0;
      std::array<double, 2> away  = {};
      for (std::size_t sign = 0; sign < 2; ++sign)
      {
        const std::array<Sure::Number, 2> &point = bows->At.at(sign);
        away.at(sign) = Distance({.X = point[0].Value.High, .Y = point[1].Value.High}, near);
        for (const Sure::Number &coordinate : point)
        {
          slack += 4.0 * (coordinate.Bound + (0x1p-50 * std::abs(coordinate.Value.High)));
        }
      }
      slack += 0x1p-30 * std::max(away[0], away[1]);
      if (!(std::abs(away[0] - away[1]) > slack))
      {
        return std::nullopt;
      }
      return Sure::Rounded(bows->At.at(away[0] < away[1] ? 0 : 1))
          .transform([](const std::array<double, 2> &point) -> Point
                     { return {.X = point[0], .Y = point[1]}; });
    }

    /**
     *  Where two curves cross near a hit, exactly, as Crossing works it out: the point alone, not
     *  yet placed on the curves.
     *
     *  @param  curves   The two curves.
     *  @param  approach The hit's place on each, and where it is.
     *  @param  roots    The curves' crossings found in double-double, as Quickly finds them.
     *  @return          The point; nothing where they have no crossing.
     */
    auto Spotting(const std::array<const Curve *, 2> &curves, const Approach &approach,
                  const Roots &roots) -> std::optional<Point>
    {
      if (approach.Touch && roots.Apart)
      {
        return std::nullopt;
      }
      const std::array<double, 2> &places = approach.Along;
      const Point                  near   = approach.At;
      if (const std::optional<Point> quick = approach.Touch ? std::nullopt : Picked(roots, near))
      {
        return quick;
      }
      if (!approach.Touch)
      {
        /**
         *  A crossing the double-double roots left open: a line's or a circle's own short sums
         *  decide it, as nearly always, the same point the exact crossing would round to.
         */
        const Touching shortly = Shortly(curves, near);
        if (shortly.Decided)
        {
          return shortly.At;
        }
      }
      if (approach.Touch)
      {
        /**
         *  A touch asked for: the circles' own short sums decide it, as nearly always, the same
         *  point the exact crossing would round to.
         */
        const Touching touching = Tangency(curves);
        if (touching.Decided)
        {
          return touching.At;
        }
      }
      const std::optional<Exact::Straightway> ruled = Ruled(*curves[0], places[0]);
      const std::optional<Exact::Straightway> ruling =
          ruled ? Ruled(*curves[1], places[1]) : std::nullopt;
      if (ruling && !Holdable({*ruled, *ruling}))
      {
        /**
         *  Two lines whose numbers part too far for sums of doubles to hold their products, as a
         *  coordinate a sliver beside the rest: worked out in whole numbers, which hold them all.
         *  Two lines cross where they are not parallel and never touch, as the sums find too.
         */
        return approach.Touch ? std::nullopt : Exact::CrossWhole({*ruled, *ruling});
      }
      return Solving(curves, approach, 0);
    }
  } // namespace

  auto Located(const Curve &curve, double along, Point point) -> std::optional<Place>
  {
    if (const auto *polyline = std::get_if<Polyline>(&curve))
    {
      return Strung(*polyline, along, point);
    }
    if (const std::optional<Place> angled = Angled(curve, along, point))
    {
      return angled;
    }
    if (const auto *line = std::get_if<Line>(&curve))
    {
      const double t      = Inside(line->Start, line->End - line->Start, point);
      const bool   within = -Slack <= t && t <= 1.0 + Slack;
      if (!within)
      {
        return std::nullopt;
      }
      if (const std::optional<Place> foot = Footed(*line, t, point))
      {
        return foot;
      }
      return Nearest(*line, point);
    }
    const Place place = Nearest(curve, point);
    if (const auto *ray = std::get_if<Ray>(&curve))
    {
      return Kept(Inside(ray->Start, ray->Way, point) >= -Slack, place);
    }
    if (std::holds_alternative<Xline>(curve))
    {
      return place;
    }
    const Range span = Domain(curve);
    if (span.Start < place.Along && place.Along < span.End)
    {
      return place;
    }
    /**
     *  A few steps of a double at the size of the curve's numbers and the point's, at any scale:
     *  held to steps of a double at 1, a point past a small curve's end by many of its own steps
     *  was taken for on it, an end hit put at a crossing far past the end, and at a scale 2^k times
     *  larger not.
     */
    const double size = std::max({std::abs(point.X), std::abs(point.Y), Met::Size(curve)});
    return Kept(Distance(place.At, point) <= 8.0 * size * std::numeric_limits<double>::epsilon(),
                place);
  }

  auto Placed(const std::array<const Curve *, 2> &curves, const std::array<double, 2> &places,
              Point point) -> std::optional<Crossed>
  {
    const std::optional<Place> one   = Located(*curves[0], places[0], point);
    const std::optional<Place> other = Located(*curves[1], places[1], point);
    if (!one || !other)
    {
      return std::nullopt;
    }
    return Crossed{.At = point, .Along = {one->Along, other->Along}, .On = {one->At, other->At}};
  }

  auto Quickly(const Line &line, const Arc &arc, bool further) -> std::optional<Point>
  {
    const Lane  lane  = {.Start = line.Start,
                         .Way = {Apart(line.End.X, line.Start.X), Apart(line.End.Y, line.Start.Y)}};
    const Chord chord = Chorded(lane, arc);
    if (!chord.Root)
    {
      return std::nullopt;
    }
    return Pointed(Sided(lane, chord, further ? 0 : 1));
  }

  auto Quickly(const std::array<const Curve *, 2> &curves) -> Roots
  {
    const std::optional<Lane> one   = Laned(curves[0]);
    const std::optional<Lane> other = one ? Laned(curves[1]) : std::nullopt;
    if (one && other)
    {
      return {.At = {Quick(*one, *other), std::nullopt}, .Twice = false, .Apart = false};
    }
    const auto *round = std::get_if<Arc>(curves[0]);
    const auto *bent  = std::get_if<Arc>(curves[1]);
    if (one && bent != nullptr)
    {
      return Quick(*one, *bent);
    }
    if (round != nullptr)
    {
      if (const std::optional<Lane> across = Laned(curves[1]))
      {
        return Quick(*across, *round);
      }
    }
    if (round != nullptr && bent != nullptr)
    {
      return Quick(*round, *bent);
    }
    return {};
  }

  auto Crossing(const std::array<const Curve *, 2> &curves, const Approach &approach,
                const Roots &roots) -> std::optional<Crossed>
  {
    return Spotting(curves, approach, roots)
        .and_then([&curves, &approach](Point point) -> std::optional<Crossed>
                  { return Placed(curves, approach.Along, point); });
  }

  auto Aiming(const std::array<const Curve *, 2> &curves, const Approach &approach,
              const Roots &roots) -> std::optional<Bearings>
  {
    const std::optional<Point> point = Spotting(curves, approach, roots);
    if (!point)
    {
      return std::nullopt;
    }
    Bearings found = {.At = *point, .Along = {}, .Ways = {}};
    for (std::size_t side = 0; side < 2; ++side)
    {
      const Curve &curve = *curves.at(side);
      const double along = approach.Along.at(side);
      const auto  *arc   = std::get_if<Arc>(&curve);
      /**
       *  On an arc placed by its angle, as Angled places it, the direction At would turn to there:
       *  Located's place and At's way, At's point left out.
       */
      const std::optional<double> swept =
          arc != nullptr ? Abeam(*arc, along, *point) : std::nullopt;
      if (swept)
      {
        found.Along.at(side) = *swept;
        found.Ways.at(side)  = Direction(arc->Start + (arc->Sweep * *swept));
        continue;
      }
      const std::optional<Place> place = Located(curve, along, *point);
      if (!place)
      {
        return std::nullopt;
      }
      found.Along.at(side) = place->Along;
    }
    return found;
  }

  auto Spotted(const std::array<const Curve *, 2> &curves, const Approach &approach, int lift)
      -> std::optional<Point>
  {
    if (lift > 0)
    {
      return Solving(curves, approach, lift);
    }
    return Spotting(curves, approach, Quickly(curves));
  }

  auto Crossing(const std::array<const Curve *, 2> &curves, const std::array<double, 2> &places,
                Point near) -> std::optional<Crossed>
  {
    return Crossing(curves, {.Along = places, .At = near}, Quickly(curves));
  }

  auto Bowed(const std::array<const Polyline *, 2> &polylines, const std::array<double, 2> &places,
             Point near) -> std::optional<Point>
  {
    std::array<const Polyline *, 2> ordered = polylines;
    std::array<double, 2>           along   = places;
    std::array<Belly, 2>            bellies = {};
    for (std::size_t side = 0; side < 2; ++side)
    {
      const std::array<const Vertex *, 2> seat = Seated(*polylines.at(side), places.at(side));
      bellies.at(side)                         = {
          .Start = seat[0]->Position, .End = seat[1]->Position, .Bulge = seat[0]->Bulge};
    }
    if (Narrower(bellies))
    {
      std::swap(bellies[0], bellies[1]);
      std::swap(ordered[0], ordered[1]);
      std::swap(along[0], along[1]);
    }
    if (const std::optional<Point> hop = Hopped(bellies, near))
    {
      return hop;
    }
    if (const std::optional<Point> swift = Bowline(bellies, near))
    {
      return swift;
    }
    const auto bulge = [](const Polyline &polyline, double place) -> Bulge
    {
      return std::get<Bulge>(Carried(polyline, place));
    };
    /**
     *  A coordinate of nothing is given the sign the polish's sums leave it, as Newton's method,
     *  settling on it, gave it: x + 0 is +0 for either zero, and x for any other x.
     */
    return Solved(bulge(*ordered[0], along[0]), bulge(*ordered[1], along[1]))
        .transform(
            [near](const Candidate &crossed) -> Point
            {
              const Point point = Chosen(crossed, near);
              return {.X = point.X + 0.0, .Y = point.Y + 0.0};
            });
  }
} // namespace Geoxy::Crossings
