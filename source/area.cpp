/***************************************************************************************************
 *  @file      source/area.cpp
 *  @brief     The area a curve encloses and its centroid, by Green's theorem: straight segments by
 *             a compensated shoelace, ellipses in closed form, polynomial pieces by their Bernstein
 *             forms, rational ones by Gauss-Legendre rules whose error is bounded.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/area.hpp"

#include "arctangent.hpp"
#include "basis.hpp"
#include "eased.hpp"
#include "elliptic.hpp"
#include "exact/two.hpp"
#include "kernel/piece.hpp"
#include "lean.hpp"
#include "rational.hpp"
#include "real.hpp"
#include "sure.hpp"
#include "trig.hpp"
#include "twofold.hpp"
#include "walk.hpp"

#include "geoxy/accuracy.hpp"
#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
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
#include <optional>
#include <ranges>
#include <span>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    using Kernel::Enclosure;

    /**
     *  How many terms of a series are summed where it is.
     */
    constexpr std::size_t SeriesTerms = 24;

    /**
     *  The series of s - sin s over s^3: 1 / 3! - s^2 / 5! + s^4 / 7! - ..., its coefficients by
     *  the power of s^2.
     */
    constexpr std::array<double, SeriesTerms> Excesses = [] -> std::array<double, SeriesTerms>
    {
      std::array<double, SeriesTerms> terms{};
      double                          factorial = 6.0;
      for (std::size_t index = 0; index < SeriesTerms; ++index)
      {
        const auto odd   = static_cast<double>((2 * index) + 3);
        terms.at(index)  = (index % 2 == 0 ? 1.0 : -1.0) / factorial;
        factorial       *= (odd + 1.0) * (odd + 2.0);
      }
      return terms;
    }();

    /**
     *  The series of what a bulge b encloses with a chord of length 1, over b: 1 / 3 + b^2 / 15 -
     *  b^4 / 105 + ..., the coefficient of b^(2 k - 2) being (-1)^k / ((2 k + 1) (2 k - 1) (2 k -
     *  3)), k from 1 on.
     */
    constexpr std::array<double, SeriesTerms> Lunes = [] -> std::array<double, SeriesTerms>
    {
      std::array<double, SeriesTerms> terms{};
      for (std::size_t index = 0; index < SeriesTerms; ++index)
      {
        const auto order = static_cast<double>(index + 1);
        terms.at(index)  = (index % 2 == 0 ? -1.0 : 1.0)
                         / (((2.0 * order) + 1.0) * ((2.0 * order) - 1.0) * ((2.0 * order) - 3.0));
      }
      return terms;
    }();

    /**
     *  The series of the moment about the chord's middle of what a bulge b encloses with a chord of
     *  length 1, over b^2: 1 / 15 + 2 b^2 / 105 - ..., the coefficient of b^(2 k - 4) being (-1)^(k
     *  + 1) (k - 1) / ((2 k - 5) (2 k - 3) (2 k - 1) (2 k + 1)), k from 2 on.
     */
    constexpr std::array<double, SeriesTerms> Pulls = [] -> std::array<double, SeriesTerms>
    {
      std::array<double, SeriesTerms> terms{};
      for (std::size_t index = 0; index < SeriesTerms; ++index)
      {
        const auto order = static_cast<double>(index + 2);
        terms.at(index)  = (index % 2 == 0 ? -1.0 : 1.0)
                         * (order - 1.0)
                         / (((2.0 * order) - 5.0)
                            * ((2.0 * order) - 3.0)
                            * ((2.0 * order) - 1.0)
                            * ((2.0 * order) + 1.0));
      }
      return terms;
    }();

    /**
     *  A power series at a value, by Horner's rule.
     *
     *  @param  terms Its coefficients, the constant first.
     *  @param  value Where.
     *  @return       Its sum there.
     */
    auto Summed(std::span<const double> terms, double value) -> double
    {
      double sum = 0.0;
      for (const double term : std::views::reverse(terms))
      {
        sum = term + (value * sum);
      }
      return sum;
    }

    /**
     *  The excess of a sweep over its sine, s - sin s: twice what the unit circle's arc of that
     *  sweep encloses with its chord, every turn counted, in double-double.
     *
     *  Below 2 in magnitude the two nearly cancel, and it is s^3 times its series, whose 24 terms
     *  reach far below a step of a double of it: s^3 and the leading 1 / 6 held in double-double,
     *  the rest of the series, a fifth of it at most, in doubles. From 2 on the sine is off by a
     *  step of a double of 1 at most, against an excess of 1.09 at least, and the difference is
     *  held exactly.
     *
     *  @param  sweep The sweep, in radians.
     *  @return       Its excess.
     */
    auto Excess(double sweep) -> Twofold::Double
    {
      if (!(std::abs(sweep) < 2.0))
      {
        const Twofold::Two excess = Twofold::TwoSum(sweep, -Trig::Sin(sweep));
        return {.High = excess.High, .Low = excess.Low};
      }
      const Twofold::Two    square = Twofold::TwoProduct(sweep, sweep);
      const Twofold::Double cube =
          Twofold::Multiply({.High = square.High, .Low = square.Low}, Twofold::Lifted(sweep));
      const double rest = square.High * Summed(std::span{Excesses}.subspan(1), square.High);
      const Twofold::Double series = Twofold::Add(
          Twofold::Divide(Twofold::Lifted(1.0), Twofold::Lifted(6.0)), Twofold::Lifted(rest));
      return Twofold::Multiply(cube, series);
    }

    /**
     *  A map's determinant in double-double, its two products exact.
     *
     *  @param  map The map.
     *  @return     The determinant: how it scales areas, negative where it mirrors them.
     */
    auto Determined(const Matrix &map) -> Twofold::Double
    {
      const Twofold::Two first  = Twofold::TwoProduct(map.X.X, map.Y.Y);
      const Twofold::Two second = Twofold::TwoProduct(map.X.Y, map.Y.X);
      return Twofold::Difference({.High = first.High, .Low = first.Low},
                                 {.High = second.High, .Low = second.Low});
    }

    /**
     *  What an arc of an ellipse encloses with its chord: det M (s - sin s) / 2, for the map M and
     *  the sweep s, whatever the ellipse's place, rounded once. It is the area of the unit circle's
     *  arc with its chord, (s - sin s) / 2 about the arc's start, as Green's theorem gives it,
     *  scaled by M.
     *
     *  @param  ellipse The arc.
     *  @return         The signed area, every turn counted.
     */
    auto Lens(const Ellipse &ellipse) -> double
    {
      const Twofold::Double area =
          Twofold::Multiply(Determined(ellipse.Map), Excess(ellipse.Sweep.InRadians()));
      return (area.High + area.Low) / 2.0;
    }

    /**
     *  The centroid of what an arc of an ellipse encloses with its chord, in closed form: on the
     *  unit circle the moment about the center is 2 sin^3(s / 2) / 3 along the bisector of the
     *  sweep, and the area (s - sin s) / 2, so the centroid lies 4 sin^3(s / 2) / (3 (s - sin s))
     *  along it; the map takes it to the ellipse's, the determinant dropping out.
     *
     *  @param  ellipse The arc.
     *  @return         The centroid; none for no area.
     */
    auto Middle(const Ellipse &ellipse) -> std::optional<Point>
    {
      const double          sweep  = ellipse.Sweep.InRadians();
      const Twofold::Double turned = Determined(ellipse.Map);
      const double          excess = Excess(sweep).High;
      if (Zero(turned.High) || Zero(excess))
      {
        return std::nullopt;
      }
      const double sine  = Trig::Sin(sweep / 2.0);
      const double reach = 4.0 * sine * sine * sine / (3.0 * excess);
      return ellipse.Map * (Point{} + (Direction(ellipse.Start + (ellipse.Sweep / 2.0)) * reach));
    }

    /**
     *  What a bulged segment encloses with its chord: the area, and its moment about the chord's
     *  middle, from the chord and the bulge alone, with no center, which a small bulge puts far
     *  off.
     *
     *  With c the chord's length, the sweep is 4 atan b and the radius c (1 + b^2) / (4 |b|), so
     *  the area is c^2 (b + 1 / b)^2 (s - sin s) / 32. Its moment points right of the chord, the
     *  way the arc bulges for a positive bulge and the area's sign turns for a negative one, and is
     *  c^2 |b + 1 / b|^3 F(|s|) / 128 times the chord turned right, F(s) being 4 sin^3(s / 2) / 3
     *  less (s less sin s) cos(s / 2). For |b| up to 1/2 both are taken from their series in b
     *  instead, which take no arctangent, its rounding times 3 in the area of a shallow arc, and
     *  which F's own cancelling never reaches.
     *
     *  @param  chord The chord, from the segment's start to its end.
     *  @param  bulge The bulge, not 0.
     *  @return       The signed area, and its moment about the chord's middle.
     */
    auto Lune(Vector chord, double bulge) -> Enclosure
    {
      const double square = Dot(chord, chord);
      const Vector right  = {.X = chord.Y, .Y = -chord.X};
      const double shrunk = bulge * bulge;
      if (!(std::abs(bulge) > 0.5))
      {
        return {.Area   = square * bulge * Summed(Lunes, shrunk),
                .Moment = right * (square * shrunk * Summed(Pulls, shrunk))};
      }
      const double sweep  = 4.0 * Trig::Atan(bulge);
      const double spread = std::abs(bulge + (1.0 / bulge));
      const double excess = Excess(sweep).High;
      const double half   = std::abs(sweep) / 2.0;
      const double sine   = Trig::Sin(half);
      const double pull   = (4.0 * sine * sine * sine / 3.0) - (std::abs(excess) * Trig::Cos(half));
      return {.Area   = square * spread * spread * excess / 32.0,
              .Moment = right * (square * spread * spread * spread * pull / 128.0)};
    }

    /**
     *  What a curve encloses, gathered about a point near it, the area carried with its rounding
     *  error.
     */
    class Gathering
    {
    public:
      /**
       *  Starts at nothing, about a point.
       *
       *  @param origin  The point areas and moments are taken about.
       *  @param moments Whether moments are gathered too, or the area alone.
       */
      explicit Gathering(Point origin, bool moments = true) : m_Origin(origin), m_Moments(moments)
      {
      }

      /**
       *  Adds what a stretch encloses.
       *
       *  @param enclosure The stretch's area and moment, about the origin.
       */
      void Add(const Enclosure &enclosure)
      {
        Carry(enclosure.Area);
        m_Moment = m_Moment + enclosure.Moment;
      }

      /**
       *  Adds what a straight segment encloses: half the cross product of its ends, and a sixth of
       *  that times their sum, the products carried exactly.
       *
       *  @param start Where it starts.
       *  @param end   Where it ends.
       */
      void Segment(Point start, Point end)
      {
        const Vector     p      = start - m_Origin;
        const Vector     q      = end - m_Origin;
        const Exact::Two first  = Exact::TwoProduct(p.X, q.Y);
        const Exact::Two second = Exact::TwoProduct(p.Y, q.X);
        Carry(first.Big / 2.0);
        Carry(-second.Big / 2.0);
        m_Error += (first.Small - second.Small) / 2.0;
        if (m_Moments)
        {
          const double cross = (first.Big - second.Big) + (first.Small - second.Small);
          m_Moment           = m_Moment + ((p + q) * (cross / 6.0));
        }
      }

      /**
       *  Adds what a bulged segment encloses: its chord's share as a straight segment's, and what
       *  the arc encloses with the chord from the bulge itself, its moment taken about the chord's
       *  middle and moved to the origin.
       *
       *  @param start Where it starts.
       *  @param end   Where it ends.
       *  @param bulge The bulge, not 0.
       */
      void Bulged(Point start, Point end, double bulge)
      {
        Segment(start, end);
        const Vector    chord = end - start;
        const Enclosure lune  = Lune(chord, bulge);
        Carry(lune.Area);
        if (m_Moments)
        {
          m_Moment = m_Moment + lune.Moment + ((Midpoint(start, end) - m_Origin) * lune.Area);
        }
      }

      /**
       *  What has been gathered.
       *
       *  @return The area and the moment, about the origin.
       */
      [[nodiscard]] auto Total() const -> Enclosure
      {
        return {.Area = m_Area + m_Error, .Moment = m_Moment};
      }

      /**
       *  The point areas and moments are taken about.
       *
       *  @return The origin.
       */
      [[nodiscard]] auto Origin() const -> Point { return m_Origin; }

    private:
      /**
       *  Adds to the area, its rounding error carried apart.
       *
       *  @param value What to add.
       */
      void Carry(double value)
      {
        const Exact::Two sum  = Exact::TwoSum(m_Area, value);
        m_Area                = sum.Big;
        m_Error              += sum.Small;
      }

      Point  m_Origin  = {};   /**< The point areas and moments are taken about. */
      bool   m_Moments = true; /**< Whether moments are gathered.                */
      double m_Area    = 0.0;  /**< The area so far, rounded.                    */
      double m_Error   = 0.0;  /**< The rounding errors of the area so far.      */
      Vector m_Moment  = {};   /**< The moment so far.                           */
    };

    /**
     *  The centroid of what has been gathered.
     *
     *  @param  gathering What has been gathered.
     *  @return           The origin plus the moment over the area; none for no area.
     */
    auto Centered(const Gathering &gathering) -> std::optional<Point>
    {
      const Enclosure total = gathering.Total();
      if (Zero(total.Area))
      {
        return std::nullopt;
      }
      return gathering.Origin() + (total.Moment / total.Area);
    }

    /**
     *  What a polyline encloses, closed when it is open.
     *
     *  @param  polyline The polyline.
     *  @param  moments  Whether moments are gathered too.
     *  @return          The gathering.
     */
    auto Encompassed(const Polyline &polyline, bool moments) -> Gathering
    {
      const std::size_t count = Count(polyline);
      if (count == 0)
      {
        return Gathering{{}};
      }
      Gathering gathering{polyline.Vertices.front().Position, moments};
      for (std::size_t index = 0; index < count; ++index)
      {
        const Step step = Locate(polyline, static_cast<double>(index));
        if (Straightness(step))
        {
          gathering.Segment(step.Start, step.End);
        }
        else
        {
          gathering.Bulged(step.Start, step.End, step.Bulge);
        }
      }
      if (!polyline.Closed)
      {
        gathering.Segment(polyline.Vertices.back().Position, polyline.Vertices.front().Position);
      }
      return gathering;
    }

    /**
     *  What a bulged segment encloses with its chord, in double-double: c^2 b S(b^2) for |b| under
     *  2^-6, S(x) = 1 / 3 + x / 15 - x^2 / 105 + ... Lunes's series, its first term exact and the
     *  rest in doubles, past 2^-81 of it at its sixth; else c^2 ((1 + b^2)^2 atan b - b (1 - b^2))
     *  / (8 b^2), the circle's segment R^2 (s - sin s) / 2 for R = c (1 + b^2) / (4 |b|), s = 4
     *  atan b and sin s = 4 b (1 - b^2) / (1 + b^2)^2, with no center and no sine.
     *
     *  The two terms of the numerator cancel to 8 b^3 / 3 of b at worst, some 2^11 times at 2^-6,
     *  so its arc tangent's 2^-68 (ArcTangent) leaves it within 2^-57 of itself; every other step
     *  is within 2^-100 (lean.hpp), the chord's square exact but for that.
     *
     *  @param  chord The chord, from the segment's start to its end, exactly.
     *  @param  bulge The bulge, not under Straight in magnitude.
     *  @return       The signed area.
     */
    auto Lunar(const Lean::Way &chord, double bulge) -> Twofold::Double
    {
      constexpr Twofold::Double Third   = {.High = 0x1.5555555555555p-2,
                                           .Low  = 0x1.5555555555555p-56};
      const Twofold::Double     one     = {.High = 1.0, .Low = 0.0};
      const Twofold::Double     square  = Lean::Dot(chord, chord).Value;
      const Twofold::Two        product = Twofold::TwoProduct(bulge, bulge);
      const Twofold::Double     shrunk  = {.High = product.High, .Low = product.Low};
      const Twofold::Double     height  = {.High = bulge, .Low = 0.0};
      if (std::abs(bulge) < 0x1p-6)
      {
        const double x = shrunk.High;
        const double tail =
            x
            * ((1.0 / 15.0)
               + (x
                  * ((-1.0 / 105.0)
                     + (x * ((1.0 / 315.0) + (x * ((-1.0 / 693.0) + (x / 1287.0))))))));
        return Lean::Times(square,
                           Lean::Times(height, Lean::Plus(Third, {.High = tail, .Low = 0.0})));
      }
      const Twofold::Double spread  = Lean::Plus(one, shrunk);
      const Twofold::Two    tangent = ArcTangent(bulge);
      const Twofold::Double turned =
          Lean::Times(Lean::Times(spread, spread), {.High = tangent.High, .Low = tangent.Low});
      const Twofold::Double chorded   = Lean::Times(height, Lean::Plus(one, Lean::Minus(shrunk)));
      const Twofold::Double numerator = Lean::Plus(turned, Lean::Minus(chorded));
      const Twofold::Double eight     = {.High = 8.0 * shrunk.High, .Low = 8.0 * shrunk.Low};
      const Twofold::Double quotient  = Lean::Over(numerator, eight, 1.0 / eight.High);
      const Twofold::Two    ratio     = Twofold::FastTwoSum(quotient.High, quotient.Low);
      return Lean::Times(square, {.High = ratio.High, .Low = ratio.Low});
    }

    /**
     *  What a polyline encloses, closed when it is open, in double-double and rounded once: each
     *  segment's half cross product of its ends less its first vertex, exactly as those differences
     *  are (Lean::Apart), within 2^-100 of its terms (Lean::Cross), and each bulge's segment of its
     *  circle (Lunar); the high parts summed by TwoSum, the rest apart. Where a difference lies
     *  past 2^450 or a nonzero one under 2^-450, whose products would leave the doubles, as
     *  Encompassed gathers it.
     *
     *  @param  polyline The polyline.
     *  @return          The signed area.
     */
    auto Polygonal(const Polyline &polyline) -> double
    {
      const std::size_t count = Count(polyline);
      if (count == 0)
      {
        return 0.0;
      }
      const std::span<const Vertex> vertices = polyline.Vertices;
      const Point                   origin   = vertices.front().Position;
      double                        sum      = 0.0;
      double                        rest     = 0.0;
      Lean::Way                     behind   = {};
      for (std::size_t index = 0; index < count; ++index)
      {
        const std::size_t next  = index + 1 < vertices.size() ? index + 1 : 0;
        const Point       start = vertices[index].Position;
        const Point       end   = vertices[next].Position;
        const Lean::Way   ahead = Lean::Apart(end, origin);
        const double      most  = std::max(std::abs(ahead[0].High), std::abs(ahead[1].High));
        if (most > 0x1p450 || (most > 0.0 && most < 0x1p-450) || !std::isfinite(most))
        {
          return Encompassed(polyline, false).Total().Area;
        }
        Twofold::Double piece = Lean::Cross(behind, ahead).Value;
        piece                 = {.High = piece.High / 2.0, .Low = piece.Low / 2.0};
        const Step step       = {
            .Start = start, .End = end, .Bulge = vertices[index].Bulge, .Fraction = 0.0};
        if (!Straightness(step))
        {
          piece = Lean::Plus(piece, Lunar(Lean::Apart(end, start), step.Bulge));
        }
        const Twofold::Two added  = Twofold::TwoSum(sum, piece.High);
        sum                       = added.High;
        rest                     += added.Low + piece.Low;
        behind                    = ahead;
      }
      return sum + rest;
    }

    /**
     *  What a polynomial piece encloses with its chord.
     *
     *  @param  piece The piece.
     *  @return       The gathering.
     */
    auto Encompassed(const Kernel::Piece &piece) -> Gathering
    {
      const Point start = Kernel::Start(piece);
      Gathering   gathering{start};
      gathering.Add(Kernel::Enclosed(piece, start));
      gathering.Segment(Kernel::End(piece), start);
      return gathering;
    }

    /**
     *  What a spline encloses with its chord: one of degree 1 as the polygon of its poles, which it
     *  runs straight between whatever its weights; a polynomial one piece by piece, in closed form;
     *  a rational one strip by strip, each strip's share of the accuracy its share of the domain;
     *  its moments only where they are asked for.
     *
     *  @param  spline   The spline, valid.
     *  @param  accuracy How far the area and each moment may be off, for a rational spline.
     *  @param  moments  Whether moments are gathered too.
     *  @return          The gathering.
     */
    auto Encompassed(const Spline &spline, Accuracy accuracy, bool moments) -> Gathering
    {
      if (spline.Degree == 1)
      {
        Gathering gathering{spline.Poles.front(), moments};
        for (std::size_t index = 0; index + 1 < spline.Poles.size(); ++index)
        {
          gathering.Segment(spline.Poles[index], spline.Poles[index + 1]);
        }
        gathering.Segment(spline.Poles.back(), spline.Poles.front());
        return gathering;
      }
      const Range domain = Domain(spline);
      const Point start  = At(spline, domain.Start);
      Gathering   gathering{start, moments};
      if (!Rational(spline))
      {
        Kernel::Sections(spline, [&gathering, start](const Kernel::Piece &piece, Range) -> void
                         { gathering.Add(Kernel::Enclosed(piece, start)); });
      }
      else
      {
        for (const Strip &strip : Stripped(spline))
        {
          const double share = accuracy.Error
                             * strip.Where.Width
                             * ((strip.Where.To - strip.Where.From) / (domain.End - domain.Start));
          gathering.Add(Enclosed(strip, {.Origin = start, .Error = share, .Moments = moments}));
        }
      }
      gathering.Segment(At(spline, domain.End), start);
      return gathering;
    }

    /**
     *  Which way of a spline's two its measures are taken along: the spline itself, or its reverse,
     *  its poles and weights turned about and its knots turned about and negated, the same curve
     *  run the other way. Of the two, the one whose poles, then weights, then knots come first,
     *  read in order, is taken, so that a spline and its reverse are measured along one way and
     *  their areas are negated exactly.
     *
     *  @param  spline The spline.
     *  @return        1 where the spline is taken, -1 where its reverse is, and 0 where the two are
     *                 the same, the curve running out and back along itself.
     */
    auto Sense(const Spline &spline) -> int
    {
      const auto order = [](double ahead, double behind) -> int
      {
        return static_cast<int>(ahead < behind) - static_cast<int>(behind < ahead);
      };
      int               way   = 0;
      const std::size_t count = spline.Poles.size();
      for (std::size_t index = 0; index < count && way == 0; ++index)
      {
        const Point ahead  = spline.Poles[index];
        const Point behind = spline.Poles[count - 1 - index];
        way                = order(ahead.X, behind.X);
        way                = way != 0 ? way : order(ahead.Y, behind.Y);
      }
      const std::size_t weights = spline.Weights.size();
      for (std::size_t index = 0; index < weights && way == 0; ++index)
      {
        way = order(spline.Weights[index], spline.Weights[weights - 1 - index]);
      }
      const std::size_t knots = spline.Knots.size();
      for (std::size_t index = 0; index < knots && way == 0; ++index)
      {
        way = order(spline.Knots[index], -spline.Knots[knots - 1 - index]);
      }
      return way;
    }

    /**
     *  A spline's reverse: the same curve run the other way, its parameter negated.
     *
     *  @param  spline The spline.
     *  @return        Its reverse.
     */
    auto Rewound(Spline spline) -> Spline
    {
      std::ranges::reverse(spline.Poles);
      std::ranges::reverse(spline.Weights);
      std::ranges::reverse(spline.Knots);
      for (double &knot : spline.Knots)
      {
        knot = -knot;
      }
      return spline;
    }

    /**
     *  A point's offset from another, each coordinate held exactly in two doubles.
     */
    struct Displacement
    {
      Twofold::Double X = {}; /**< Along x. */
      Twofold::Double Y = {}; /**< Along y. */
    };

    /**
     *  A point's offset from another, exactly, then times a power of two.
     *
     *  @param  point  The point.
     *  @param  origin The other.
     *  @param  scale  The power of two.
     *  @return        The offset, scaled.
     */
    auto Translated(Point point, Point origin, double scale) -> Displacement
    {
      const Twofold::Two x = Twofold::TwoSum(point.X, -origin.X);
      const Twofold::Two y = Twofold::TwoSum(point.Y, -origin.Y);
      return {.X = {.High = x.High * scale, .Low = x.Low * scale},
              .Y = {.High = y.High * scale, .Low = y.Low * scale}};
    }

    /**
     *  The cross product of two offsets, in double-double.
     *
     *  @param  lhs One offset.
     *  @param  rhs The other.
     *  @return     Their cross product.
     */
    auto Wedged(const Displacement &lhs, const Displacement &rhs) -> Twofold::Double
    {
      return Twofold::Difference(Twofold::Multiply(lhs.X, rhs.Y), Twofold::Multiply(lhs.Y, rhs.X));
    }

    /**
     *  What a polynomial Bézier curve of degree 2 or 3 encloses with its chord, in closed form
     *  about its start, p_i its control points less the start: p1 x p2 / 3 for a quadratic, and 3
     *  (p1 x p2 + p1 x p3 + 2 p2 x p3) / 20 for a cubic.
     *
     *  Each offset is held exactly in two doubles, after a power of two that brings the largest
     *  near 1 so that no product overflows, and the products and their sum in double-double:
     *  rounded once, the area is off by a step of a double but where the products cancel to far
     *  below themselves.
     *
     *  @param  points The control points, three or four.
     *  @return        The signed area.
     */
    auto Chordal(std::span<const Point> points) -> double
    {
      double largest = 0.0;
      for (const Point point : points)
      {
        largest = std::max(largest, Largest(point - points.front()));
      }
      if (!(largest > 0.0))
      {
        return 0.0;
      }
      const int          power = -std::ilogb(largest);
      const double       scale = std::ldexp(1.0, power);
      const Displacement first = Translated(points[1], points[0], scale);
      const Displacement last  = Translated(points.back(), points[0], scale);
      Twofold::Double    area{};
      if (points.size() == 3)
      {
        area = Twofold::Divide(Wedged(first, last), Twofold::Lifted(3.0));
      }
      else
      {
        const Displacement    second = Translated(points[2], points[0], scale);
        const Twofold::Double twice  = Wedged(second, last);
        const Twofold::Double sum    = Twofold::Add(
            Twofold::Add(Wedged(first, second), Wedged(first, last)), Twofold::Add(twice, twice));
        area = Twofold::Divide(Twofold::Multiply(sum, Twofold::Lifted(3.0)), Twofold::Lifted(20.0));
      }
      return std::ldexp(area.High + area.Low, -2 * power);
    }

    /**
     *  A coordinate pair in double-double, each with its bound.
     */
    using Paired = std::array<Sure::Number, 2>;

    /**
     *  The highest degree whose area Swathed works out: binomials of twice it less one are exact in
     *  doubles.
     */
    constexpr std::size_t Loftiest = 26;

    /**
     *  A polynomial span's Bézier control points less a point, in double-double: each the span's
     *  blossom at its start d - k times and at its end k times, de Boor's mixes worked out from the
     *  poles' exact offsets with every rounding bounded, the rounds at the start one triangle's. A
     *  round's fraction at a cell is the same in every triangle that takes it at the same end, and
     *  is worked out once.
     *
     *  @tparam CountT How many points to take room for: more than the degree.
     *  @param  spline The spline, valid, of degree at most Loftiest.
     *  @param  span   The index of the knot that starts the span, which is not empty.
     *  @param  origin The point.
     *  @return        The control points, the first degree and one.
     */
    template <std::size_t CountT>
    auto Draped(const Spline &spline, std::size_t span, Point origin) -> std::array<Paired, CountT>
    {
      const auto                 degree = static_cast<std::size_t>(spline.Degree);
      const std::vector<double> &knots  = spline.Knots;
      std::array<Paired, CountT> cells  = {};
      std::array<Paired, CountT> out    = {};
      for (std::size_t index = 0; index <= degree; ++index)
      {
        const Point pole = spline.Poles[span - degree + index];
        cells.at(index)  = {Sure::Apart(pole.X, origin.X), Sure::Apart(pole.Y, origin.Y)};
      }
      /**
       *  The fractions of each round at each cell, at the span's start and at its end, and the
       *  rests: the round r's at the cell i at (r - 1) d + i.
       */
      std::array<std::array<Sure::Number, 2>, CountT * CountT> starts = {};
      std::array<std::array<Sure::Number, 2>, CountT * CountT> ends   = {};
      for (std::size_t which = 1; which <= degree; ++which)
      {
        for (std::size_t index = which; index <= degree; ++index)
        {
          const double       from   = knots[span - degree + index];
          const Sure::Number width  = Sure::Apart(knots[span + 1 + index - which], from);
          const Sure::Number ahead  = Sure::Quotient(Sure::Apart(knots[span], from), width);
          const Sure::Number behind = Sure::Quotient(Sure::Apart(knots[span + 1], from), width);
          const std::size_t  cell   = ((which - 1) * degree) + index;
          starts.at(cell)           = {ahead, Sure::Exact(1.0) - ahead};
          ends.at(cell)             = {behind, Sure::Exact(1.0) - behind};
        }
      }
      const auto round =
          [degree](
              std::array<Paired, CountT> &work, std::size_t which,
              const std::array<std::array<Sure::Number, 2>, CountT * CountT> &fractions) -> void
      {
        for (std::size_t index = degree; index >= which; --index)
        {
          const auto &[fraction, rest] = fractions.at(((which - 1) * degree) + index);
          for (std::size_t axis = 0; axis < 2; ++axis)
          {
            work.at(index).at(axis) =
                (work.at(index - 1).at(axis) * rest) + (work.at(index).at(axis) * fraction);
          }
        }
      };
      for (std::size_t taken = 0; taken <= degree; ++taken)
      {
        std::array<Paired, CountT> work = cells;
        for (std::size_t which = taken + 1; which <= degree; ++which)
        {
          round(work, which, ends);
        }
        out.at(degree - taken) = work.at(degree);
        if (taken < degree)
        {
          round(cells, taken + 1, starts);
        }
      }
      return out;
    }

    /**
     *  The cross product of two offsets in double-double.
     *
     *  @param  lhs One offset.
     *  @param  rhs The other.
     *  @return     lhs x rhs, with its bound.
     */
    auto Swept(const Paired &lhs, const Paired &rhs) -> Sure::Number
    {
      return (lhs[0] * rhs[1]) - (lhs[1] * rhs[0]);
    }

    /**
     *  The binomial coefficient C(n, k), exact in doubles for n to 2 Loftiest - 1.
     *
     *  @param  count  n.
     *  @param  chosen k, at most n.
     *  @return        C(n, k).
     */
    auto Chosen(std::size_t count, std::size_t chosen) -> double
    {
      double value = 1.0;
      for (std::size_t index = 0; index < chosen; ++index)
      {
        value = value * static_cast<double>(count - index) / static_cast<double>(index + 1);
      }
      return value;
    }

    /**
     *  The area a polynomial spline encloses with its chord, the double nearest it: by Green's
     *  theorem about its first pole, each span's share in closed form from its Bézier control
     *  points q_i, 1/2 sum over i and j of C(n, i) C(n - 1, j) / (2 C(2 n - 1, i + j)) q_i x (q_j+1
     *  - q_j), the integral of q x q' as products of Bernstein polynomials give it, and the chord's
     *    share 1/2 (end x start); every number in double-double with its bound, the sum rounded
     *    where its bound settles it.
     *
     *  @tparam CountT How many control points to take room for: more than the degree.
     *  @param  spline The spline, valid and polynomial, of degree from 2 to Loftiest.
     *  @return        The area; nothing where its bound does not settle it.
     */
    template <std::size_t CountT>
    auto Swathed(const Spline &spline) -> std::optional<double>
    {
      const auto   degree = static_cast<std::size_t>(spline.Degree);
      const Point  origin = spline.Poles.front();
      Sure::Number twice  = Sure::Exact(0.0);
      Paired       start  = {};
      Paired       end    = {};
      bool         first  = true;
      for (std::size_t span = degree; span < spline.Poles.size(); ++span)
      {
        if (!(spline.Knots[span] < spline.Knots[span + 1]))
        {
          continue;
        }
        const std::array<Paired, CountT> points = Draped<CountT>(spline, span, origin);
        start                                   = first ? points[0] : start;
        end                                     = points.at(degree);
        first                                   = false;
        for (std::size_t total = 0; total < 2 * degree; ++total)
        {
          Sure::Number share = Sure::Exact(0.0);
          for (std::size_t index = total > degree - 1 ? total - (degree - 1) : 0;
               index <= std::min(total, degree); ++index)
          {
            const std::size_t other = total - index;
            const Paired      step  = {points.at(other + 1)[0] - points.at(other)[0],
                                       points.at(other + 1)[1] - points.at(other)[1]};
            share = share
                  + (Sure::Product(Chosen(degree, index), Chosen(degree - 1, other))
                     * Swept(points.at(index), step));
          }
          twice = twice + Sure::Quotient(share, Sure::Exact(2.0 * Chosen((2 * degree) - 1, total)));
        }
      }
      twice = twice + Swept(end, start);
      return Sure::Rounded(
          Sure::Number{.Value = {.High = twice.Value.High / 2.0, .Low = twice.Value.Low / 2.0},
                       .Bound = twice.Bound / 2.0});
    }
  } // namespace

  auto Area(const Line & /*line*/) -> double
  {
    return 0.0;
  }

  auto Area(const Ray & /*ray*/) -> double
  {
    return 0.0;
  }

  auto Area(const Xline & /*xline*/) -> double
  {
    return 0.0;
  }

  auto Area(const Arc &arc) -> double
  {
    return Lens(Elliptic(arc));
  }

  auto Area(const Ellipse &ellipse) -> double
  {
    return Lens(ellipse);
  }

  auto Area(const Polyline &polyline) -> double
  {
    return Polygonal(polyline);
  }

  auto Area(const Spline &spline, Accuracy accuracy) -> double
  {
    if (spline.Degree >= 2
        && std::cmp_less_equal(spline.Degree, Loftiest)
        && Valid(spline)
        && !Rational(spline))
    {
      const std::optional<double> area =
          spline.Degree < 8 ? Swathed<8>(spline) : Swathed<Loftiest + 1>(spline);
      if (area)
      {
        return *area;
      }
    }
    const int way = Valid(spline) ? Sense(spline) : 0;
    if (way == 0)
    {
      return 0.0;
    }
    const double area =
        Encompassed(way > 0 ? spline : Rewound(spline), accuracy, false).Total().Area;
    return way > 0 ? area : -area;
  }

  auto Area(const Bezier::Quadratic &quadratic) -> double
  {
    return Chordal(std::array{quadratic.Start, quadratic.Control, quadratic.End});
  }

  auto Area(const Bezier::Cubic &cubic) -> double
  {
    return Chordal(std::array{cubic.Start, cubic.Out, cubic.In, cubic.End});
  }

  auto Area(const Curve &curve, Accuracy accuracy) -> double
  {
    return std::visit(
        [accuracy](const auto &shape) -> double
        {
          if constexpr (requires { Area(shape); })
          {
            return Area(shape);
          }
          else
          {
            return Area(shape, accuracy);
          }
        },
        curve);
  }

  auto Centroid(const Line & /*line*/) -> std::optional<Point>
  {
    return std::nullopt;
  }

  auto Centroid(const Ray & /*ray*/) -> std::optional<Point>
  {
    return std::nullopt;
  }

  auto Centroid(const Xline & /*xline*/) -> std::optional<Point>
  {
    return std::nullopt;
  }

  auto Centroid(const Arc &arc) -> std::optional<Point>
  {
    return Middle(Elliptic(arc));
  }

  auto Centroid(const Ellipse &ellipse) -> std::optional<Point>
  {
    return Middle(ellipse);
  }

  auto Centroid(const Polyline &polyline) -> std::optional<Point>
  {
    return Centered(Encompassed(polyline, true));
  }

  auto Centroid(const Spline &spline, Accuracy accuracy) -> std::optional<Point>
  {
    const int way = Valid(spline) ? Sense(spline) : 0;
    if (way == 0)
    {
      return std::nullopt;
    }
    return Centered(Encompassed(way > 0 ? spline : Rewound(spline), accuracy, true));
  }

  auto Centroid(const Bezier::Quadratic &quadratic) -> std::optional<Point>
  {
    return Centered(Encompassed(Kernel::Polynomial(quadratic)));
  }

  auto Centroid(const Bezier::Cubic &cubic) -> std::optional<Point>
  {
    return Centered(Encompassed(Kernel::Polynomial(cubic)));
  }

  auto Centroid(const Curve &curve, Accuracy accuracy) -> std::optional<Point>
  {
    return std::visit(
        [accuracy](const auto &shape) -> std::optional<Point>
        {
          if constexpr (requires { Centroid(shape); })
          {
            return Centroid(shape);
          }
          else
          {
            return Centroid(shape, accuracy);
          }
        },
        curve);
  }
} // namespace Geoxy
