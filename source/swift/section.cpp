/***************************************************************************************************
 *  @file      source/swift/section.cpp
 *  @brief     Quadratic and cubic pieces of Bézier curves and polynomial splines as the swift
 *             kernels take them: a span's control points by its blossom in doubles, within a slack
 *             of the true ones, and its true points by de Boor's algorithm in double-double.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "swift/section.hpp"

#include "basis.hpp"
#include "scratch.hpp"
#include "sure.hpp"
#include "swift/casteljau.hpp"
#include "swift/lane.hpp"
#include "swift/roots.hpp"
#include "trig.hpp"

#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/box.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <type_traits>
#include <utility>

namespace Geoxy::Swift
{
  namespace
  {
    /**
     *  The least a fraction of knot differences may be and not be nothing: a fraction below the
     *  normal doubles is rounded to a step of them, not to a share of itself, which the slack does
     *  not charge.
     */
    constexpr double Faintest = 0x1p-1000;

    /**
     *  A section of a Bézier curve's control points, exact.
     *
     *  @param  poles The control points, three or four.
     *  @return       The section.
     */
    auto Exacted(std::span<const Point> poles) -> Section
    {
      Section section = {.Degree = poles.size() - 1};
      for (std::size_t index = 0; index < poles.size(); ++index)
      {
        section.Poles.at(index) = poles[index];
        section.Hull            = Grown(section.Hull, poles[index]);
      }
      return section;
    }

    /**
     *  Where a span's blossom is taken: the span, how many of its rounds take its end, and the
     *  point its poles are taken from.
     */
    struct Blossom
    {
      std::size_t Span = 0;  /**< The index of the knot that starts the span.         */
      std::size_t Ends = 0;  /**< How many rounds take the span's end, the last ones. */
      Point       Base = {}; /**< The point the poles are taken from.                 */
    };

    /**
     *  A span's Bézier control point: the blossom at the span's start taken so many times, then at
     *  its end, by de Boor's rounds in doubles.
     *
     *  @tparam CountT  Room for the span's poles: its degree and one at the least.
     *  @param  spline  The spline, valid and polynomial.
     *  @param  blossom The span, its ends' rounds, and the base its poles are taken from, each gap
     *                  rounded once.
     *  @return         The control point less the base; nothing where a fraction falls below the
     *                  normal doubles.
     */
    auto Blossomed(const Spline &spline, const Blossom &blossom) -> std::optional<Point>
    {
      const auto                      degree = static_cast<std::size_t>(spline.Degree);
      const std::span<const double>   knots  = spline.Knots;
      const std::size_t               span   = blossom.Span;
      const std::size_t               ends   = blossom.Ends;
      std::array<Point, Broadest + 1> cells  = {};
      for (std::size_t index = 0; index <= degree; ++index)
      {
        const Point pole = spline.Poles[span - degree + index];
        cells.at(index)  = {.X = pole.X - blossom.Base.X, .Y = pole.Y - blossom.Base.Y};
      }
      for (std::size_t round = 1; round <= degree; ++round)
      {
        const double value = round + ends > degree ? knots[span + 1] : knots[span];
        for (std::size_t index = degree; index >= round; --index)
        {
          const double start    = knots[span - degree + index];
          const double end      = knots[span + 1 + index - round];
          const double ahead    = value - start;
          const double fraction = ahead / (end - start);
          if (ahead > 0.0 && !(fraction >= Faintest))
          {
            return std::nullopt;
          }
          const double rest = 1.0 - fraction;
          const Point  lhs  = cells.at(index - 1);
          const Point  rhs  = cells.at(index);
          cells.at(index)   = {.X = (lhs.X * rest) + (rhs.X * fraction),
                               .Y = (lhs.Y * rest) + (rhs.Y * fraction)};
        }
      }
      return cells.at(degree);
    }

    /**
     *  The point a span's poles are taken from: a broad section's first pole, so that a spline far
     *  from nothing keeps the digits of its size in its control points; nothing for a narrow one,
     *  whose control points are the true ones' rounded.
     *
     *  @tparam SectionT The kind of section.
     *  @param  spline   The spline.
     *  @param  span     The index of the knot that starts the span.
     *  @return          The base.
     */
    template <typename SectionT>
    auto Based(const Spline &spline, std::size_t span) -> Point
    {
      if constexpr (std::is_same_v<SectionT, Broad>)
      {
        return spline.Poles[span - static_cast<std::size_t>(spline.Degree)];
      }
      else
      {
        return {};
      }
    }

    /**
     *  A span of a polynomial spline as a section.
     *
     *  @tparam SectionT The kind of section, with room for the span's control points.
     *  @param  spline   The spline, valid and polynomial, of a degree from two the section has room
     *                   for, its numbers within the range the kernels take.
     *  @param  span     The index of the knot that starts the span, which is not empty.
     *  @return          The section; nothing where a fraction falls below the normal doubles.
     */
    template <typename SectionT>
    auto Spanned(const Spline &spline, std::size_t span) -> std::optional<SectionT>
    {
      const auto  degree  = static_cast<std::size_t>(spline.Degree);
      SectionT    section = {.Degree = degree,
                             .Span   = {.Start = spline.Knots[span], .End = spline.Knots[span + 1]},
                             .Source = &spline,
                             .Knot   = span,
                             .Base   = Based<SectionT>(spline, span)};
      const Point base    = section.Base;
      double      most    = 0.0;
      for (std::size_t index = span - degree; index <= span; ++index)
      {
        const Point pole = spline.Poles[index];
        most             = std::max({most, std::abs(pole.X - base.X), std::abs(pole.Y - base.Y)});
      }
      Box hull = {.Min = {.X = 0.0, .Y = 0.0}, .Max = {.X = 0.0, .Y = 0.0}};
      for (std::size_t index = 0; index <= degree; ++index)
      {
        const std::optional<Point> pole =
            Blossomed(spline, {.Span = span, .Ends = index, .Base = base});
        if (!pole)
        {
          return std::nullopt;
        }
        section.Poles.at(index) = *pole;
        hull                    = index == 0 ? Box{.Min = *pole, .Max = *pole} : Grown(hull, *pole);
      }
      /**
       *  Each cell the pole less the base, rounded once, within u of the largest coordinate M of
       *  the gaps; each round within ten steps u of the largest coordinate of the cells it mixes,
       *  which are within M grown by the slack, and what it carries as it was: under 11 n u M after
       *  n rounds. The hull is padded by twice the slack, so that its edges' own rounding keeps the
       *  true piece within, and moved back by the base, each edge's rounding padded.
       */
      section.Slack =
          (11.0 * static_cast<double>(degree) * Unit * most * Sure::Grow * Sure::Grow) + Sure::Tiny;
      const double   pad   = 2.0 * section.Slack;
      const double   shift = Unit
                           * (std::abs(base.X)
                              + std::abs(base.Y)
                              + std::abs(hull.Min.X)
                              + std::abs(hull.Min.Y)
                              + std::abs(hull.Max.X)
                              + std::abs(hull.Max.Y))
                           * 2.0;
      constexpr auto Wide  = static_cast<double>(std::is_same_v<SectionT, Broad>);
      const double   room  = pad + (Wide * shift);
      section.Hull         = {
          .Min = {.X = (hull.Min.X + base.X) - room, .Y = (hull.Min.Y + base.Y) - room},
          .Max = {.X = (hull.Max.X + base.X) + room, .Y = (hull.Max.Y + base.Y) + room}};
      return section;
    }

    /**
     *  A rational span's homogeneous Bézier control point: Blossomed on the cells (w x, w y, w),
     *  the products rounded once.
     *
     *  @tparam CountT  Room for the span's poles: its degree and one at the least.
     *  @param  spline  The spline, valid, of degree two at the least.
     *  @param  blossom The span, its ends' rounds, and the base its poles are taken from: each cell
     *                  w (P - B), the gap and the product rounded once each.
     *  @return         w (P - B)'s x and y and w; nothing where a fraction falls below the normal
     *                  doubles.
     */
    auto Lifted(const Spline &spline, const Blossom &blossom)
        -> std::optional<std::array<double, 3>>
    {
      const std::size_t             span   = blossom.Span;
      const std::size_t             ends   = blossom.Ends;
      const auto                    degree = static_cast<std::size_t>(spline.Degree);
      const std::span<const double> knots  = spline.Knots;
      std::array<std::array<double, 3>, Broadest + 1> cells = {};
      for (std::size_t index = 0; index <= degree; ++index)
      {
        const Point  pole   = spline.Poles[span - degree + index];
        const double weight = spline.Weights[span - degree + index];
        cells.at(index) = {weight * (pole.X - blossom.Base.X), weight * (pole.Y - blossom.Base.Y),
                           weight};
      }
      for (std::size_t round = 1; round <= degree; ++round)
      {
        const double value = round + ends > degree ? knots[span + 1] : knots[span];
        for (std::size_t index = degree; index >= round; --index)
        {
          const double start    = knots[span - degree + index];
          const double end      = knots[span + 1 + index - round];
          const double ahead    = value - start;
          const double fraction = ahead / (end - start);
          if (ahead > 0.0 && !(fraction >= Faintest))
          {
            return std::nullopt;
          }
          const double rest = 1.0 - fraction;
          for (std::size_t axis = 0; axis < 3; ++axis)
          {
            cells.at(index).at(axis) =
                (cells.at(index - 1).at(axis) * rest) + (cells.at(index).at(axis) * fraction);
          }
        }
      }
      return cells.at(degree);
    }

    /**
     *  A span of a rational spline as a section: its homogeneous control points and weights within
     *  the slack, the hull about the true control points, each within (s + |P| s) / (w - s) of w P
     *  / w worked out, the division's rounding aside.
     *
     *  @tparam SectionT The kind of section, with room for the span's control points.
     *  @param  spline   The spline, valid, of a degree from two the section has room for, its
     *                   numbers within the range the kernels take.
     *  @param  span     The index of the knot that starts the span, which is not empty.
     *  @return          The section; nothing where a fraction falls below the normal doubles, or a
     *                   weight may be nothing.
     */
    auto Rated(const Spline &spline, std::size_t span) -> std::optional<Broad>
    {
      const auto  degree  = static_cast<std::size_t>(spline.Degree);
      Broad       section = {.Degree = degree,
                             .Span   = {.Start = spline.Knots[span], .End = spline.Knots[span + 1]},
                             .Source = &spline,
                             .Knot   = span,
                             .Rated  = true,
                             .Base   = Based<Broad>(spline, span)};
      const Point base    = section.Base;
      double      most    = 0.0;
      for (std::size_t index = span - degree; index <= span; ++index)
      {
        const Point  pole   = spline.Poles[index];
        const double weight = spline.Weights[index];
        most                = std::max({most, std::abs(weight * (pole.X - base.X)),
                                        std::abs(weight * (pole.Y - base.Y)), weight});
      }
      section.Slack =
          (11.0 * static_cast<double>(degree) * Unit * most * Sure::Grow * Sure::Grow) + Sure::Tiny;
      for (std::size_t index = 0; index <= degree; ++index)
      {
        const std::optional<std::array<double, 3>> cell =
            Lifted(spline, {.Span = span, .Ends = index, .Base = base});
        if (!cell || !((*cell)[2] - section.Slack > 0.0))
        {
          return std::nullopt;
        }
        const auto [across, upward, weight] = *cell;
        const Point  near                   = {.X = across / weight, .Y = upward / weight};
        const Point  pole                   = {.X = near.X + base.X, .Y = near.Y + base.Y};
        const double size                   = std::max(std::abs(near.X), std::abs(near.Y));
        /**
         *  The pole w P / w taken back from its base, rounded once more: within a step u of its own
         *  size past what the base leaves.
         */
        const double moved        = std::max(std::abs(pole.X), std::abs(pole.Y));
        const double off          = ((section.Slack * (1.0 + size) / (weight - section.Slack))
                                     + (Unit * size)
                                     + (Unit * moved))
                                  * Sure::Grow;
        section.Poles.at(index)   = {.X = across, .Y = upward};
        section.Weights.at(index) = weight;
        section.Hull              = Grown(Grown(section.Hull, pole - Vector{.X = off, .Y = off}),
                                          pole + Vector{.X = off, .Y = off});
      }
      return section;
    }

    /**
     *  Whether a spline is cut into sections with room for so many control points: valid, of a
     *  degree from two the room holds, its numbers within the range the kernels take, and a
     *  rational one's weights above 2^-200.
     *
     *  @param  spline The spline.
     *  @param  room   The room.
     *  @return        Whether it is.
     */
    auto Cuttable(const Spline &spline, std::size_t room) -> bool
    {
      if (spline.Degree < 2 || !std::cmp_less(spline.Degree, room) || !Valid(spline))
      {
        return false;
      }
      return std::ranges::all_of(spline.Poles,
                                 [](Point pole) -> bool { return Held({pole.X, pole.Y}); })
          && std::ranges::all_of(spline.Knots, [](double knot) -> bool { return Held({knot}); })
          && (!Rational(spline)
              || std::ranges::all_of(spline.Weights, [](double weight) -> bool
                                     { return weight > Narrowest && Held({weight}); }));
    }

    /**
     *  Whether a span is passed over: empty, or its poles' box, which it lies in, missing a box.
     *
     *  @param  spline The spline.
     *  @param  span   The index of the knot that starts the span.
     *  @param  near   The box.
     *  @return        Whether it is.
     */
    auto Unspanned(const Spline &spline, std::size_t span, const Box &near) -> bool
    {
      const auto degree = static_cast<std::size_t>(spline.Degree);
      Box        poles  = {.Min = spline.Poles[span], .Max = spline.Poles[span]};
      for (std::size_t index = span - degree; index < span; ++index)
      {
        poles = Grown(poles, spline.Poles[index]);
      }
      return !(spline.Knots[span] < spline.Knots[span + 1])
          || poles.Min.X > near.Max.X
          || poles.Max.X < near.Min.X
          || poles.Min.Y > near.Max.Y
          || poles.Max.Y < near.Min.Y;
    }

    /**
     *  A spline cut into sections, one per span that is not empty, in order: a polynomial one's
     *  spans, or a rational one's.
     *
     *  @tparam SectionT The kind of section: its room bounds the degree.
     *  @param  spline   The spline.
     *  @param  near     The box outside which no span is cut: a span lies in the hull of its poles,
     *                   its weights above zero, and one whose poles' box misses it is left out.
     *  @return          Its sections; nothing where it is not valid, its degree is below two or
     *                   past the room, a number lies outside the range the kernels take, a rational
     *                   one's weight is not above 2^-200, or a span cannot be cut.
     */
    template <typename SectionT>
    auto Cut(const Spline &spline, const Box &near) -> std::optional<List<SectionT>>
    {
      constexpr std::size_t Room = std::tuple_size_v<decltype(SectionT::Poles)>;
      if (!Cuttable(spline, Room))
      {
        return std::nullopt;
      }
      const auto     degree = static_cast<std::size_t>(spline.Degree);
      List<SectionT> sections;
      for (std::size_t span = degree; span < spline.Poles.size(); ++span)
      {
        if (Unspanned(spline, span, near))
        {
          continue;
        }
        std::optional<SectionT> section = std::nullopt;
        if constexpr (std::is_same_v<SectionT, Broad>)
        {
          section = Rational(spline) ? Rated(spline, span) : Spanned<SectionT>(spline, span);
        }
        else
        {
          section = Spanned<SectionT>(spline, span);
        }
        if (!section)
        {
          return std::nullopt;
        }
        sections.push_back(*section);
      }
      return sections;
    }

    /**
     *  A Bézier curve's section's true point at a place on the grid: de Casteljau's algorithm on
     *  its exact control points in double-double, every rounding bounded.
     *
     *  @param  section The section, its control points exact.
     *  @param  along   The place, from 0 to 1, on the grid.
     *  @return         The point's x and y, each with its bound.
     */
    auto Casted(const Section &section, double along) -> std::array<Sure::Number, 2>
    {
      const std::size_t degree    = section.Degree;
      Bounded           abscissas = {};
      Bounded           ordinates = {};
      for (std::size_t index = 0; index <= degree; ++index)
      {
        abscissas.at(index) = Sure::Exact(section.Poles.at(index).X);
        ordinates.at(index) = Sure::Exact(section.Poles.at(index).Y);
      }
      return {Valued(abscissas, degree, along), Valued(ordinates, degree, along)};
    }

    /**
     *  A spline's true point at a place of a span on the grid: de Boor's algorithm on its own poles
     *  in double-double, every rounding bounded. The place is u = a + w t on the span from a to b,
     *  w = b - a; u less a knot k is w t + (a - k), each part exact or bounded, so that a place a
     *  sliver from the span's start keeps its digits. Every fraction is then (u - start) / (end -
     *  start), its rest 1 less it.
     *
     *  @tparam CountT  Room for the span's poles: its degree and one at the least.
     *  @param  section The span's section.
     *  @param  along   The place, from 0 to 1, on the grid.
     *  @return         The point's x and y, each with its bound.
     */
    template <std::size_t CountT, typename SectionT>
    auto Boored(const SectionT &section, double along) -> std::array<Sure::Number, 2>
    {
      const std::size_t             degree = section.Degree;
      const Spline                 &spline = *section.Source;
      const std::span<const double> knots  = spline.Knots;
      const std::size_t             span   = section.Knot;
      const double                  start  = knots[span];
      const Sure::Number ahead = Sure::Apart(knots[span + 1], start) * Sure::Exact(along);
      std::array<Sure::Number, CountT> abscissas = {};
      std::array<Sure::Number, CountT> ordinates = {};
      for (std::size_t index = 0; index <= degree; ++index)
      {
        const Point pole    = spline.Poles[span - degree + index];
        abscissas.at(index) = Sure::Exact(pole.X);
        ordinates.at(index) = Sure::Exact(pole.Y);
      }
      for (std::size_t round = 1; round <= degree; ++round)
      {
        for (std::size_t index = degree; index >= round; --index)
        {
          const double       from  = knots[span - degree + index];
          const double       until = knots[span + 1 + index - round];
          const Sure::Number fraction =
              Sure::Quotient(ahead + Sure::Apart(start, from), Sure::Apart(until, from));
          const Sure::Number rest = Sure::Exact(1.0) - fraction;
          abscissas.at(index) = (abscissas.at(index - 1) * rest) + (abscissas.at(index) * fraction);
          ordinates.at(index) = (ordinates.at(index - 1) * rest) + (ordinates.at(index) * fraction);
        }
      }
      return {abscissas.at(degree), ordinates.at(degree)};
    }

    /**
     *  The curve's parameter at a place of a section, as Spread gives it.
     *
     *  @tparam SectionT The kind of section.
     *  @param  section  The section.
     *  @param  along    The place, with its bound.
     *  @return          The curve's parameter, with its bound.
     */
    template <typename SectionT>
    auto Spreading(const SectionT &section, const Sure::Number &along) -> Sure::Number
    {
      if (section.Source == nullptr)
      {
        return along;
      }
      return Sure::Exact(section.Span.Start)
           + (Sure::Apart(section.Span.End, section.Span.Start) * along);
    }
  } // namespace

  auto Sectioned(const Bezier::Quadratic &quadratic) -> Section
  {
    return Exacted(std::array{quadratic.Start, quadratic.Control, quadratic.End});
  }

  auto Sectioned(const Bezier::Cubic &cubic) -> Section
  {
    return Exacted(std::array{cubic.Start, cubic.Out, cubic.In, cubic.End});
  }

  auto Sectioned(const Spline &spline) -> std::optional<List<Section>>
  {
    constexpr double Endless = std::numeric_limits<double>::infinity();
    if (Rational(spline))
    {
      return std::nullopt;
    }
    return Cut<Section>(
        spline, {.Min = {.X = -Endless, .Y = -Endless}, .Max = {.X = Endless, .Y = Endless}});
  }

  auto Weighed(const Spline &spline, const Box &near) -> std::optional<List<Broad>>
  {
    return Cut<Broad>(spline, near);
  }

  auto Homogeneously(const Broad &section, const Sure::Number &along) -> std::array<Sure::Number, 3>
  {
    const std::size_t             degree = section.Degree;
    const Spline                 &spline = *section.Source;
    const std::span<const double> knots  = spline.Knots;
    const std::size_t             span   = section.Knot;
    const double                  start  = knots[span];
    const Sure::Number            ahead  = Sure::Apart(knots[span + 1], start) * along;
    std::array<std::array<Sure::Number, 3>, Broadest + 1> cells = {};
    for (std::size_t index = 0; index <= degree; ++index)
    {
      const Point  pole   = spline.Poles[span - degree + index];
      const double weight = spline.Weights[span - degree + index];
      cells.at(index)     = {Sure::Product(weight, pole.X), Sure::Product(weight, pole.Y),
                             Sure::Exact(weight)};
    }
    for (std::size_t round = 1; round <= degree; ++round)
    {
      for (std::size_t index = degree; index >= round; --index)
      {
        const double       from  = knots[span - degree + index];
        const double       until = knots[span + 1 + index - round];
        const Sure::Number fraction =
            Sure::Quotient(ahead + Sure::Apart(start, from), Sure::Apart(until, from));
        const Sure::Number rest = Sure::Exact(1.0) - fraction;
        for (std::size_t axis = 0; axis < 3; ++axis)
        {
          cells.at(index).at(axis) =
              (cells.at(index - 1).at(axis) * rest) + (cells.at(index).at(axis) * fraction);
        }
      }
    }
    return cells.at(degree);
  }

  auto Fastest(const Broad &section) -> double
  {
    const auto   order    = static_cast<double>(section.Degree);
    const auto   weights  = std::span{section.Weights}.first(section.Degree + 1);
    const double heaviest = std::ranges::max(weights) + section.Slack;
    const double lightest = std::ranges::min(weights) - section.Slack;
    const Box   &hull     = section.Hull;
    const double across   = Trig::Hypot(hull.Max.X - hull.Min.X, hull.Max.Y - hull.Min.Y);
    return 2.0 * order * heaviest * across / lightest * (1.0 + 0x1p-48);
  }

  auto Exactly(const Section &section, double along) -> std::array<Sure::Number, 2>
  {
    return section.Source == nullptr ? Casted(section, along) : Boored<Highest + 1>(section, along);
  }

  auto Exactly(const Broad &section, double along) -> std::array<Sure::Number, 2>
  {
    if (section.Degree > Highest)
    {
      return Boored<Broadest + 1>(section, along);
    }
    /**
     *  A quadratic or cubic piece is worked out in a narrow section's cells, the same numbers in
     *  fewer: its control points and its span are all it reads.
     */
    Section narrow = {.Degree = section.Degree,
                      .Slack  = section.Slack,
                      .Span   = section.Span,
                      .Hull   = section.Hull,
                      .Source = section.Source,
                      .Knot   = section.Knot,
                      .Rated  = section.Rated,
                      .Base   = section.Base};
    std::ranges::copy(std::span{section.Poles}.first(Highest + 1), narrow.Poles.begin());
    return Exactly(narrow, along);
  }

  auto Spread(const Section &section, const Sure::Number &along) -> Sure::Number
  {
    return Spreading(section, along);
  }

  auto Spread(const Broad &section, const Sure::Number &along) -> Sure::Number
  {
    return Spreading(section, along);
  }
} // namespace Geoxy::Swift
