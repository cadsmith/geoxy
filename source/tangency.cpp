/***************************************************************************************************
 *  @file      source/tangency.cpp
 *  @brief     Whether a line and a circle, or two circles, touch exactly, by short exact sums; and
 *             where such two near tangent come nearest, in expansions.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "tangency.hpp"

#include "expansion.hpp"
#include "real.hpp"
#include "scratch.hpp"
#include "twofold.hpp"

#include "geoxy/arc.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/vector.hpp"
#include "geoxy/xline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <variant>

namespace Geoxy::Crossings
{
  namespace
  {
    using Twofold::Double;
    using Twofold::FastTwoSum;
    using Twofold::Two;
    using Twofold::TwoProduct;
    using Twofold::TwoSum;

    /**
     *  A number held exactly as a sum of doubles, from the least, none overlapping another, none
     *  zero: Shewchuk's expansion, in the thread's arena.
     */
    using Terms = List<double>;

    /**
     *  A part kept where it is not zero.
     *
     *  @param terms The sum.
     *  @param part  The part.
     */
    void Banked(Terms &terms, double part)
    {
      if (!Zero(part))
      {
        terms.push_back(part);
      }
    }

    /**
     *  A double as a sum.
     *
     *  @param  value The double.
     *  @return       The sum.
     */
    auto Termed(double value) -> Terms
    {
      Terms terms;
      Banked(terms, value);
      return terms;
    }

    /**
     *  A difference of two doubles, exactly.
     *
     *  @param  head The double.
     *  @param  tail The double taken from it.
     *  @return      The difference.
     */
    auto Differenced(double head, double tail) -> Terms
    {
      const Two difference = TwoSum(head, -tail);
      Terms     terms;
      Banked(terms, difference.Low);
      Banked(terms, difference.High);
      return terms;
    }

    /**
     *  A sum with a double added, exactly: Shewchuk's grow-expansion, zeros left out.
     *
     *  @param  terms The sum.
     *  @param  value The double.
     *  @return       The sum.
     */
    auto Grown(std::span<const double> terms, double value) -> Terms
    {
      Terms  grown;
      double carried = value;
      for (const double part : terms)
      {
        const Two sum = TwoSum(carried, part);
        Banked(grown, sum.Low);
        carried = sum.High;
      }
      Banked(grown, carried);
      return grown;
    }

    /**
     *  A sum of two sums, exactly.
     *
     *  @param  first  One.
     *  @param  second The other.
     *  @return        The sum.
     */
    auto Plus(std::span<const double> first, std::span<const double> second) -> Terms
    {
      Terms sum(first.begin(), first.end());
      for (const double part : second)
      {
        sum = Grown(sum, part);
      }
      return sum;
    }

    /**
     *  A sum with its sign changed.
     *
     *  @param  terms The sum.
     *  @return       Its negative.
     */
    auto Negated(std::span<const double> terms) -> Terms
    {
      Terms negated;
      negated.reserve(terms.size());
      for (const double part : terms)
      {
        negated.push_back(-part);
      }
      return negated;
    }

    /**
     *  A sum times a double, exactly: Shewchuk's scale-expansion, zeros left out.
     *
     *  @param  terms The sum.
     *  @param  scale The double.
     *  @return       The product.
     */
    auto Scaled(std::span<const double> terms, double scale) -> Terms
    {
      Terms scaled;
      if (terms.empty())
      {
        return scaled;
      }
      const Two first   = TwoProduct(terms.front(), scale);
      double    carried = first.High;
      Banked(scaled, first.Low);
      for (const double part : terms.subspan(1))
      {
        const Two product = TwoProduct(part, scale);
        const Two sum     = TwoSum(carried, product.Low);
        Banked(scaled, sum.Low);
        const Two next = FastTwoSum(product.High, sum.High);
        Banked(scaled, next.Low);
        carried = next.High;
      }
      Banked(scaled, carried);
      return scaled;
    }

    /**
     *  A product of two sums, exactly.
     *
     *  @param  first  One.
     *  @param  second The other.
     *  @return        The product.
     */
    auto Times(std::span<const double> first, std::span<const double> second) -> Terms
    {
      Terms product;
      for (const double part : second)
      {
        product = Plus(product, Scaled(first, part));
      }
      return product;
    }

    /**
     *  The sign of a sum: its largest part's.
     *
     *  @param  terms The sum.
     *  @return       1, -1 or 0.
     */
    auto Signum(std::span<const double> terms) -> int
    {
      if (terms.empty())
      {
        return 0;
      }
      return terms.back() > 0.0 ? 1 : -1;
    }

    /**
     *  A sum in two doubles: its parts added from the least in double-double, each step within
     *  about 2^-104 of the running sum, which its parts, none overlapping another, keep within a
     *  factor of two of the whole; so within about 2^-99 of it for the few parts these sums have.
     *
     *  @param  terms The sum.
     *  @return       The sum, in two doubles.
     */
    auto Doubled(std::span<const double> terms) -> Double
    {
      Double sum = {};
      for (const double part : terms)
      {
        sum = Twofold::Add(sum, Twofold::Lifted(part));
      }
      return sum;
    }

    /**
     *  The double nearest a quotient of two sums, where their double-double quotient, within 2^-100
     *  or so of the value, leaves no halfway point between doubles within 2^-90 of it, as nearly
     *  always; zero for a numerator of zero; nothing otherwise, or for a value so large or small
     *  the bound does not hold.
     *
     *  @param  top    The numerator.
     *  @param  bottom The denominator, not zero.
     *  @return        The double nearest the quotient.
     */
    auto Quotient(std::span<const double> top, std::span<const double> bottom)
        -> std::optional<double>
    {
      if (top.empty())
      {
        return 0.0;
      }
      const Double value = Twofold::Divide(Doubled(top), Doubled(bottom));
      const double size  = std::abs(value.High);
      const double bound = 0x1p-90 * size;
      const bool   held  = size > 0x1p-900 && size < 0x1p900;
      if (held && Spacing(value.High) - (2.0 * std::abs(value.Low)) > 2.0 * bound)
      {
        return value.High;
      }
      return std::nullopt;
    }

    /**
     *  A straight curve as the short way takes it: a point on it and its way, exactly.
     */
    struct Ruler
    {
      Point                Through = {}; /**< Where its parameter is 0. */
      std::array<Terms, 2> Way     = {}; /**< Its way, exactly.         */
    };

    /**
     *  A line, a ray or an endless line as the short way takes it; nothing for any other curve.
     *
     *  @param  curve The curve.
     *  @return       It, straight.
     */
    auto Straightened(const Curve &curve) -> std::optional<Ruler>
    {
      if (const auto *line = std::get_if<Line>(&curve))
      {
        return Ruler{.Through = line->Start,
                     .Way     = {Differenced(line->End.X, line->Start.X),
                                 Differenced(line->End.Y, line->Start.Y)}};
      }
      if (const auto *ray = std::get_if<Ray>(&curve))
      {
        return Ruler{.Through = ray->Start, .Way = {Termed(ray->Way.X), Termed(ray->Way.Y)}};
      }
      if (const auto *xline = std::get_if<Xline>(&curve))
      {
        return Ruler{.Through = xline->Through,
                     .Way     = {Termed(xline->Way.X), Termed(xline->Way.Y)}};
      }
      return std::nullopt;
    }

    /**
     *  A point from its coordinates, where both settled.
     *
     *  @param  x Its x.
     *  @param  y Its y.
     *  @return   The touch: decided, there; undecided where either did not settle.
     */
    auto Settled(std::optional<double> x, std::optional<double> y) -> Touching
    {
      if (!x || !y)
      {
        return {};
      }
      return {.Decided = true, .At = Point{.X = *x, .Y = *y}};
    }

    /**
     *  Whether a straight curve's carrier touches a circle exactly: (d x g)^2 = r^2 |d|^2, g = P -
     *  O, and where: P - d (g . d) / |d|^2.
     *
     *  @param  straight The straight curve.
     *  @param  arc      The arc whose circle.
     *  @return          The touch.
     */
    auto Lined(const Ruler &straight, const Arc &arc) -> Touching
    {
      const auto &[wayX, wayY] = straight.Way;
      const Terms gapX         = Differenced(straight.Through.X, arc.Center.X);
      const Terms gapY         = Differenced(straight.Through.Y, arc.Center.Y);
      const Terms square       = Plus(Times(wayX, wayX), Times(wayY, wayY));
      if (square.empty())
      {
        return {};
      }
      const Terms cross   = Plus(Times(wayX, gapY), Negated(Times(wayY, gapX)));
      const Two   radius  = TwoProduct(arc.Radius, arc.Radius);
      const Terms squared = Plus(Termed(radius.Low), Termed(radius.High));
      const Terms gap     = Plus(Times(cross, cross), Negated(Times(squared, square)));
      if (Signum(gap) != 0)
      {
        return {.Decided = true, .At = std::nullopt};
      }
      const Terms along = Plus(Times(gapX, wayX), Times(gapY, wayY));
      const auto  axis  = [&square, &along](double start, const Terms &way) -> std::optional<double>
      {
        return Quotient(Plus(Scaled(square, start), Negated(Times(way, along))), square);
      };
      return Settled(axis(straight.Through.X, wayX), axis(straight.Through.Y, wayY));
    }

    /**
     *  Whether two circles touch exactly: |g|^2 = s^2, g = c2 - c1, s the radii's sum or
     *  difference, and where: c1 + g r1 / s.
     *
     *  @param  first  One arc's circle.
     *  @param  second The other's.
     *  @return        The touch.
     */
    auto Circled(const Arc &first, const Arc &second) -> Touching
    {
      const Terms gapX   = Differenced(second.Center.X, first.Center.X);
      const Terms gapY   = Differenced(second.Center.Y, first.Center.Y);
      const Terms square = Plus(Times(gapX, gapX), Times(gapY, gapY));
      if (square.empty())
      {
        return {};
      }
      for (const double sign : {1.0, -1.0})
      {
        const Terms reach = Differenced(first.Radius, -sign * second.Radius);
        if (Signum(Plus(square, Negated(Times(reach, reach)))) == 0)
        {
          const auto axis = [&reach, &first](double       centre,
                                             const Terms &gap) -> std::optional<double>
          {
            return Quotient(Plus(Scaled(reach, centre), Scaled(gap, first.Radius)), reach);
          };
          return Settled(axis(first.Center.X, gapX), axis(first.Center.Y, gapY));
        }
      }
      return {.Decided = true, .At = std::nullopt};
    }

    /**
     *  A quadratic's root over a denominator, exactly: (top + side bottom sqrt(root)) / below, as
     *  the exact crossing holds each coordinate of a root.
     */
    struct Radical
    {
      Terms Top    = {}; /**< What the root's share is added to. */
      Terms Bottom = {}; /**< What the square root is times.     */
      Terms Root   = {}; /**< Under the square root, above zero. */
      Terms Below  = {}; /**< What all is divided by.            */
    };

    /**
     *  The double nearest a root over a denominator, where its double-double estimate settles it:
     *  where the top and the root's share have one sign, their sum; where they cancel, its
     *  conjugate's quotient, (t^2 - s^2) / (t - s), the numerator exact, so that no digit is lost
     *  however small the value, and a numerator of zero a value of zero exactly. Nothing where it
     *  does not settle, or the terms lie outside 2^-300 to 2^300, where these bounds may not hold.
     *
     *  @param  radical The root.
     *  @param  side    Which root: 1 or -1.
     *  @return         The double nearest it.
     */
    auto Surd(const Radical &radical, double side) -> std::optional<double>
    {
      const Double top   = Doubled(radical.Top);
      const Double share = Twofold::Multiply(
          Twofold::Multiply(Doubled(radical.Bottom), Twofold::Root(Doubled(radical.Root))),
          Twofold::Lifted(side));
      const Double below = Doubled(radical.Below);
      const bool   cancels =
          (top.High > 0.0 && share.High < 0.0) || (top.High < 0.0 && share.High > 0.0);
      const double larger = std::max(std::abs(top.High), std::abs(share.High));
      const bool   inside = larger > 0x1p-300 && larger < 0x1p300;
      if (!inside)
      {
        return std::nullopt;
      }
      Double value = Twofold::Divide(Twofold::Add(top, share), below);
      if (cancels)
      {
        const Terms gap = Plus(Times(radical.Top, radical.Top),
                               Negated(Times(Times(radical.Bottom, radical.Bottom), radical.Root)));
        if (gap.empty())
        {
          return 0.0;
        }
        value = Twofold::Divide(Doubled(gap),
                                Twofold::Multiply(Twofold::Difference(top, share), below));
      }
      /**
       *  Terms well inside the doubles' range keep the value there, or past it, where the step is
       *  not a number and nothing settles.
       */
      const double size = std::abs(value.High);
      if (Spacing(value.High) - (2.0 * std::abs(value.Low)) > 2.0 * 0x1p-88 * size)
      {
        return value.High;
      }
      return std::nullopt;
    }

    /**
     *  The crossing a hit is at, of a quadratic's two roots each held as a radical per coordinate:
     *  the nearer the hit, the first where they are as near, as the exact crossing chooses; nothing
     *  where any coordinate does not settle.
     *
     *  @param  axes The two coordinates' radicals.
     *  @param  near Where the hit is.
     *  @return      The crossing.
     */
    auto Chosen(const std::array<Radical, 2> &axes, Point near) -> std::optional<Point>
    {
      const auto root = [&axes](double side) -> std::optional<Point>
      {
        const std::optional<double> x = Surd(axes[0], side);
        const std::optional<double> y = x ? Surd(axes[1], side) : std::nullopt;
        if (!y)
        {
          return std::nullopt;
        }
        return Point{.X = *x, .Y = *y};
      };
      const std::optional<Point> one   = root(1.0);
      const std::optional<Point> other = one ? root(-1.0) : std::nullopt;
      if (!other)
      {
        return std::nullopt;
      }
      return Distance(*one, near) <= Distance(*other, near) ? *one : *other;
    }

    /**
     *  Where a straight curve's carrier crosses a circle near a hit, as the exact crossing finds
     *  it: (P a - d b +- d sqrt(b^2 - a c)) / a, a = d . d, b = d . (P - O), c = |P - O|^2 - r^2.
     *
     *  @param  straight The straight curve.
     *  @param  arc      The arc whose circle.
     *  @param  near     Where the hit is.
     *  @return          The crossing, decided; undecided where the carriers do not cross, or a
     *                   coordinate does not settle.
     */
    auto Transected(const Ruler &straight, const Arc &arc, Point near) -> Touching
    {
      const auto &[wayX, wayY] = straight.Way;
      const Terms gapX         = Differenced(straight.Through.X, arc.Center.X);
      const Terms gapY         = Differenced(straight.Through.Y, arc.Center.Y);
      const Terms square       = Plus(Times(wayX, wayX), Times(wayY, wayY));
      const Terms along        = Plus(Times(gapX, wayX), Times(gapY, wayY));
      const Two   radius       = TwoProduct(arc.Radius, arc.Radius);
      const Terms rest         = Plus(Plus(Times(gapX, gapX), Times(gapY, gapY)),
                                      Negated(Plus(Termed(radius.Low), Termed(radius.High))));
      const Terms root         = Plus(Times(along, along), Negated(Times(square, rest)));
      if (Signum(square) == 0 || Signum(root) <= 0)
      {
        return {};
      }
      const auto axis = [&](double start, const Terms &way) -> Radical
      {
        return {.Top    = Plus(Scaled(square, start), Negated(Times(way, along))),
                .Bottom = way,
                .Root   = root,
                .Below  = square};
      };
      const std::optional<Point> point =
          Chosen({axis(straight.Through.X, wayX), axis(straight.Through.Y, wayY)}, near);
      return {.Decided = point.has_value(), .At = point};
    }

    /**
     *  Where two circles cross near a hit, as the exact crossing finds it: c1 + (k g +- g' sqrt(4
     *  |g|^2 r1^2 - k^2)) / (2 |g|^2), g = c2 - c1, g' it turned a quarter, k = |g|^2 + r1^2
     *  - r2^2.
     *
     *  @param  first  One arc's circle.
     *  @param  second The other's.
     *  @param  near   Where the hit is.
     *  @return        The crossing, decided; undecided where they do not cross, or a coordinate
     *                 does not settle.
     */
    auto Transected(const Arc &first, const Arc &second, Point near) -> Touching
    {
      const Terms gapX   = Differenced(second.Center.X, first.Center.X);
      const Terms gapY   = Differenced(second.Center.Y, first.Center.Y);
      const Terms square = Plus(Times(gapX, gapX), Times(gapY, gapY));
      const Two   inner  = TwoProduct(first.Radius, first.Radius);
      const Two   outer  = TwoProduct(second.Radius, second.Radius);
      const Terms one    = Plus(Termed(inner.Low), Termed(inner.High));
      const Terms chord =
          Plus(Plus(square, one), Negated(Plus(Termed(outer.Low), Termed(outer.High))));
      const Terms root = Plus(Scaled(Times(square, one), 4.0), Negated(Times(chord, chord)));
      if (Signum(square) == 0 || Signum(root) <= 0)
      {
        return {};
      }
      const Terms twice = Scaled(square, 2.0);
      const auto  axis  = [&](double centre, const Terms &gap, const Terms &across) -> Radical
      {
        return {.Top    = Plus(Scaled(twice, centre), Times(chord, gap)),
                .Bottom = across,
                .Root   = root,
                .Below  = twice};
      };
      const std::optional<Point> point = Chosen(
          {axis(first.Center.X, gapX, Negated(gapY)), axis(first.Center.Y, gapY, gapX)}, near);
      return {.Decided = point.has_value(), .At = point};
    }
  } // namespace

  auto Tangency(const std::array<const Curve *, 2> &curves) -> Touching
  {
    const auto *round  = std::get_if<Arc>(curves[0]);
    const auto *across = std::get_if<Arc>(curves[1]);
    if (round != nullptr && across != nullptr)
    {
      return Circled(*round, *across);
    }
    if (across != nullptr)
    {
      return Straightened(*curves[0])
          .transform([across](const Ruler &straight) -> Touching
                     { return Lined(straight, *across); })
          .value_or(Touching{});
    }
    if (round != nullptr)
    {
      return Straightened(*curves[1])
          .transform([round](const Ruler &straight) -> Touching { return Lined(straight, *round); })
          .value_or(Touching{});
    }
    return {};
  }

  auto Shortly(const std::array<const Curve *, 2> &curves, Point near) -> Touching
  {
    const auto *round  = std::get_if<Arc>(curves[0]);
    const auto *across = std::get_if<Arc>(curves[1]);
    if (round != nullptr && across != nullptr)
    {
      return Transected(*round, *across, near);
    }
    if (across != nullptr)
    {
      return Straightened(*curves[0])
          .transform([across, near](const Ruler &straight) -> Touching
                     { return Transected(straight, *across, near); })
          .value_or(Touching{});
    }
    if (round != nullptr)
    {
      return Straightened(*curves[1])
          .transform([round, near](const Ruler &straight) -> Touching
                     { return Transected(straight, *round, near); })
          .value_or(Touching{});
    }
    return {};
  }

  auto Nearing(const std::array<const Arc *, 2> &arcs, const std::array<double, 2> &ways) -> Point
  {
    using Expansion::Number;
    const Arc   &first  = *arcs[0];
    const Arc   &second = *arcs[1];
    const Number gapX =
        Expansion::Difference(Expansion::Of(second.Center.X), Expansion::Of(first.Center.X));
    const Number gapY =
        Expansion::Difference(Expansion::Of(second.Center.Y), Expansion::Of(first.Center.Y));
    const Number square =
        Expansion::Sum(Expansion::Product(gapX, gapX), Expansion::Product(gapY, gapY));
    const Number reach = Expansion::Sum(Expansion::Of(ways[0] * first.Radius),
                                        Expansion::Of(ways[1] * second.Radius));
    const Number twice = Expansion::Scaled(square, 2.0);
    const auto   axis  = [&square, &reach, &twice](std::array<double, 2> centers,
                                                   const Number         &gap) -> double
    {
      const Number middle = Expansion::Sum(Expansion::Of(centers[0]), Expansion::Of(centers[1]));
      return Expansion::Rounded({.Top    = Expansion::Product(middle, square),
                                 .Bottom = Expansion::Product(reach, gap),
                                 .Root   = square,
                                 .Below  = twice},
                                1);
    };
    return {.X = axis({first.Center.X, second.Center.X}, gapX),
            .Y = axis({first.Center.Y, second.Center.Y}, gapY)};
  }

  namespace
  {
    /**
     *  A straight curve's carrier exactly: a point it runs through, and its way, a line's the
     *  difference of its ends.
     */
    struct Tramline
    {
      Point                            Through = {}; /**< A point it runs through. */
      std::array<Expansion::Number, 2> Way     = {}; /**< Its way, exactly.        */
    };

    /**
     *  A line, a ray or an endless line's carrier exactly.
     *
     *  @pre          The curve is a line, a ray or an endless line.
     *  @param  curve The curve.
     *  @return       Its carrier.
     */
    auto Trammed(const Curve &curve) -> Tramline
    {
      const auto exactly = [](Vector way) -> std::array<Expansion::Number, 2>
      {
        return {Expansion::Of(way.X), Expansion::Of(way.Y)};
      };
      if (const auto *line = std::get_if<Line>(&curve))
      {
        return {
            .Through = line->Start,
            .Way     = {
                Expansion::Difference(Expansion::Of(line->End.X), Expansion::Of(line->Start.X)),
                Expansion::Difference(Expansion::Of(line->End.Y), Expansion::Of(line->Start.Y))}};
      }
      if (const auto *ray = std::get_if<Ray>(&curve))
      {
        return {.Through = ray->Start, .Way = exactly(ray->Way)};
      }
      const auto &xline = std::get<Xline>(curve);
      return {.Through = xline.Through, .Way = exactly(xline.Way)};
    }

  } // namespace

  auto Nearing(const Curve &straight, const Arc &arc) -> Point
  {
    using Expansion::Number;
    const Tramline               carrier = Trammed(straight);
    const std::array<Number, 2> &way     = carrier.Way;
    const std::array<Number, 2>  gap     = {
        Expansion::Difference(Expansion::Of(carrier.Through.X), Expansion::Of(arc.Center.X)),
        Expansion::Difference(Expansion::Of(carrier.Through.Y), Expansion::Of(arc.Center.Y))};
    const Number square =
        Expansion::Sum(Expansion::Product(way[0], way[0]), Expansion::Product(way[1], way[1]));
    const Number cross = Expansion::Difference(Expansion::Product(way[0], gap[1]),
                                               Expansion::Product(way[1], gap[0]));
    const Number off = Expansion::Sign(cross) < 0 ? Expansion::Difference(Number{}, cross) : cross;
    const Number along =
        Expansion::Sum(Expansion::Product(gap[0], way[0]), Expansion::Product(gap[1], way[1]));
    const Number below = Expansion::Scaled(Expansion::Product(square, off), 2.0);
    const auto   axis  = [&](std::size_t index, double centre) -> double
    {
      const Number aside = Expansion::Difference(Expansion::Product(gap.at(index), square),
                                                 Expansion::Product(way.at(index), along));
      const Number top   = Expansion::Sum(Expansion::Scaled(square, 2.0 * centre), aside);
      return Expansion::Rounded({.Top    = Expansion::Product(top, off),
                                 .Bottom = Expansion::Scaled(aside, arc.Radius),
                                 .Root   = square,
                                 .Below  = below},
                                1);
    };
    return {.X = axis(0, arc.Center.X), .Y = axis(1, arc.Center.Y)};
  }
} // namespace Geoxy::Crossings
