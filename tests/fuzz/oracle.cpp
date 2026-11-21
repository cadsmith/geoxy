/***************************************************************************************************
 *  @file      tests/fuzz/oracle.cpp
 *  @brief     The exact predicates on an input of the fuzzer, against Shewchuk's own.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "oracle.hpp"

#include "../geoxy/exact/original.hpp"
#include "../geoxy/exact/same.hpp"
#include "draw.hpp"
#include "exact/batch.hpp"
#include "exact/circle.hpp"
#include "exact/grid.hpp"
#include "exact/orient.hpp"
#include "exact/triangle.hpp"
#include "exact/two.hpp"
#include "exact/whole.hpp"

#include "geoxy/point.hpp"

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace Geoxy::Fuzz
{
  namespace
  {
    /**
     *  The magnitudes the orientation test and the exact sums of two doubles are exact within, as a
     *  power of two.
     */
    constexpr int OrientRange = 400;

    /**
     *  The magnitudes the circle test is exact within, as a power of two.
     */
    constexpr int CircleRange = 150;

    /**
     *  The most points a batch draws: past three blocks of 256, and an end in the middle of one.
     */
    constexpr unsigned BatchMost = 1000;

    /**
     *  A value the oracle computed, as components smallest first.
     */
    using Components = std::vector<double>;

    /**
     *  What a predicate gave, and what it should have.
     */
    struct Law
    {
      std::string_view What;     /**< How the points were given.   */
      int              Got;      /**< The sign the predicate gave. */
      int              Expected; /**< The sign it should give.     */
    };

    /**
     *  Whether one of the arithmetic's routines gave the original's value.
     */
    struct Check
    {
      std::string_view What; /**< The routine.        */
      bool             Held; /**< Whether it gave it. */
    };

    /**
     *  Sets the original's constants up once.
     */
    void Start()
    {
      static const bool started = [] -> bool
      {
        exactinit();
        return true;
      }();
      static_cast<void>(started);
    }

    /**
     *  Whether a double is zero or of a magnitude within a range.
     *
     *  @param  value    The double.
     *  @param  exponent The range, from 2^-exponent to 2^exponent.
     *  @return          Whether it is in the range; false for NaN and infinity.
     */
    auto Within(double value, int exponent) -> bool
    {
      if (std::fpclassify(value) == FP_ZERO)
      {
        return true;
      }
      const double magnitude = std::abs(value);
      return std::ldexp(1.0, -exponent) <= magnitude && magnitude <= std::ldexp(1.0, exponent);
    }

    /**
     *  Whether both coordinates of every point are within a range.
     *
     *  @param  points   The points.
     *  @param  exponent The range, as Within's.
     *  @return          Whether they are.
     */
    auto Within(std::span<const Point> points, int exponent) -> bool
    {
      return std::ranges::all_of(
          points, [exponent](Point point) -> bool
          { return Within(point.X, exponent) && Within(point.Y, exponent); });
    }

    /**
     *  Whether every component is within a range.
     *
     *  @param  components The components.
     *  @param  exponent   The range, as Within's.
     *  @return            Whether they are.
     */
    auto Within(std::span<const double> components, int exponent) -> bool
    {
      return std::ranges::all_of(components, [exponent](double component) -> bool
                                 { return Within(component, exponent); });
    }

    /**
     *  Points as text that reads back to the same doubles.
     *
     *  @param  points The points.
     *  @return        Each as (x, y) in hexadecimal.
     */
    auto Text(std::span<const Point> points) -> std::string
    {
      std::string text;
      for (const Point point : points)
      {
        text += std::format(" ({:a}, {:a})", point.X, point.Y);
      }
      return text;
    }

    /**
     *  The first law that does not hold, told.
     *
     *  @param  predicate The predicate's name.
     *  @param  laws      The laws.
     *  @param  points    The points they were checked on.
     *  @return           What broke, empty when nothing did.
     */
    auto Told(std::string_view predicate, std::span<const Law> laws, std::span<const Point> points)
        -> std::string
    {
      for (const Law &law : laws)
      {
        if (law.Got != law.Expected)
        {
          return std::format("{}: {} gave {}, the original's sign {}, for{}", predicate, law.What,
                             law.Got, law.Expected, Text(points));
        }
      }
      return {};
    }

    /**
     *  A point mirrored in the y axis, which swaps the sign of both determinants.
     *
     *  @param  point The point.
     *  @return       The point mirrored, exactly.
     */
    auto Mirrored(Point point) -> Point
    {
      return {.X = -point.X, .Y = point.Y};
    }

    /**
     *  A point turned a half turn about the origin, which keeps the sign of both determinants.
     *
     *  @param  point The point.
     *  @return       The point turned, exactly.
     */
    auto Negated(Point point) -> Point
    {
      return {.X = -point.X, .Y = -point.Y};
    }

    /**
     *  The original's orientation sign.
     *
     *  @param  corners The triangle.
     *  @return         The sign of orient2d.
     */
    auto Original(const Exact::Triangle &corners) -> int
    {
      std::array<double, 2> a = {corners.First.X, corners.First.Y};
      std::array<double, 2> b = {corners.Second.X, corners.Second.Y};
      std::array<double, 2> c = {corners.Third.X, corners.Third.Y};
      return Exact::Sign(orient2d(a.data(), b.data(), c.data()));
    }

    /**
     *  The original's circle sign.
     *
     *  @param  quad The circle and the point asked about.
     *  @return      The sign of incircle.
     */
    auto Original(const Quad &quad) -> int
    {
      std::array<double, 2> a = {quad.Circle.First.X, quad.Circle.First.Y};
      std::array<double, 2> b = {quad.Circle.Second.X, quad.Circle.Second.Y};
      std::array<double, 2> c = {quad.Circle.Third.X, quad.Circle.Third.Y};
      std::array<double, 2> d = {quad.Asked.X, quad.Asked.Y};
      return Exact::Sign(incircle(a.data(), b.data(), c.data(), d.data()));
    }

    /**
     *  Whether every coordinate is finite.
     *
     *  @param  points The points.
     *  @return        Whether none is NaN or infinite.
     */
    auto Finite(std::span<const Point> points) -> bool
    {
      return std::ranges::all_of(points, [](Point point) -> bool
                                 { return std::isfinite(point.X) && std::isfinite(point.Y); });
    }

    /**
     *  Finite points scaled by one power of two, the largest coordinate to 2^(exponent / 2), where
     *  the scaling is exact and leaves every coordinate within a range: both determinants keep
     *  their sign, so the original, exact only within its range, still judges them.
     *
     *  @param  points   The points, finite.
     *  @param  exponent The range, as Within's.
     *  @return          The points scaled, or none.
     */
    auto Into(std::span<const Point> points, int exponent) -> std::optional<std::vector<Point>>
    {
      int top = INT_MIN;
      for (const Point point : points)
      {
        for (const double value : {point.X, point.Y})
        {
          if (std::fpclassify(value) != FP_ZERO)
          {
            top = std::max(top, std::ilogb(value));
          }
        }
      }
      const int          shift = top == INT_MIN ? 0 : (exponent / 2) - top;
      std::vector<Point> scaled;
      for (const Point point : points)
      {
        const Point moved = {.X = std::ldexp(point.X, shift), .Y = std::ldexp(point.Y, shift)};
        const bool  exact = Exact::Same(std::ldexp(moved.X, -shift), point.X)
                         && Exact::Same(std::ldexp(moved.Y, -shift), point.Y);
        if (!exact)
        {
          return std::nullopt;
        }
        scaled.push_back(moved);
      }
      if (!Within(scaled, exponent))
      {
        return std::nullopt;
      }
      return scaled;
    }

    /**
     *  The sign three points turn by the original, scaled into its range where they lie outside.
     *
     *  @param  corners The triangle, finite.
     *  @return         The original's sign, or none where no scaling brings it into its range.
     */
    auto Reference(const Exact::Triangle &corners) -> std::optional<int>
    {
      const std::array points = {corners.First, corners.Second, corners.Third};
      const auto       scaled = Into(points, OrientRange);
      if (!scaled)
      {
        return std::nullopt;
      }
      return Original({.First = scaled->at(0), .Second = scaled->at(1), .Third = scaled->at(2)});
    }

    /**
     *  The side of a circle by the original, scaled into its range where the points lie outside.
     *
     *  @param  quad The circle and the point asked about, finite.
     *  @return      The original's sign, or none where no scaling brings them into its range.
     */
    auto Reference(const Quad &quad) -> std::optional<int>
    {
      const std::array points = {quad.Circle.First, quad.Circle.Second, quad.Circle.Third,
                                 quad.Asked};
      const auto       scaled = Into(points, CircleRange);
      if (!scaled)
      {
        return std::nullopt;
      }
      return Original(
          Quad{.Circle = {.First = scaled->at(0), .Second = scaled->at(1), .Third = scaled->at(2)},
               .Asked  = scaled->at(3)});
    }

    /**
     *  The orientation test, its adaptive stages run even where the first decides, and its laws:
     *  the same sign when the corners turn or the plane turns a half turn, the other when two swap
     *  or the plane is mirrored.
     *
     *  @param  draw The draw.
     *  @return      What broke.
     */
    auto BrokenOrient(Draw &draw) -> std::string
    {
      const Exact::Triangle t = draw.Corners();
      const Point           a = t.First;
      const Point           b = t.Second;
      const Point           c = t.Third;

      const int given    = Exact::Orient(a, b, c);
      const int adaptive = Exact::OrientAdapt(t);
      const int turned   = Exact::Orient(b, c, a);
      const int twice    = Exact::Orient(c, a, b);
      const int swapped  = Exact::Orient(b, a, c);
      const int last     = Exact::Orient(a, c, b);
      const int mirrored = Exact::Orient(Mirrored(a), Mirrored(b), Mirrored(c));
      const int negated  = Exact::Orient(Negated(a), Negated(b), Negated(c));
      const int whole    = Exact::OrientWhole(t);

      const std::array points = {a, b, c};
      const int        sign   = Finite(points) ? Reference(t).value_or(whole) : 0;
      const std::array laws   = {
          Law{.What = "as given"             , .Got = given   , .Expected = sign },
          Law{.What = "its adaptive stages"  , .Got = adaptive, .Expected = sign },
          Law{.What = "turned once"          , .Got = turned  , .Expected = sign },
          Law{.What = "turned twice"         , .Got = twice   , .Expected = sign },
          Law{.What = "the first two swapped", .Got = swapped , .Expected = -sign},
          Law{.What = "the last two swapped" , .Got = last    , .Expected = -sign},
          Law{.What = "mirrored"             , .Got = mirrored, .Expected = -sign},
          Law{.What = "turned a half turn"   , .Got = negated , .Expected = sign },
          Law{.What = "in whole numbers"     , .Got = whole   , .Expected = sign },
      };
      return Told("orient", laws, points);
    }

    /**
     *  The circle test, its adaptive stages run even where the first decides, its sum in integers
     *  wherever the places lie on one grid, and the laws of its determinant, which swaps sign when
     *  any two of its four points swap.
     *
     *  @param  draw The draw.
     *  @return      What broke.
     */
    auto BrokenCircle(Draw &draw) -> std::string
    {
      const Quad  quad = draw.Cocircular();
      const Point a    = quad.Circle.First;
      const Point b    = quad.Circle.Second;
      const Point c    = quad.Circle.Third;
      const Point d    = quad.Asked;

      const Exact::Triangle mirrored = {
          .First = Mirrored(a), .Second = Mirrored(b), .Third = Mirrored(c)};
      const Exact::Triangle negated = {
          .First = Negated(a), .Second = Negated(b), .Third = Negated(c)};

      const int given = Exact::InCircle(quad.Circle, d);
      const int adaptive =
          Exact::CircleAdapt(quad.Circle, d, Exact::Lift(quad.Circle, d).Permanent);
      const int                turned  = Exact::InCircle({.First = b, .Second = c, .Third = a}, d);
      const int                swapped = Exact::InCircle({.First = b, .Second = a, .Third = c}, d);
      const int                asked   = Exact::InCircle({.First = a, .Second = b, .Third = d}, c);
      const int                reflected = Exact::InCircle(mirrored, Mirrored(d));
      const int                halved    = Exact::InCircle(negated, Negated(d));
      const std::optional<int> grid      = Exact::OnGrid(quad.Circle, d);
      const int                whole     = Exact::CircleWhole(quad.Circle, d);

      const std::array points  = {a, b, c, d};
      const int        sign    = Finite(points) ? Reference(quad).value_or(whole) : 0;
      const int        gridded = grid.value_or(sign);
      const std::array laws    = {
          Law{.What = "as given"                      , .Got = given    , .Expected = sign },
          Law{.What = "its adaptive stages"           , .Got = adaptive , .Expected = sign },
          Law{.What = "turned once"                   , .Got = turned   , .Expected = sign },
          Law{.What = "two corners swapped"           , .Got = swapped  , .Expected = -sign},
          Law{.What = "a corner and the point swapped", .Got = asked    , .Expected = -sign},
          Law{.What = "mirrored"                      , .Got = reflected, .Expected = -sign},
          Law{.What = "turned a half turn"            , .Got = halved   , .Expected = sign },
          Law{.What = "in integers on its grid"       , .Got = gridded  , .Expected = sign },
          Law{.What = "in whole numbers"              , .Got = whole    , .Expected = sign },
      };
      return Told("circle", laws, points);
    }

    /**
     *  Points against a line, every one OrientSigns gave the sign Orient gives, wherever the input
     *  lies, and the original's, scaled into its range where they lie outside it. The points lie
     *  near the line, near its ends, near the point before, or anywhere.
     *
     *  @param  draw The draw.
     *  @return      What broke.
     */
    auto BrokenBatch(Draw &draw) -> std::string
    {
      const Exact::Triangle corners = draw.Corners();
      const Exact::Line     line    = {.First = corners.First, .Second = corners.Second};

      const unsigned low   = draw.Byte();
      const unsigned high  = draw.Byte();
      const unsigned count = ((high << 8U) | low) % (BatchMost + 1);

      std::vector<Point> points;
      points.reserve(count);
      Point before = corners.Third;
      for (unsigned index = 0; index < count; ++index)
      {
        Point point;
        switch (draw.Byte() % 4U)
        {
        case 0:
          point = draw.Free();
          break;
        case 1:
          point = draw.Near(line.First, line.Second);
          break;
        case 2:
        {
          const Point  end = (draw.Byte() & 1U) == 0 ? line.First : line.Second;
          const double x   = draw.Nudged(end.X);
          point            = {.X = x, .Y = draw.Nudged(end.Y)};
          break;
        }
        default:
        {
          const double x = draw.Nudged(before.X);
          point          = {.X = x, .Y = draw.Nudged(before.Y)};
          break;
        }
        }
        points.push_back(point);
        before = point;
      }

      std::vector<std::int8_t> signs(points.size());
      Exact::OrientSigns(line, points, signs);

      for (std::size_t index = 0; index < points.size(); ++index)
      {
        const Point           point = points[index];
        const std::array      three = {line.First, line.Second, point};
        const int             one   = Exact::Orient(line.First, line.Second, point);
        const std::int8_t     got   = signs[index];
        const Exact::Triangle t     = {.First = line.First, .Second = line.Second, .Third = point};
        if (got != static_cast<std::int8_t>(one))
        {
          return std::format("batch: point {} of {} gave {}, Orient {}, for{}", index,
                             points.size(), got, one, Text(three));
        }
        const int expected = Finite(three) ? Reference(t).value_or(one) : 0;
        if (one != expected)
        {
          return std::format("batch: point {} of {} gave {}, the original's sign {}, for{}", index,
                             points.size(), one, expected, Text(three));
        }
      }
      return {};
    }

    /**
     *  The original's sum, given copies one longer than they are, since it reads one past the end
     *  of each input.
     *
     *  @param  first  One expansion.
     *  @param  second The other.
     *  @return        fast_expansion_sum_zeroelim's sum.
     */
    auto OriginalSum(Components first, Components second) -> Components
    {
      Components sum(first.size() + second.size());
      const auto firstLength  = static_cast<int>(first.size());
      const auto secondLength = static_cast<int>(second.size());
      first.push_back(0.0);
      second.push_back(0.0);
      const int length = fast_expansion_sum_zeroelim(firstLength, first.data(), secondLength,
                                                     second.data(), sum.data());
      sum.resize(static_cast<std::size_t>(length));
      return sum;
    }

    /**
     *  The original's product, given a copy one longer than it is.
     *
     *  @param  expansion The expansion.
     *  @param  factor    The double.
     *  @return           scale_expansion_zeroelim's product.
     */
    auto OriginalScale(Components expansion, double factor) -> Components
    {
      Components product(2 * expansion.size());
      const auto length = static_cast<int>(expansion.size());
      expansion.push_back(0.0);
      product.resize(static_cast<std::size_t>(
          scale_expansion_zeroelim(length, expansion.data(), factor, product.data())));
      return product;
    }

    /**
     *  Components with their zeros dropped, as the original's routines leave them.
     *
     *  @param  components The components, smallest first.
     *  @return            The ones that are not zero, or a single zero when none is left.
     */
    auto Dropped(std::span<const double> components) -> Components
    {
      Components kept;
      for (const double component : components)
      {
        if (!Exact::IsZero(component))
        {
          kept.push_back(component);
        }
      }
      if (kept.empty())
      {
        kept.push_back(0.0);
      }
      return kept;
    }

    /**
     *  A value of two parts as components, its zeros dropped.
     *
     *  @param  two The value.
     *  @return     Its roundoff, then its rounded result.
     */
    auto Dropped(Exact::Two two) -> Components
    {
      const std::array components = {two.Small, two.Big};
      return Dropped(components);
    }

    /**
     *  The exact sums and products of two doubles, each the value the original's routines give.
     *
     *  @param  draw The draw.
     *  @return      What broke.
     */
    auto BrokenTwo(Draw &draw) -> std::string
    {
      const double x = draw.Coordinate();
      const double y = draw.Coordinate();
      if (!Within(std::array{x, y}, OrientRange))
      {
        return {};
      }

      const Exact::Two added  = Exact::TwoSum(x, y);
      const Exact::Two taken  = Exact::TwoDiff(x, y);
      const Exact::Two times  = Exact::TwoProduct(x, y);
      const Exact::Two square = Exact::Square(x);

      const std::array checks = {
          Check{.What = "TwoSum", .Held = Exact::Same(Dropped(added), OriginalSum({x}, {y}))},
          Check{.What = "TwoDiff", .Held = Exact::Same(Dropped(taken), OriginalSum({x}, {-y}))},
          Check{.What = "TwoProduct", .Held = Exact::Same(Dropped(times), OriginalScale({x}, y))},
          Check{.What = "Square", .Held = Exact::Same(Dropped(square), OriginalScale({x}, x))},
      };
      for (const Check &check : checks)
      {
        if (!check.Held)
        {
          return std::format("two: {} is not the original's for {:a} {:a}", check.What, x, y);
        }
      }
      return {};
    }
  } // namespace

  auto Broken(Target target, std::span<const std::uint8_t> bytes) -> std::string
  {
    Start();
    Draw draw{bytes};
    switch (target)
    {
    case Target::Batch:
      return BrokenBatch(draw);
    case Target::Circle:
      return BrokenCircle(draw);
    case Target::Orient:
      return BrokenOrient(draw);
    case Target::Two:
      return BrokenTwo(draw);
    }
    return {};
  }
} // namespace Geoxy::Fuzz
