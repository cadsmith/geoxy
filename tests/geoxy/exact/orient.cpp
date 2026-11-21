/***************************************************************************************************
 *  @file      tests/geoxy/exact/orient.cpp
 *  @brief     Tests of the orientation test's sign, against Shewchuk's own.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "exact/orient.hpp"

#include "../draws.hpp"
#include "exact/bounds.hpp"
#include "exact/two.hpp"
#include "original.hpp"
#include "same.hpp"

#include "geoxy/point.hpp"

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstddef>

namespace Geoxy::Exact
{
  namespace
  {
    auto OrientOriginal(Point first, Point second, Point third) -> double
    {
      std::array<double, 2> one{first.X, first.Y};
      std::array<double, 2> two{second.X, second.Y};
      std::array<double, 2> three{third.X, third.Y};
      return orient2d(one.data(), two.data(), three.data());
    }

    /**
     *  Orient gives the original's sign, and keeps it when the corners turn, and flips it when two
     *  swap.
     */
    auto OrientAgrees(Point one, Point two, Point three) -> bool
    {
      const int sign = Sign(OrientOriginal(one, two, three));

      const bool same        = Exact::Orient(one, two, three) == sign;
      const bool turnedOnce  = Exact::Orient(two, three, one) == sign;
      const bool turnedTwice = Exact::Orient(three, one, two) == sign;
      const bool swapped     = Exact::Orient(two, one, three) == -sign;
      return same && turnedOnce && turnedTwice && swapped;
    }

    /**
     *  The double the given number of steps away.
     */
    auto OrientNudged(double value, int steps) -> double
    {
      for (; steps > 0; --steps)
      {
        value = std::nextafter(value, HUGE_VAL);
      }
      for (; steps < 0; ++steps)
      {
        value = std::nextafter(value, -HUGE_VAL);
      }
      return value;
    }

    struct Scale
    {
      double Offset;
      double Spread;
    };

    constexpr std::array Scales{
        Scale{.Offset =    0, .Spread =     1},
        Scale{.Offset =    0, .Spread = 1e-20},
        Scale{.Offset =    0, .Spread =  1e20},
        Scale{.Offset =  5e6, .Spread =   1e3},
        Scale{.Offset = -5e6, .Spread =     1},
        Scale{.Offset =  1e8, .Spread =  1e-3},
        Scale{.Offset =  0.5, .Spread = 1e-12},
    };

    struct OrientFixture
    {
      OrientFixture() { exactinit(); }
    };
  } // namespace

  TEST_CASE_FIXTURE(OrientFixture, "orient: its constants are the original's")
  {
    CHECK(Bits(Epsilon)      == Bits(epsilon));
    CHECK(Bits(Splitter)     == Bits(splitter));
    CHECK(Bits(ResultBound)  == Bits(resulterrbound));
    CHECK(Bits(OrientBoundA) == Bits(ccwerrboundA));
  }

  TEST_CASE_FIXTURE(OrientFixture,
                    "orient: counterclockwise is positive, clockwise negative, collinear zero")
  {
    CHECK(Exact::Orient({0, 0}, {1, 0}, {0, 1}) == 1);
    CHECK(Exact::Orient({0, 0}, {0, 1}, {1, 0}) == -1);
    CHECK(Exact::Orient({0, 0}, {1, 1}, {2, 2}) == 0);
    CHECK(Exact::Orient({1, 1}, {1, 1}, {1, 1}) == 0);
  }

  /**
   *  Kettner et al., "Classroom examples of robustness problems in geometric computations": points
   *  a step of a double apart near the diagonal through q and r, where the plain determinant gives
   *  scattered wrong signs.
   */
  TEST_CASE_FIXTURE(OrientFixture, "orient: Kettner's grid near a line agrees with the original")
  {
    constexpr double Step = 0x1.0p-53;
    constexpr int    Side = 256;

    const Point q = {.X = 12, .Y = 12};
    const Point r = {.X = 24, .Y = 24};

    int agreed     = 0;
    int naiveWrong = 0;
    for (int i = 0; i < Side; ++i)
    {
      for (int j = 0; j < Side; ++j)
      {
        const Point p       = {.X = 0.5 + (i * Step), .Y = 0.5 + (j * Step)};
        agreed             += static_cast<int>(OrientAgrees(p, q, r));
        const double naive  = ((q.X - p.X) * (r.Y - p.Y)) - ((q.Y - p.Y) * (r.X - p.X));
        naiveWrong         += static_cast<int>(Sign(naive) != Exact::Orient(p, q, r));
      }
    }
    CHECK(agreed == Side * Side);
    CHECK(naiveWrong > 0);
  }

  TEST_CASE_FIXTURE(OrientFixture, "orient: random corners agree with the original at every scale")
  {
    Draws draws{0x6F72'6965'6E74U};

    int agreed = 0;
    int count  = 0;
    for (const Scale &scale : Scales)
    {
      for (int index = 0; index < 20'000; ++index)
      {
        const Point a  = {.X = draws.Near(scale.Offset, scale.Spread),
                          .Y = draws.Near(-scale.Offset, scale.Spread)};
        const Point b  = {.X = draws.Near(scale.Offset, scale.Spread),
                          .Y = draws.Near(-scale.Offset, scale.Spread)};
        const Point c  = {.X = draws.Near(scale.Offset, scale.Spread),
                          .Y = draws.Near(-scale.Offset, scale.Spread)};
        agreed        += static_cast<int>(OrientAgrees(a, b, c));
        ++count;
      }
    }
    CHECK(agreed == count);
  }

  /**
   *  Three corners on a line, worked in doubles, then the third moved a few steps off it: the
   *  inputs that need the later stages.
   */
  TEST_CASE_FIXTURE(OrientFixture,
                    "orient: corners steps off a line agree with the original at every scale")
  {
    Draws draws{0x6C69'6E65U};

    int agreed = 0;
    int count  = 0;
    for (const Scale &scale : Scales)
    {
      for (int index = 0; index < 20'000; ++index)
      {
        const Point  a      = {.X = draws.Near(scale.Offset, scale.Spread),
                               .Y = draws.Near(-scale.Offset, scale.Spread)};
        const Point  b      = {.X = draws.Near(scale.Offset, scale.Spread),
                               .Y = draws.Near(-scale.Offset, scale.Spread)};
        const double along  = draws.Near(0.5, 2);
        const auto   steps  = static_cast<int>(draws.Next() % 7U) - 3;
        const Point  c      = {.X = OrientNudged(a.X + (along * (b.X - a.X)), steps),
                               .Y = a.Y + (along * (b.Y - a.Y))};
        agreed             += static_cast<int>(OrientAgrees(a, b, c));
        ++count;
      }
    }
    CHECK(agreed == count);
  }

  /**
   *  Corners on the line y = x are collinear exactly, though no difference of them need be: the
   *  answer is zero after every stage.
   */
  TEST_CASE_FIXTURE(OrientFixture, "orient: corners exactly on a diagonal are collinear")
  {
    constexpr int Count = 20'000;

    Draws draws{0x6469'6167U};

    int zero   = 0;
    int agreed = 0;
    for (int index = 0; index < Count; ++index)
    {
      const double x  = draws.Near(0, 1);
      const double y  = draws.Near(0.25, 1e-3);
      const double z  = draws.Near(5e6, 1);
      const Point  a  = {.X = x, .Y = x};
      const Point  b  = {.X = y, .Y = y};
      const Point  c  = {.X = z, .Y = z};
      zero           += static_cast<int>(Exact::Orient(a, b, c) == 0);
      agreed         += static_cast<int>(OrientAgrees(a, b, c));
    }
    CHECK(zero   == Count);
    CHECK(agreed == Count);
  }

  /**
   *  The same diagonal with one corner a step off it: a determinant so small for its inputs that
   *  only the last, exact stage finds its sign.
   */
  TEST_CASE_FIXTURE(OrientFixture, "orient: corners a step off a diagonal agree with the original")
  {
    constexpr int Count = 20'000;

    Draws draws{0x7374'6570U};

    int agreed  = 0;
    int nonzero = 0;
    for (int index = 0; index < Count; ++index)
    {
      const double x  = draws.Near(0, 1);
      const double y  = draws.Near(0.25, 1e-3);
      const double z  = draws.Near(5e6, 1);
      const Point  a  = {.X = x, .Y = x};
      const Point  b  = {.X = y, .Y = y};
      const Point  c  = {.X = z, .Y = OrientNudged(z, (index % 2 == 0) ? 1 : -1)};
      agreed         += static_cast<int>(OrientAgrees(a, b, c));
      nonzero        += static_cast<int>(Exact::Orient(a, b, c) != 0);
    }
    CHECK(agreed  == Count);
    CHECK(nonzero == Count);
  }

  /**
   *  Integers to 2^40 on one line through a random direction: every difference exact, so the sign
   *  is the cross product's, and that is zero, which the fused comparison leaves to the exact one.
   */
  TEST_CASE_FIXTURE(OrientFixture, "orient: integers on a line are collinear")
  {
    constexpr int Count = 20'000;

    Draws draws{0x696E'7473U};

    const auto integer = [&draws](unsigned bits) -> double
    {
      const auto magnitude = static_cast<double>(draws.Next() >> (64U - bits));
      return (draws.Next() % 2U == 0) ? magnitude : -magnitude;
    };

    int zero   = 0;
    int agreed = 0;
    for (int index = 0; index < Count; ++index)
    {
      const Point  origin  = {.X = integer(40), .Y = integer(40)};
      const double dx      = integer(10);
      const double dy      = integer(10);
      const double far     = integer(10);
      const double near    = integer(10);
      const Point  a       = origin;
      const Point  b       = {.X = origin.X + (far * dx), .Y = origin.Y + (far * dy)};
      const Point  c       = {.X = origin.X + (near * dx), .Y = origin.Y + (near * dy)};
      zero                += static_cast<int>(Exact::Orient(a, b, c) == 0);
      agreed              += static_cast<int>(OrientAgrees(a, b, c));
    }
    CHECK(zero   == Count);
    CHECK(agreed == Count);
  }

  /**
   *  Integers whose differences, u = (base + step, base) and v = (base, base - step) with the base
   *  near 2^52, make a cross product of -step^2 beside products near 2^104: the first stage cannot
   *  tell, the differences are exact, and one fused comparison decides.
   */
  TEST_CASE_FIXTURE(OrientFixture,
                    "orient: integers just off a line take their sign from one comparison")
  {
    constexpr int Count = 20'000;

    Draws draws{0x6F66'66U};

    int agreed = 0;
    int signs  = 0;
    for (int index = 0; index < Count; ++index)
    {
      const double base   = 0x1.0p51 + static_cast<double>(draws.Next() >> 13U) - 0x1.0p21;
      const double step   = 1.0 + static_cast<double>(draws.Next() >> 44U);
      const Point  origin = {.X = static_cast<double>(draws.Next() >> 24U),
                             .Y = static_cast<double>(draws.Next() >> 24U)};
      const Point  u      = {.X = origin.X + base + step, .Y = origin.Y + base};
      const Point  v      = {.X = origin.X + base, .Y = origin.Y + base - step};
      const bool   turn   = index % 2 == 0;
      agreed += static_cast<int>(turn ? OrientAgrees(u, v, origin) : OrientAgrees(v, u, origin));
      signs  += turn ? Exact::Orient(u, v, origin) : Exact::Orient(v, u, origin);
    }
    CHECK(agreed == Count);
    CHECK(signs  == 0);
  }

  /**
   *  Every triangle over a few values, signed zeros among them: repeated corners, corners on the
   *  axes, products of zero.
   */
  TEST_CASE_FIXTURE(OrientFixture,
                    "orient: every triangle over a few values agrees with the original")
  {
    constexpr std::array  Values{0.0, -0.0, 1.0, -1.0, 0.1, 0x1.0p52};
    constexpr std::size_t Count = Values.size()
                                * Values.size()
                                * Values.size()
                                * Values.size()
                                * Values.size()
                                * Values.size();

    std::size_t agreed = 0;
    for (std::size_t code = 0; code < Count; ++code)
    {
      std::array<double, 6> pick = {};
      std::size_t           rest = code;
      for (double &value : pick)
      {
        value  = Values.at(rest % Values.size());
        rest  /= Values.size();
      }
      agreed += static_cast<std::size_t>(OrientAgrees({.X = pick[0], .Y = pick[1]},
                                                      {.X = pick[2], .Y = pick[3]},
                                                      {.X = pick[4], .Y = pick[5]}));
    }
    CHECK(agreed == Count);
  }
} // namespace Geoxy::Exact
