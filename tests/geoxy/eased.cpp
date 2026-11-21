/***************************************************************************************************
 *  @file      tests/geoxy/eased.cpp
 *  @brief     Tests of the spline that stands for one whose weights make it stiff.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "eased.hpp"

#include "basis.hpp"
#include "kernel/piece.hpp"
#include "twofold.hpp"

#include "geoxy/at.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <ranges>
#include <span>

namespace Geoxy
{
  namespace
  {
    /**
     *  A quadratic of three poles, its middle weight a power of two times its ends'.
     *
     *  @param  power The power.
     *  @return       The spline.
     */
    auto Conic(int power) -> Spline
    {
      return {.Degree  = 2,
              .Knots   = {0, 0, 0, 1, 1, 1},
              .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 1}, {.X = 2, .Y = 0}},
              .Weights = {1, std::ldexp(1.0, power), 1}};
    }
  } // namespace

  TEST_CASE("eased: a spline of even weights, of none, or that is none, needs no stand-in")
  {
    CHECK_FALSE(Eased(Conic(0), Mild).has_value());
    CHECK_FALSE(Eased(Conic(4), Mild).has_value());
    CHECK_FALSE(Eased(Conic(-40), Mild).has_value());
    Spline plain  = Conic(0);
    plain.Weights = {};
    CHECK_FALSE(Eased(plain, Mild).has_value());
    Spline broken = Conic(40);
    broken.Knots.pop_back();
    CHECK_FALSE(Eased(broken, Mild).has_value());
  }

  TEST_CASE("eased: a heavy weight's spline, cut until each piece is even, and the same curve")
  {
    for (const int power : {6, 10, 40, 300})
    {
      CAPTURE(power);
      const Spline                spline = Conic(power);
      const std::optional<Easing> easing = Eased(spline, Mild);
      REQUIRE(easing.has_value());
      REQUIRE(Valid(easing->Shape));
      const std::size_t pieces = easing->Pieces.size();
      CHECK(pieces >= static_cast<std::size_t>(power) / 2);
      CHECK(Domain(easing->Shape) == Range{.Start = 0, .End = static_cast<double>(pieces)});
      CHECK(At(easing->Shape, 0)  == Point{.X = 0, .Y = 0});
      CHECK(At(easing->Shape, static_cast<double>(pieces)) == Point{.X = 2, .Y = 0});

      /**
       *  Each piece's weights within sixteen times the lighter of its ends, and its places going on
       *  along the spline, from its start to its end.
       */
      double last = 0;
      for (std::size_t piece = 0; piece < pieces; ++piece)
      {
        const auto   from  = static_cast<std::ptrdiff_t>(2 * piece);
        const auto   three = std::ranges::subrange(easing->Shape.Weights.begin() + from,
                                                   easing->Shape.Weights.begin() + from + 3);
        const double light = std::min(three.front(), three.back());
        CHECK(std::ranges::max(three) <= 16 * light);
        const double here = Homed(*easing, static_cast<double>(piece) + 0.5);
        CHECK(here >= last);
        last = here;
      }
      CHECK(Homed(*easing, 0) == 0);
      CHECK(Homed(*easing, static_cast<double>(pieces)) == 1);
      CHECK(Homed(*easing, -3) <= 0);
    }
  }

  TEST_CASE("eased: a place of the stand-in is the spline's own, where a double can tell it")
  {
    const Spline                spline = Conic(10);
    const std::optional<Easing> easing = Eased(spline, Mild);
    REQUIRE(easing.has_value());
    const auto pieces = static_cast<double>(easing->Pieces.size());
    for (int step = 0; step <= 400; ++step)
    {
      const double along = pieces * step / 400;
      const Point  stood = At(easing->Shape, along);
      const Point  own   = At(spline, Homed(*easing, along));
      CAPTURE(along);
      CHECK(Distance(stood, own) <= 1e-12);
    }

    /**
     *  Told from both ends of its span, a place a sliver short of the span's end keeps the sliver:
     *  beside a weight of 2^300 the last piece is 2^-300 or so of the span.
     */
    const std::optional<Easing> stiff = Eased(Conic(300), Mild);
    REQUIRE(stiff.has_value());
    const auto  count = static_cast<double>(stiff->Pieces.size());
    const Sited site  = Back(*stiff, Twofold::Lifted(count - 0.5));
    CHECK(site.Span == 2);
    CHECK(site.Short.High > 0);
    CHECK(site.Short.High < 0x1p-200);
    CHECK(site.Ahead.High            == 1);
    CHECK(Homed(*stiff, count - 0.5) == 1);
  }

  TEST_CASE("eased: a paving's pieces are gentle, carry their derivatives, and tell their places")
  {
    CHECK_FALSE(Paved(Conic(0)).has_value());
    CHECK_FALSE(Paved(Conic(1)).has_value());
    CHECK_FALSE(Paved(Conic(-40)).has_value());
    Spline plain  = Conic(0);
    plain.Weights = {};
    CHECK_FALSE(Paved(plain).has_value());
    Spline broken = Conic(40);
    broken.Knots.pop_back();
    CHECK_FALSE(Paved(broken).has_value());
    Spline beyond  = Conic(0);
    beyond.Weights = {0x1p-1000, 0x1p1000, 0x1p-1000};
    CHECK_FALSE(Paved(beyond).has_value());

    for (const int power : {2, 10, 60, 1000})
    {
      CAPTURE(power);
      const Spline                spline = Conic(power);
      const std::optional<Paving> paving = Paved(spline);
      REQUIRE(paving.has_value());
      REQUIRE(paving->Treads.size() == paving->Pieces.size());
      CHECK(paving->Treads.size() > 1);
      for (const Kernel::Tread &tread : paving->Treads)
      {
        const auto   poles = std::span{tread.Shape.Poles}.first(3);
        const double light = std::min(poles.front().W, poles.back().W);
        CHECK(std::ranges::all_of(poles, [light](const Homogeneous &pole) -> bool
                                  { return pole.W <= Gentle * light; }));
        CHECK(tread.Speed.Degree == 2);
      }
      /**
       *  The pieces run on from each other, from the spline's start to its end, and a place on one
       *  is the spline's parameter there: the span's end itself for the last piece's.
       */
      CHECK(Homed(paving->Pieces.front(), 0.0) == 0.0);
      CHECK(Homed(paving->Pieces.back(), 1.0)  == 1.0);
      for (std::size_t index = 0; index + 1 < paving->Pieces.size(); ++index)
      {
        CHECK(Homed(paving->Pieces[index], 1.0) == Homed(paving->Pieces[index + 1], 0.0));
        CHECK(Distance(Kernel::End(paving->Treads[index].Shape),
                       Kernel::Start(paving->Treads[index + 1].Shape))
              <= 1e-15);
      }
    }

    /**
     *  A span of no length is passed by.
     */
    const Spline                doubled = {.Degree  = 2,
                                           .Knots   = {0, 0, 0, 1, 1, 2, 2, 2},
                                           .Poles   = {{.X = 0, .Y = 0},
                                                       {.X = 1, .Y = 0},
                                                       {.X = 1, .Y = 1},
                                                       {.X = 2, .Y = 1},
                                                       {.X = 2, .Y = 2}},
                                           .Weights = {1, 0x1p40, 1, 0x1p40, 1}};
    const std::optional<Paving> paving  = Paved(doubled);
    REQUIRE(paving.has_value());
    CHECK(paving->Pieces.front().Span == 2);
    CHECK(paving->Pieces.back().Span  == 4);
  }

  TEST_CASE("eased: a span of no length is passed by, and weights beyond doubles are not stood for")
  {
    const Spline doubled = {.Degree  = 2,
                            .Knots   = {0, 0, 0, 1, 1, 2, 2, 2},
                            .Poles   = {{.X = 0, .Y = 0},
                                        {.X = 1, .Y = 0},
                                        {.X = 1, .Y = 1},
                                        {.X = 2, .Y = 1},
                                        {.X = 2, .Y = 2}},
                            .Weights = {1, 0x1p40, 1, 0x1p40, 1}};
    REQUIRE(Valid(doubled));
    const std::optional<Easing> easing = Eased(doubled, Mild);
    REQUIRE(easing.has_value());
    CHECK(Valid(easing->Shape));
    CHECK(easing->Pieces.front().Span == 2);
    CHECK(easing->Pieces.back().Span  == 4);

    Spline beyond                    = Conic(0);
    beyond.Weights                   = {0x1p-1000, 0x1p1000, 0x1p-1000};
    const std::optional<Easing> none = Eased(beyond, Mild);
    REQUIRE(none.has_value());
    CHECK(none->Pieces.empty());
  }
} // namespace Geoxy
