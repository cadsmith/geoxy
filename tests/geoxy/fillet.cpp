/***************************************************************************************************
 *  @file      tests/geoxy/fillet.cpp
 *  @brief     Tests of FILLET: two curves of every kind rounded, their touches against mpmath's at
 *             120 digits, ties to even, every misfit; and a polyline's corners rounded.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/fillet.hpp"

#include "corner/forms.hpp"
#include "corner/guide.hpp"
#include "corner/smooth.hpp"
#include "corner/surd.hpp"
#include "expansion.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/line.hpp"
#include "geoxy/nearest.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <expected>
#include <limits>
#include <optional>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  A case of the sheet: two lines or circles, a radius and picks, and the fillet's center and
     *  touches as mpmath worked them out at 120 digits and rounded them once.
     */
    struct Sheet
    {
      std::array<int, 2>    Kinds    = {};  /**< 0 for a line, 1 for a circle.                   */
      std::array<double, 4> First    = {};  /**< A line's ends, or a circle's center and radius. */
      std::array<double, 4> Second   = {};  /**< Likewise.                                       */
      double                Radius   = 0.0; /**< The fillet's radius.                            */
      std::array<Point, 2>  Picks    = {};  /**< The picks.                                      */
      std::array<Point, 3>  Expected = {};  /**< The center, then each touch.                    */
    };

    /**
     *  A sheet's curve.
     *
     *  @param  kind    0 for a line, 1 for a circle.
     *  @param  numbers Its numbers.
     *  @return         The curve.
     */
    auto Sheeted(int kind, const std::array<double, 4> &numbers) -> Curve
    {
      if (kind == 0)
      {
        return Line{.Start = {.X = numbers[0], .Y = numbers[1]},
                    .End   = {.X = numbers[2], .Y = numbers[3]}};
      }
      return Circle({.X = numbers[0], .Y = numbers[1]}, numbers[2]);
    }

    /**
     *  Whether a point is within a distance of another.
     *
     *  @param  lhs      One point.
     *  @param  rhs      The other.
     *  @param  distance The distance.
     *  @return          Whether it is.
     */
    auto Near(Point lhs, Point rhs, double distance) -> bool
    {
      return Distance(lhs, rhs) <= distance;
    }

    /**
     *  A fillet's arc's ends.
     *
     *  @param  filleted The fillet, with an arc.
     *  @return          At its start and its end.
     */
    auto Ends(const Filleted &filleted) -> std::array<Point, 2>
    {
      return {At(*filleted.Round, 0.0), At(*filleted.Round, 1.0)};
    }

    /**
     *  A fillet's misfit, or none.
     *
     *  @param  filleted What Fillet gave.
     *  @return          Its misfit.
     */
    auto Why(const std::expected<Filleted, Misfit> &filleted) -> std::optional<Misfit>
    {
      if (filleted)
      {
        return std::nullopt;
      }
      return filleted.error();
    }

    constexpr Tolerance Exacting = {.Length = 1e-9, .Radians = 1e-9};
  } // namespace

  TEST_CASE("fillet: two lines at right angles, each picked side its own corner")
  {
    const Line across = {.Start = {.X = 0, .Y = 0}, .End = {.X = 10, .Y = 0}};
    const Line upward = {.Start = {.X = 12, .Y = -2}, .End = {.X = 12, .Y = 10}};
    const std::expected<Filleted, Misfit> rising = Fillet(
        across, upward, {.Radius = 1, .Picks = {Point{.X = 5, .Y = 0}, Point{.X = 12, .Y = 5}}});
    REQUIRE(rising);
    CHECK(std::get<Line>(rising->First)
          == Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 11, .Y = 0}});
    CHECK(std::get<Line>(rising->Second)
          == Line{.Start = {.X = 12, .Y = 1}, .End = {.X = 12, .Y = 10}});
    CHECK(rising->Round->Center == Point{.X = 11, .Y = 1});
    CHECK(rising->Touches[0]    == Point{.X = 11, .Y = 0});
    CHECK(rising->Touches[1]    == Point{.X = 12, .Y = 1});
    CHECK(rising->Round->Sweep.InDegrees() == doctest::Approx(90));
    const std::array<Point, 2> ends = Ends(*rising);
    CHECK(Near(ends[0], rising->Touches[0], 1e-15));
    CHECK(Near(ends[1], rising->Touches[1], 1e-15));

    const std::expected<Filleted, Misfit> down = Fillet(
        across, upward, {.Radius = 1, .Picks = {Point{.X = 5, .Y = 0}, Point{.X = 12, .Y = -1}}});
    REQUIRE(down);
    CHECK(std::get<Line>(down->Second)
          == Line{.Start = {.X = 12, .Y = -2}, .End = {.X = 12, .Y = -1}});
    CHECK(down->Round->Center            == Point{.X = 11, .Y = -1});
    CHECK(down->Round->Sweep.InDegrees() == doctest::Approx(-90));

    const std::expected<Filleted, Misfit> corner = Fillet(
        across, upward, {.Radius = 0, .Picks = {Point{.X = 5, .Y = 0}, Point{.X = 12, .Y = 5}}});
    REQUIRE(corner);
    CHECK(!corner->Round);
    CHECK(std::get<Line>(corner->First)
          == Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 12, .Y = 0}});
    CHECK(std::get<Line>(corner->Second)
          == Line{.Start = {.X = 12, .Y = 0}, .End = {.X = 12, .Y = 10}});

    CHECK(Why(Fillet(across, upward,
                     {.Radius = 20, .Picks = {Point{.X = 5, .Y = 0}, Point{.X = 12, .Y = 5}}}))
          == Misfit::TooLarge);
    CHECK(Why(Fillet(Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 20, .Y = 0}},
                     Line{.Start = {.X = 12, .Y = 10}, .End = {.X = 12, .Y = 1}},
                     {.Radius = 11, .Picks = {Point{.X = 5, .Y = 0}, Point{.X = 12, .Y = 5}}}))
          == Misfit::TooLarge);
  }

  TEST_CASE("fillet: lines, a line and a circle, two circles, every touch mpmath's double")
  {
    const std::array sheets = {
        Sheet{.Kinds    = {0, 0},
              .First    = {-0x1.5cc68d1d6ba8cp-3, 0x1.09beb6fa35b23p+2, -0x1.9a7f1cbeee014p-2,
                           0x1.4575fad4715d8p-3},
              .Second   = {-0x1.9a7f1cbeee014p-2, 0x1.4575fad4715d8p-3, -0x1.15383660f736cp+2,
                           -0x1.2a42236f0ec64p-1},
              .Radius   = 0x1.4ca3e37a053b3p-2,
              .Picks    = {Point{.X = -0x1.7cfba1f8e6fbap-2, .Y = 0x1.50f085d8ae7eap-1},
                           Point{.X = -0x1.c8cfd2f47f575p-1, .Y = 0x1.0f4c5384b79d6p-4}},
              .Expected = {Point{.X = -0x1.6bccb33adaeb1p-1, .Y = 0x1.b96bd200c9505p-2},
                           Point{.X = -0x1.8b831b7cd8063p-2, .Y = 0x1.a63f1debf7273p-2},
                           Point{.X = -0x1.4cf85256f0b24p-1, .Y = 0x1.ca2eb2ac0b9b3p-4}}},
        Sheet{
            .Kinds  = {1, 1},
            .First  = {-0x1.a2ff0d7102afbp+0, 0x1.0877a7c31dd3ap+1, 0x1.250a8d7bc5772p+1, 0x0p+0},
            .Second = {-0x1.d68ca6f33edfep+1, -0x1.d28fe8fbda986p+0, 0x1.c95e82c0e4b84p+1, 0x0p+0},
            .Radius = 0x1.d126c8f80a83p-4,
            .Picks  = {Point{.X = -0x1.12f66c07d1cep+0, .Y = -0x1.39631838093bp-3},
                       Point{.X = -0x1.cbdaa78addd14p-1, .Y = 0x1.b2fd3dba94c98p-2}},
            .Expected = {Point{.X = -0x1.8e5bf550be7b5p-1, .Y = 0x1.121ff834b7b55p-4},
                         Point{.X = -0x1.7769da65f428fp-1, .Y = -0x1.328c1edd402fbp-5},
                         Point{.X = -0x1.5da6fad607d8ep-1, .Y = 0x1.0814a25b34d97p-3}}},
        Sheet{.Kinds    = {0, 0},
              .First    = {0x1.312cf15e170b5p+22, 0x1.312cf6adefb58p+22, 0x1.312cfe4cf16f4p+22,
                           0x1.312d001976e7dp+22},
              .Second   = {0x1.312cfe4cf16f4p+22, 0x1.312d001976e7dp+22, 0x1.312cf92f3798cp+22,
                           0x1.312d0f426c1edp+22},
              .Radius   = 0x1.c73656389bd3dp-4,
              .Picks    = {Point{.X = 0x1.312cfcaf1622cp+22, .Y = 0x1.312cfeec06018p+22},
                           Point{.X = 0x1.312cfda93a347p+22, .Y = 0x1.312d01fe958ebp+22}},
              .Expected = {Point{.X = 0x1.312cfdc664de2p+22, .Y = 0x1.312d004440ea8p+22},
                           Point{.X = 0x1.312cfe096552cp+22, .Y = 0x1.312cffe843d6fp+22},
                           Point{.X = 0x1.312cfe3238ec7p+22, .Y = 0x1.312d0068a497bp+22}}},
        Sheet{
            .Kinds  = {0, 1},
            .First  = {0x1.312d018e1980dp+22, 0x1.312d0ffdb8aafp+22, 0x1.312cfc90be2bfp+22,
                       0x1.312d00c9fd59ap+22},
            .Second = {0x1.312cef790b4bfp+22, 0x1.312cf12666ab6p+22, 0x1.4655c4d385df7p+2, 0x0p+0},
            .Radius = 0x1.85b0889993261p-2,
            .Picks  = {Point{.X = 0x1.312cfd3069969p+22, .Y = 0x1.312d02b074c3dp+22},
                       Point{.X = 0x1.312cfafa08018p+22, .Y = 0x1.312d01fdf9359p+22}},
            .Expected = {Point{.X = 0x1.312cfbc852bd6p+22, .Y = 0x1.312d0348f4f0ap+22},
                         Point{.X = 0x1.312cfd3a94326p+22, .Y = 0x1.312d02cf6e26ap+22},
                         Point{.X = 0x1.312cfaed7735ap+22, .Y = 0x1.312d0206878eap+22}}},
        Sheet{
            .Kinds  = {1, 1},
            .First  = {0x1.312d113b571aap+22, 0x1.312d0f70777dcp+22, 0x1.3cf8b70d9b7dep+2, 0x0p+0},
            .Second = {0x1.312d0b538468cp+22, 0x1.312cf437f0ddfp+22, 0x1.26b7756f6b64dp+2, 0x0p+0},
            .Radius = 0x1.d9479d1a12a06p-3,
            .Picks  = {Point{.X = 0x1.312d003204201p+22, .Y = 0x1.312d055424bc2p+22},
                       Point{.X = 0x1.312cffb01d578p+22, .Y = 0x1.312d027ee8553p+22}},
            .Expected = {Point{.X = 0x1.312d0000b80d7p+22, .Y = 0x1.312d03e6f3fc7p+22},
                         Point{.X = 0x1.312d00c55864ap+22, .Y = 0x1.312d046a9f0b8p+22},
                         Point{.X = 0x1.312d008b3d753p+22, .Y = 0x1.312d032717c7ap+22}}},
        Sheet{.Kinds    = {0, 0},
              .First    = {-0x1.17f8baff26cd4p+0, -0x1.b1cb5cb74fb8ap+1, 0x1.7bf347e0d63bcp-2,
                           0x1.552aee2f9cb38p-2},
              .Second   = {0x1.7bf347e0d63bcp-2, 0x1.552aee2f9cb38p-2, 0x1.14306defcbf06p+2,
                           -0x1.53546d1904e14p-2},
              .Radius   = 0x1.1374febb988ep-2,
              .Picks    = {Point{.X = 0x1.80f102ca501b5p-3, .Y = -0x1.0e8b989b4d372p-3},
                           Point{.X = 0x1.ba6add6229aa8p-1, .Y = 0x1.001b02c68880ep-2}},
              .Expected = {Point{.X = 0x1.153ce86068cb8p-1, .Y = 0x1.03660881a0b03p-5},
                           Point{.X = 0x1.2a266601cbe8fp-2, .Y = 0x1.0a94138667a91p-3},
                           Point{.X = 0x1.2c1ea6cf94185p-1, .Y = 0x1.300dc044bf748p-2}}},
        Sheet{
            .Kinds  = {0, 1},
            .First  = {0x1.35dbca7877622p+1, 0x1.1fd3a6e74c203p+1, -0x1.234875efff3cp-2,
                       -0x1.654a956a2f4f9p-1},
            .Second = {-0x1.54db00b2a575fp+2, -0x1.6e8a9b7ec41eap+1, 0x1.5f29d5d3fdefcp+2, 0x0p+0},
            .Radius = 0x1.ed410ab1862c4p-3,
            .Picks  = {Point{.X = 0x1.b7e31a33c06dp-5, .Y = -0x1.516ede9286ab1p-2},
                       Point{.X = -0x1.00ec13a2f76ep-1, .Y = -0x1.fe7eabc236cap-3}},
            .Expected = {Point{.X = -0x1.c0296819272aap-3, .Y = -0x1.14b844ed3be13p-2},
                         Point{.X = -0x1.534879ba39379p-5, .Y = -0x1.bb83e8fd3269p-2},
                         Point{.X = -0x1.bbfaeb0f9caebp-2, .Y = -0x1.8462a886a1611p-2}}},
        Sheet{.Kinds  = {1, 1},
              .First  = {0x1.604e3bf02b8f4p+1, -0x1.09e1f2fcdca5bp-1, 0x1.78f7e9deeb66ap+1, 0x0p+0},
              .Second = {-0x1.21a1922ef6bd7p-1, 0x1.ec7a1263841b2p+1, 0x1.ba98cb296f73dp+1, 0x0p+0},
              .Radius = 0x1.1ea8a194bdff7p-3,
              .Picks  = {Point{.X = 0x1.4deb130546e4p-3, .Y = 0x1.c482b8be3b2cdp-1},
                         Point{.X = 0x1.cb0452f40c62ap-2, .Y = 0x1.155873af23d18p-1}},
              .Expected = {Point{.X = 0x1.7fd60a49fcd76p-3, .Y = 0x1.3b868c7a9db94p-1},
                           Point{.X = 0x1.e6ddc7b3c32e5p-5, .Y = 0x1.5889778751b1bp-1},
                           Point{.X = 0x1.c0e77dedcd96bp-3, .Y = 0x1.eb768b23084d8p-2}}},
        Sheet{.Kinds    = {0, 0},
              .First    = {0x1.312d12126126fp+22, 0x1.312cf9a082eeep+22, 0x1.312d029094f64p+22,
                           0x1.312cfd916edf6p+22},
              .Second   = {0x1.312d029094f64p+22, 0x1.312cfd916edf6p+22, 0x1.312cf487960b5p+22,
                           0x1.312cf5e2cff27p+22},
              .Radius   = 0x1.da7ef10668785p-4,
              .Picks    = {Point{.X = 0x1.312d0480ce7c5p+22, .Y = 0x1.312cfd1351615p+22},
                           Point{.X = 0x1.312d00cf7518ep+22, .Y = 0x1.312cfc9b9b01cp+22}},
              .Expected = {Point{.X = 0x1.312d02a098b79p+22, .Y = 0x1.312cfd12f7c7p+22},
                           Point{.X = 0x1.312d02bdd0dc3p+22, .Y = 0x1.312cfd85efd94p+22},
                           Point{.X = 0x1.312d0267a431dp+22, .Y = 0x1.312cfd7b063d9p+22}}},
        Sheet{
            .Kinds  = {0, 1},
            .First  = {0x1.312cfffda69abp+22, 0x1.312d0f2d0beedp+22, 0x1.312d03469ad26p+22,
                       0x1.312cff844dfb3p+22},
            .Second = {0x1.312cf78dac418p+22, 0x1.312d0d9d80668p+22, 0x1.255d0854ee1e6p+2, 0x0p+0},
            .Radius = 0x1.412c1d0dbe69fp-3,
            .Picks  = {Point{.X = 0x1.312d02dd7c4b7p+22, .Y = 0x1.312d017965b9ap+22},
                       Point{.X = 0x1.312d04bc472dp+22, .Y = 0x1.312d00deff468p+22}},
            .Expected = {Point{.X = 0x1.312d03ab210c9p+22, .Y = 0x1.312d00b3469aep+22},
                         Point{.X = 0x1.312d030df6c4ep+22, .Y = 0x1.312d00924e41dp+22},
                         Point{.X = 0x1.312d0418feee1p+22, .Y = 0x1.312d003e27c01p+22}}},
        Sheet{
            .Kinds  = {1, 1},
            .First  = {0x1.312cf5030ba6fp+22, 0x1.312cfe9d17478p+22, 0x1.dd72711366d56p+1, 0x0p+0},
            .Second = {0x1.312d167ec43d4p+22, 0x1.312d0a4734942p+22, 0x1.6a412973013adp+2, 0x0p+0},
            .Radius = 0x1.28574c3a21304p-2,
            .Picks  = {Point{.X = 0x1.312d03964b7aep+22, .Y = 0x1.312cfb6c6d0f6p+22},
                       Point{.X = 0x1.312d05168fcfdp+22, .Y = 0x1.312cfbcce680bp+22}},
            .Expected = {Point{.X = 0x1.312d04953e83bp+22, .Y = 0x1.312cfa9be907ep+22},
                         Point{.X = 0x1.312d03763e265p+22, .Y = 0x1.312cfae5b8d6ep+22},
                         Point{.X = 0x1.312d057449dcep+22, .Y = 0x1.312cfb5f060c6p+22}}},
        Sheet{.Kinds    = {0, 0},
              .First    = {-0x1.931aaf7ab6972p+1, -0x1.31ddf39a9ed12p+0, 0x1.122ad28f1675cp-2,
                           0x1.c4f328419e0ap-1},
              .Second   = {0x1.122ad28f1675cp-2, 0x1.c4f328419e0ap-1, 0x1.f4ec09de9cbcbp+0,
                           -0x1.5ee181ddf3e69p+1},
              .Radius   = 0x1.657195ff9a478p-3,
              .Picks    = {Point{.X = -0x1.466a6e7b05e02p-3, .Y = 0x1.3fdd4652c2948p-1},
                           Point{.X = 0x1.ea5b7d2c82056p-2, .Y = 0x1.b9c80494e0aafp-2}},
              .Expected = {Point{.X = 0x1.88f7460aff571p-3, .Y = 0x1.44b4473919167p-1},
                           Point{.X = 0x1.9e496aab1eb9ep-4, .Y = 0x1.910a66e3a0bc4p-1},
                           Point{.X = 0x1.667dae527c1c5p-2, .Y = 0x1.6a6fc07bb9a81p-1}}},
        Sheet{.Kinds    = {0, 1},
              .First    = {0x1.77d7172b1cc35p+1, 0x1.8b64b342c6b1cp+1, -0x1.263c4ade3e456p-2,
                           0x1.710b96e32aaa2p-1},
              .Second   = {0x1.294fd9f3ec729p-2, 0x1.53f08831da24p+1, 0x1.027affd6a92f8p+1, 0x0p+0},
              .Radius   = 0x1.01ddeb31def09p-4,
              .Picks    = {Point{.X = 0x1.d98956a2991a8p-4, .Y = 0x1.044e3ed41456ep+0},
                           Point{.X = 0x1.8ea6d52532e9ep-3, .Y = 0x1.46ff753f16198p-1}},
              .Expected = {Point{.X = -0x1.267bece068194p-3, .Y = 0x1.7f0b01dac8f45p-1},
                           Point{.X = -0x1.72d1da61aaea4p-3, .Y = 0x1.990510cda8e37p-1},
                           Point{.X = -0x1.43183ed2a634cp-3, .Y = 0x1.5f9cfe07f42a4p-1}}},
        Sheet{.Kinds = {1, 1},
              .First = {-0x1.65e73415a3ce2p-1, -0x1.67bbef02cba5cp+1, 0x1.754e0d1c9a026p+1, 0x0p+0},
              .Second = {-0x1.3f7d0228838cdp+0, 0x1.084a6d743b61ep+2, 0x1.2e0e023795293p+2, 0x0p+0},
              .Radius = 0x1.9a80ebd7a48c1p-5,
              .Picks  = {Point{.X = 0x1.093cdaa97f804p+0, .Y = -0x1.dd77567ae2f9p-2},
                         Point{.X = 0x1.0fbe7dc923abbp+0, .Y = 0x1.be15c1d01f6p-7}},
              .Expected = {Point{.X = 0x1.75bc72a0bf6bcp-1, .Y = -0x1.af8df18868b3bp-3},
                           Point{.X = 0x1.6960a0f7c9269p-1, .Y = -0x1.04bea3fc65fe1p-2},
                           Point{.X = 0x1.6b18bcf4f76ecp-1, .Y = -0x1.522b322a64b08p-3}}},
        Sheet{.Kinds    = {0, 0},
              .First    = {0x1.312d05e54f4ap+22, 0x1.312cf03d13297p+22, 0x1.312cfc9cca44dp+22,
                           0x1.312cfd45255eap+22},
              .Second   = {0x1.312cfc9cca44dp+22, 0x1.312cfd45255eap+22, 0x1.312ced2de03a9p+22,
                           0x1.312d017da6dc4p+22},
              .Radius   = 0x1.74717b98c8f95p-4,
              .Picks    = {Point{.X = 0x1.312cfdc5dae57p+22, .Y = 0x1.312cfba42318p+22},
                           Point{.X = 0x1.312cfaaeed038p+22, .Y = 0x1.312cfdcc358e5p+22}},
              .Expected = {Point{.X = 0x1.312cfc643646bp+22, .Y = 0x1.312cfcf416cbdp+22},
                           Point{.X = 0x1.312cfcb00c4bfp+22, .Y = 0x1.312cfd2a1cc6fp+22},
                           Point{.X = 0x1.312cfc7cc6335p+22, .Y = 0x1.312cfd4de6d9dp+22}}},
        Sheet{.Kinds  = {0, 1},
              .First  = {0x1.312d113df6a85p+22, 0x1.312cff2d4d448p+22, 0x1.312d013e343bbp+22,
                         0x1.312cff00e9bbdp+22},
              .Second = {0x1.312d0cbcac0d7p+22, 0x1.312cef044123cp+22, 0x1.3b0a1d3e7d48p+2, 0x0p+0},
              .Radius = 0x1.adb2a6492d896p-4,
              .Picks  = {Point{.X = 0x1.312d033e2c894p+22, .Y = 0x1.312cff06762cep+22},
                         Point{.X = 0x1.312cffafaf61cp+22, .Y = 0x1.312cfdc29cc94p+22}},
              .Expected = {Point{.X = 0x1.312d016188fe1p+22, .Y = 0x1.312cfe95dd7b2p+22},
                           Point{.X = 0x1.312d01605ef6fp+22, .Y = 0x1.312cff0148875p+22},
                           Point{.X = 0x1.312d01223add7p+22, .Y = 0x1.312cfeeca7ab3p+22}}},
        Sheet{.Kinds = {1, 1},
              .First = {0x1.312d154157849p+22, 0x1.312cf56621c8fp+22, 0x1.699b7bd065028p+2, 0x0p+0},
              .Second = {0x1.312d049fd9a65p+22, 0x1.312d1513ab91ap+22, 0x1.55e3f37450bdp+2, 0x0p+0},
              .Radius = 0x1.c1a7a54e8e4c3p-3,
              .Picks  = {Point{.X = 0x1.312d024a4b214p+22, .Y = 0x1.312d01b10a6d9p+22},
                         Point{.X = 0x1.312d034302e59p+22, .Y = 0x1.312cffc08eaeap+22}},
              .Expected = {Point{.X = 0x1.312d02b0e617ep+22, .Y = 0x1.312d00ada786fp+22},
                           Point{.X = 0x1.312d01f0c1c9ap+22, .Y = 0x1.312d012265b34p+22},
                           Point{.X = 0x1.312d029baf322p+22, .Y = 0x1.312cffcdd480ep+22}}},
        Sheet{.Kinds    = {0, 0},
              .First    = {-0x1.33bb8154a041ap+0, -0x1.be0d18af58e12p+1, 0x1.98aa67d78582p-4,
                           0x1.30954c3246844p-2},
              .Second   = {0x1.98aa67d78582p-4, 0x1.30954c3246844p-2, 0x1.e242ecc974259p+1,
                           0x1.e4aaf41e6e54ep+0},
              .Radius   = 0x1.db932015854fap-3,
              .Picks    = {Point{.X = -0x1.01e1e7ccabb18p-4, .Y = -0x1.6714ec06b65acp-3},
                           Point{.X = 0x1.1dd419c04cadp-1, .Y = 0x1.fcd81cbb34de2p-2}},
              .Expected = {Point{.X = 0x1.24f59a51560c3p-2, .Y = 0x1.008c3bfedebb5p-3},
                           Point{.X = 0x1.10794102a5637p-4, .Y = 0x1.9b54693afd5fdp-3},
                           Point{.X = 0x1.8c30863000b24p-3, .Y = 0x1.5a5195761d9d6p-2}}},
        Sheet{.Kinds  = {0, 1},
              .First  = {0x1.082260c4ace81p+2, -0x1.8e5d15c2fc58ep+0, 0x1.55ff09ce88fcap-1,
                         0x1.cf5375ae00a58p-2},
              .Second = {-0x1.5cbcdf04274a1p+0, 0x1.e150b8914ed53p+1, 0x1.f0c9c93fac544p+1, 0x0p+0},
              .Radius = 0x1.26aecfe5a5344p-2,
              .Picks  = {Point{.X = 0x1.19b0c4acb2629p+0, .Y = 0x1.9c74f82d84c8cp-3},
                         Point{.X = 0x1.12ef5cc00e9dfp+0, .Y = 0x1.7a5bfadde8ce4p-1}},
              .Expected = {Point{.X = 0x1.3411f95f25eefp+0, .Y = 0x1.e5b78a58e8961p-2},
                           Point{.X = 0x1.0f13b19240c35p+0, .Y = 0x1.cdc27e470562dp-3},
                           Point{.X = 0x1.06bb498d7225cp+0, .Y = 0x1.66fe1b31570fep-1}}},
        Sheet{.Kinds    = {0, 0},
              .First    = {0x1.312cfe13a5e3dp+22, 0x1.312d12028864p+22, 0x1.312cfd68a45ccp+22,
                           0x1.312d02061aa26p+22},
              .Second   = {0x1.312cfd68a45ccp+22, 0x1.312d02061aa26p+22, 0x1.312d095ce993fp+22,
                           0x1.312cf7639cc01p+22},
              .Radius   = 0x1.558e855f5da54p-3,
              .Picks    = {Point{.X = 0x1.312cfd7e048dap+22, .Y = 0x1.312d0405a85a9p+22},
                           Point{.X = 0x1.312cfee72d03ap+22, .Y = 0x1.312d00b1cae61p+22}},
              .Expected = {Point{.X = 0x1.312cfe16a6f2p+22, .Y = 0x1.312d024fe0ac8p+22},
                           Point{.X = 0x1.312cfd6c05cdap+22, .Y = 0x1.312d025701ee6p+22},
                           Point{.X = 0x1.312cfda524073p+22, .Y = 0x1.312d01d0486d2p+22}}},
        Sheet{
            .Kinds  = {0, 1},
            .First  = {0x1.312cf77e0a8b7p+22, 0x1.312d0ad38efbbp+22, 0x1.312d02fac3fcfp+22,
                       0x1.312cffb05b205p+22},
            .Second = {0x1.312d1092d5191p+22, 0x1.312cf6cecd7c3p+22, 0x1.03ce93c82dddap+2, 0x0p+0},
            .Radius = 0x1.306048509ef7ap-2,
            .Picks  = {Point{.X = 0x1.312d018b2ccecp+22, .Y = 0x1.312d0114c19bcp+22},
                       Point{.X = 0x1.312d01feefcp+22, .Y = 0x1.312cfdf5e295p+22}},
            .Expected = {Point{.X = 0x1.312d0177d82adp+22, .Y = 0x1.312cff7f8c047p+22},
                         Point{.X = 0x1.312d024bb809ap+22, .Y = 0x1.312d005a12e83p+22},
                         Point{.X = 0x1.312d027fac7ebp+22, .Y = 0x1.312cfee7c2e9cp+22}}},
        Sheet{.Kinds = {1, 1},
              .First = {0x1.312d0b98ef85fp+22, 0x1.312cfd02dfd6bp+22, 0x1.b3ed44993537ap+1, 0x0p+0},
              .Second = {0x1.312ce928b7938p+22, 0x1.312cfc5fecefp+22, 0x1.702144089c459p+2, 0x0p+0},
              .Radius = 0x1.46d738614441ap-3,
              .Picks  = {Point{.X = 0x1.312d0039e8dd8p+22, .Y = 0x1.312d048349b77p+22},
                         Point{.X = 0x1.312cfe9e8b116p+22, .Y = 0x1.312d04abd0a01p+22}},
              .Expected = {Point{.X = 0x1.312cff66b2d3p+22, .Y = 0x1.312d0467063aap+22},
                           Point{.X = 0x1.312cfff275219p+22, .Y = 0x1.312d041253b95p+22},
                           Point{.X = 0x1.312cfeccfbd5fp+22, .Y = 0x1.312d042f8b7bcp+22}}},
    };
    for (const Sheet &sheet : sheets)
    {
      const std::expected<Filleted, Misfit> filleted =
          Fillet(Sheeted(sheet.Kinds[0], sheet.First), Sheeted(sheet.Kinds[1], sheet.Second),
                 {.Radius = sheet.Radius, .Picks = sheet.Picks});
      REQUIRE(filleted);
      CHECK(filleted->Round->Center == sheet.Expected[0]);
      CHECK(filleted->Touches[0]    == sheet.Expected[1]);
      CHECK(filleted->Touches[1]    == sheet.Expected[2]);
    }
  }

  TEST_CASE("fillet: a touch halfway between two doubles rounds to the even one")
  {
    const double huge   = 0x1p53;
    const Line   across = {.Start = {.X = 0, .Y = 0}, .End = {.X = huge + 4, .Y = 0}};
    const Line   upward = {.Start = {.X = huge + 4, .Y = 0}, .End = {.X = huge + 4, .Y = 8}};
    const std::expected<Filleted, Misfit> filleted =
        Fillet(across, upward,
               {.Radius = 3, .Picks = {Point{.X = 1, .Y = 0}, Point{.X = huge + 4, .Y = 4}}});
    REQUIRE(filleted);
    CHECK(filleted->Touches[0]    == Point{.X = huge, .Y = 0});
    CHECK(filleted->Touches[1]    == Point{.X = huge + 4, .Y = 3});
    CHECK(filleted->Round->Center == Point{.X = huge, .Y = 3});
  }

  TEST_CASE("fillet: a radius that just fits between a line and a circle, a double root")
  {
    const Line across = {.Start = {.X = -5, .Y = 0}, .End = {.X = 5, .Y = 0}};
    const Arc  above  = Circle({.X = 0, .Y = 2}, 1);
    const std::expected<Filleted, Misfit> filleted = Fillet(
        across, above, {.Radius = 0.5, .Picks = {Point{.X = -4, .Y = 0}, Point{.X = -1, .Y = 2}}});
    REQUIRE(filleted);
    CHECK(filleted->Round->Center == Point{.X = 0, .Y = 0.5});
    CHECK(filleted->Touches[0]    == Point{.X = 0, .Y = 0});
    CHECK(filleted->Touches[1]    == Point{.X = 0, .Y = 1});
    CHECK(std::holds_alternative<Arc>(filleted->Second));
  }

  TEST_CASE("fillet: parallel lines joined by a half circle, as AutoCAD joins them")
  {
    const Line low  = {.Start = {.X = 0, .Y = 0}, .End = {.X = 10, .Y = 0}};
    const Line high = {.Start = {.X = 0, .Y = 3}, .End = {.X = 8, .Y = 3}};
    const std::expected<Filleted, Misfit> right =
        Fillet(low, high, {.Radius = 1, .Picks = {Point{.X = 9, .Y = 0}, Point{.X = 5, .Y = 3}}});
    REQUIRE(right);
    CHECK(std::get<Line>(right->First) == low);
    CHECK(std::get<Line>(right->Second)
          == Line{.Start = {.X = 0, .Y = 3}, .End = {.X = 10, .Y = 3}});
    CHECK(right->Round->Center            == Point{.X = 10, .Y = 1.5});
    CHECK(right->Round->Radius            == 1.5);
    CHECK(right->Round->Sweep.InDegrees() == doctest::Approx(180));

    const std::expected<Filleted, Misfit> left =
        Fillet(low, Line{.Start = {.X = 8, .Y = 3}, .End = {.X = 1, .Y = 3}},
               {.Radius = 1, .Picks = {Point{.X = 1, .Y = 0}, Point{.X = 5, .Y = 3}}});
    REQUIRE(left);
    CHECK(left->Round->Center == Point{.X = 0, .Y = 1.5});
    CHECK(std::get<Line>(left->Second) == Line{.Start = {.X = 8, .Y = 3}, .End = {.X = 0, .Y = 3}});
    CHECK(left->Round->Sweep.InDegrees() == doctest::Approx(-180));

    const std::expected<Filleted, Misfit> ray =
        Fillet(Ray{.Start = {.X = 2, .Y = 0}, .Way = {.X = 1, .Y = 0}},
               Ray{.Start = {.X = 5, .Y = -3}, .Way = {.X = 1, .Y = 0}},
               {.Radius = 1, .Picks = {Point{.X = 9, .Y = 0}, Point{.X = 9, .Y = -3}}});
    REQUIRE(ray);
    CHECK(ray->Touches[1]            == Point{.X = 2, .Y = -3});
    CHECK(std::get<Ray>(ray->Second) == Ray{.Start = {.X = 2, .Y = -3}, .Way = {.X = 1, .Y = 0}});

    const Polyline                        path = {.Vertices = {{.Position = {.X = 0, .Y = 3}},
                                                               {.Position = {.X = 6, .Y = 3}},
                                                               {.Position = {.X = 6, .Y = 9}}}};
    const std::expected<Filleted, Misfit> polyline =
        Fillet(low, path, {.Radius = 1, .Picks = {Point{.X = 1, .Y = 0}, Point{.X = 3, .Y = 3}}});
    REQUIRE(polyline);
    CHECK(polyline->Touches[1] == Point{.X = 0, .Y = 3});
    const std::expected<Filleted, Misfit> segment =
        Fillet(path, low, {.Radius = 1, .Picks = {Point{.X = 5, .Y = 3}, Point{.X = 3, .Y = 0}}});
    REQUIRE(segment);
    CHECK(segment->Touches[0] == Point{.X = 6, .Y = 3});

    CHECK(Why(Fillet(Xline{.Through = {.X = 0, .Y = 0}, .Way = {.X = 1, .Y = 0}}, high,
                     {.Radius = 1, .Picks = {Point{.X = 1, .Y = 0}, Point{.X = 1, .Y = 3}}}))
          == Misfit::Parallel);
    CHECK(Why(Fillet(low, Line{.Start = {.X = 20, .Y = 0}, .End = {.X = 30, .Y = 0}},
                     {.Radius = 1, .Picks = {Point{.X = 1, .Y = 0}, Point{.X = 25, .Y = 0}}}))
          == Misfit::Parallel);
    CHECK(Why(Fillet(low, Ray{.Start = {.X = 20, .Y = 3}, .Way = {.X = 1, .Y = 0}},
                     {.Radius = 1, .Picks = {Point{.X = 9, .Y = 0}, Point{.X = 25, .Y = 3}}}))
          == Misfit::TooLarge);
  }

  TEST_CASE("fillet: rays and endless lines cut to rays and lines")
  {
    const Ray   ray   = {.Start = {.X = 0, .Y = 0}, .Way = {.X = 1, .Y = 0}};
    const Xline xline = {.Through = {.X = 5, .Y = 0}, .Way = {.X = 0, .Y = 2}};
    const std::expected<Filleted, Misfit> ahead =
        Fillet(ray, xline, {.Radius = 1, .Picks = {Point{.X = 9, .Y = 0}, Point{.X = 5, .Y = 7}}});
    REQUIRE(ahead);
    CHECK(std::get<Ray>(ahead->First)  == Ray{.Start = {.X = 6, .Y = 0}, .Way = {.X = 1, .Y = 0}});
    CHECK(std::get<Ray>(ahead->Second) == Ray{.Start = {.X = 5, .Y = 1}, .Way = {.X = 0, .Y = 2}});
    const std::expected<Filleted, Misfit> back =
        Fillet(ray, xline, {.Radius = 1, .Picks = {Point{.X = 1, .Y = 0}, Point{.X = 5, .Y = -7}}});
    REQUIRE(back);
    CHECK(std::get<Line>(back->First) == Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 4, .Y = 0}});
    CHECK(std::get<Ray>(back->Second)
          == Ray{.Start = {.X = 5, .Y = -1}, .Way = {.X = -0.0, .Y = -2}});
    CHECK(Why(Fillet(ray, xline,
                     {.Radius = 6, .Picks = {Point{.X = 1, .Y = 0}, Point{.X = 5, .Y = -7}}}))
          == Misfit::TooLarge);
  }

  TEST_CASE("fillet: arcs cut back or run on round their circle, circles kept whole")
  {
    const Line across = {.Start = {.X = 0, .Y = 0}, .End = {.X = 10, .Y = 0}};
    const Arc  bottom = {
        .Center = {.X = 5, .Y = 5}, .Radius = 3, .Start = Degrees(180), .Sweep = Degrees(180)};
    const std::expected<Filleted, Misfit> apart = Fillet(
        across, bottom, {.Radius = 1, .Picks = {Point{.X = 1, .Y = 0}, Point{.X = 2, .Y = 5}}});
    REQUIRE(apart);
    CHECK(apart->Touches[0] == Point{.X = 5, .Y = 0});
    CHECK(apart->Touches[1] == Point{.X = 5, .Y = 2});
    CHECK(std::get<Arc>(apart->Second).Sweep.InDegrees() == doctest::Approx(90));

    const Line through = {.Start = {.X = 0, .Y = 4}, .End = {.X = 10, .Y = 4}};
    const Arc  right   = {
        .Center = {.X = 5, .Y = 4}, .Radius = 3, .Start = Degrees(-60), .Sweep = Degrees(120)};
    const std::expected<Filleted, Misfit> crossed =
        Fillet(through, right, {.Radius = 0.5, .Picks = {Point{.X = 9.5, .Y = 4}, At(right, 0.9)}});
    REQUIRE(crossed);
    CHECK(std::abs(Distance(crossed->Round->Center, right.Center) - 3.5) < 1e-14);
    CHECK(std::get<Arc>(crossed->Second).Sweep.InDegrees() < 60.0);
    const std::expected<Filleted, Misfit> lower =
        Fillet(through, right, {.Radius = 0.5, .Picks = {Point{.X = 9.5, .Y = 4}, At(right, 0.1)}});
    REQUIRE(lower);
    CHECK(std::abs(std::get<Arc>(lower->Second).Start.InDegrees() + 60) < 1e-12);

    const Arc brief = {
        .Center = {.X = 5, .Y = 4}, .Radius = 3, .Start = Degrees(30), .Sweep = Degrees(40)};
    const std::expected<Filleted, Misfit> run =
        Fillet(through, brief, {.Radius = 0.5, .Picks = {Point{.X = 9.5, .Y = 4}, At(brief, 0.5)}});
    REQUIRE(run);
    CHECK(std::get<Arc>(run->Second).Start.InDegrees() < 30.0);
    const Arc clockwise = {
        .Center = {.X = 5, .Y = 4}, .Radius = 3, .Start = Degrees(70), .Sweep = Degrees(-40)};
    const std::expected<Filleted, Misfit> backward =
        Fillet(through, clockwise,
               {.Radius = 0.5, .Picks = {Point{.X = 9.5, .Y = 4}, At(clockwise, 0.5)}});
    REQUIRE(backward);
    CHECK(std::get<Arc>(backward->Second).Sweep.InDegrees() < -40.0);
    const Arc beyond = {
        .Center = {.X = 5, .Y = 4}, .Radius = 3, .Start = Degrees(181), .Sweep = Degrees(1)};
    CHECK(Why(Fillet(Line{.Start = {.X = 0, .Y = 4}, .End = {.X = 6, .Y = 4}}, beyond,
                     {.Radius = 0.5, .Picks = {Point{.X = 1, .Y = 4}, At(beyond, 0.5)}}))
          == Misfit::TooLarge);

    const Arc                             circle = Circle({.X = 5, .Y = 4}, 3);
    const std::expected<Filleted, Misfit> whole =
        Fillet(through, circle,
               {.Radius = 0.5, .Picks = {Point{.X = 9.5, .Y = 4}, Point{.X = 8, .Y = 5}}});
    REQUIRE(whole);
    CHECK(std::get<Arc>(whole->Second) == circle);

    const Arc                             other   = Circle({.X = 10, .Y = 4}, 3);
    const std::expected<Filleted, Misfit> circles = Fillet(
        circle, other, {.Radius = 0.5, .Picks = {Point{.X = 5, .Y = 7}, Point{.X = 10, .Y = 7}}});
    REQUIRE(circles);
    CHECK(std::abs(Distance(circles->Round->Center, circle.Center) - 3.5) < 1e-14);
    CHECK(std::abs(Distance(circles->Round->Center, other.Center) - 3.5) < 1e-14);
    const std::expected<Filleted, Misfit> inner =
        Fillet(circle, Circle({.X = 5, .Y = 4}, 1),
               {.Radius = 0.5, .Picks = {Point{.X = 5, .Y = 7}, Point{.X = 5, .Y = 5}}});
    CHECK(Why(inner) == Misfit::Apart);
    const std::expected<Filleted, Misfit> enclosing =
        Fillet(Circle({.X = 0, .Y = 0}, 1), Circle({.X = 5, .Y = 0}, 1),
               {.Radius = 4, .Picks = {Point{.X = 0, .Y = 1}, Point{.X = 5, .Y = 1}}});
    REQUIRE(enclosing);
    CHECK(std::abs(Distance(enclosing->Round->Center, Point{.X = 0, .Y = 0}) - 5) < 1e-14);
    const std::expected<Filleted, Misfit> same =
        Fillet(Circle({.X = 0, .Y = 0}, 1), Circle({.X = 0, .Y = 0}, 2),
               {.Radius = 0.5, .Picks = {Point{.X = 0, .Y = 1}, Point{.X = 0, .Y = 2}}});
    CHECK(Why(same) == Misfit::Apart);
  }

  TEST_CASE("fillet: a polyline by its segment picked, cut, run on, or kept whole")
  {
    const Polyline path    = {.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                           {.Position = {.X = 5, .Y = 0}, .Bulge = 0.5},
                                           {.Position = {.X = 10, .Y = 0}}}};
    const Line     upright = {.Start = {.X = 2, .Y = -5}, .End = {.X = 2, .Y = 5}};
    const std::expected<Filleted, Misfit> start = Fillet(
        path, upright, {.Radius = 1, .Picks = {Point{.X = 4, .Y = 0}, Point{.X = 2, .Y = 3}}});
    REQUIRE(start);
    const auto &kept = std::get<Polyline>(start->First);
    REQUIRE(kept.Vertices.size()    == 3);
    CHECK(kept.Vertices[0].Position == Point{.X = 3, .Y = 0});
    const std::expected<Filleted, Misfit> end = Fillet(
        path, upright, {.Radius = 1, .Picks = {Point{.X = 1, .Y = 0}, Point{.X = 2, .Y = 3}}});
    REQUIRE(end);
    CHECK(std::get<Polyline>(end->First).Vertices.size() == 2);

    const Line late = {.Start = {.X = 8, .Y = -5}, .End = {.X = 8, .Y = 5}};
    const std::expected<Filleted, Misfit> bow = Fillet(
        path, late, {.Radius = 0.5, .Picks = {Point{.X = 6, .Y = -1}, Point{.X = 8, .Y = -4}}});
    REQUIRE(bow);
    CHECK(std::get<Polyline>(bow->First).Vertices.size() == 3);
    const std::expected<Filleted, Misfit> past =
        Fillet(path, Line{.Start = {.X = 12, .Y = -5}, .End = {.X = 12, .Y = 5}},
               {.Radius = 0.5, .Picks = {Point{.X = 1, .Y = 0}, Point{.X = 12, .Y = 4}}});
    REQUIRE(past);
    CHECK(std::get<Polyline>(past->First).Vertices.size()      == 2);
    CHECK(std::get<Polyline>(past->First).Vertices[1].Position == Point{.X = 11.5, .Y = 0});
    CHECK(Why(Fillet(path, Line{.Start = {.X = 3, .Y = -5}, .End = {.X = 3, .Y = 5}},
                     {.Radius = 2.5, .Picks = {Point{.X = 4, .Y = 0}, Point{.X = 3, .Y = 4}}}))
          == Misfit::TooLarge);
    CHECK(Why(Fillet(path, Line{.Start = {.X = 3, .Y = -5}, .End = {.X = 3, .Y = 5}},
                     {.Radius = 3.5, .Picks = {Point{.X = 1, .Y = 0}, Point{.X = 3, .Y = 4}}}))
          == Misfit::TooLarge);

    Polyline loop                                = path;
    loop.Closed                                  = true;
    const std::expected<Filleted, Misfit> closed = Fillet(
        loop, upright, {.Radius = 1, .Picks = {Point{.X = 4, .Y = 0}, Point{.X = 2, .Y = 3}}});
    REQUIRE(closed);
    CHECK(std::get<Polyline>(closed->First) == loop);
  }

  TEST_CASE("fillet: ellipses, splines and Bézier curves, by Newton's method on their offsets")
  {
    const Ellipse oval = {.Map   = {.X = {.X = 4, .Y = 0}, .Y = {.X = 0, .Y = 2}},
                          .Start = Degrees(0),
                          .Sweep = Degrees(360)};
    const Ellipse egg  = {
        .Map   = {.X = {.X = 3, .Y = 0}, .Y = {.X = 0, .Y = 1.5}, .Origin = {.X = 6, .Y = 0}},
        .Start = Degrees(0),
        .Sweep = Degrees(360)};
    const std::expected<Filleted, Misfit> ovals = Fillet(
        oval, egg, {.Radius = 0.5, .Picks = {Point{.X = 3, .Y = 1.3}, Point{.X = 4, .Y = 1}}});
    REQUIRE(ovals);
    CHECK(std::get<Ellipse>(ovals->First) == oval);
    CHECK(std::abs(Distance(ovals->Touches[0], ovals->Round->Center) - 0.5) < 1e-12);
    CHECK(std::abs(Distance(ovals->Touches[1], ovals->Round->Center) - 0.5) < 1e-12);
    CHECK(Distance(At(oval, Nearest(oval, ovals->Touches[0]).Along), ovals->Touches[0]) < 1e-12);

    const Ellipse half  = {.Map = oval.Map, .Start = Degrees(0), .Sweep = Degrees(180)};
    const Line    floor = {.Start = {.X = -6, .Y = 1}, .End = {.X = 6, .Y = 1}};
    const std::expected<Filleted, Misfit> cut = Fillet(
        half, floor, {.Radius = 0.3, .Picks = {Point{.X = 0, .Y = 2}, Point{.X = 5, .Y = 1}}});
    REQUIRE(cut);
    CHECK(std::get<Ellipse>(cut->First).Sweep.InDegrees() < 180.0);
    const std::expected<Filleted, Misfit> corner =
        Fillet(half, floor, {.Radius = 0, .Picks = {Point{.X = 0, .Y = 2}, Point{.X = 5, .Y = 1}}});
    REQUIRE(corner);
    CHECK(std::abs(corner->Touches[0].Y - 1) < 1e-12);

    const Arc                             ring  = Circle({.X = 0, .Y = 5}, 2);
    const std::expected<Filleted, Misfit> round = Fillet(
        ring, oval, {.Radius = 0.75, .Picks = {Point{.X = -2, .Y = 5}, Point{.X = -3, .Y = 1.3}}});
    REQUIRE(round);
    CHECK(std::abs(Distance(round->Round->Center, ring.Center) - 2.75) < 1e-12);

    const Bezier::Cubic                   arch  = {.Start = {.X = 0, .Y = 0},
                                                   .Out   = {.X = 2, .Y = 4},
                                                   .In    = {.X = 6, .Y = 4},
                                                   .End   = {.X = 8, .Y = 0}};
    const std::expected<Filleted, Misfit> cubic = Fillet(
        arch, floor, {.Radius = 0.3, .Picks = {Point{.X = 4, .Y = 3}, Point{.X = -2, .Y = 1}}});
    REQUIRE(cubic);
    CHECK(std::holds_alternative<Bezier::Cubic>(cubic->First));
    CHECK(std::abs(cubic->Round->Center.Y - 1.3) < 1e-12);
    const Bezier::Quadratic hump = {
        .Start = {.X = 0, .Y = 0}, .Control = {.X = 4, .Y = 6}, .End = {.X = 8, .Y = 0}};
    const std::expected<Filleted, Misfit> quadratic = Fillet(
        floor, hump, {.Radius = 0.3, .Picks = {Point{.X = 7, .Y = 1}, Point{.X = 6, .Y = 2}}});
    REQUIRE(quadratic);
    CHECK(std::holds_alternative<Bezier::Quadratic>(quadratic->Second));

    const Spline                          wave = {.Degree = 3,
                                                  .Knots  = {0, 0, 0, 0, 1, 2, 2, 2, 2},
                                                  .Poles  = {{.X = 0, .Y = 0},
                                                             {.X = 2, .Y =  3},
                                                             {.X = 4, .Y = -3},
                                                             {.X = 6, .Y =  3},
                                                             {.X = 8, .Y = 0}}};
    const std::expected<Filleted, Misfit> spline =
        Fillet(wave, Line{.Start = {.X = 0, .Y = 0.5}, .End = {.X = 1, .Y = 0.5}},
               {.Radius = 0.2, .Picks = {At(wave, 0.5), Point{.X = -1, .Y = 0.5}}});
    REQUIRE(spline);
    CHECK(std::holds_alternative<Spline>(spline->First));
    CHECK(std::abs(spline->Round->Center.Y - 0.7) < 1e-12);
    const Spline loop = {
        .Degree = 2,
        .Knots  = {0, 0, 0, 1, 2, 2, 2},
        .Poles  = {{.X = 0, .Y = 0}, {.X = 4, .Y = 4}, {.X = 8, .Y = 0}, {.X = 0, .Y = 0}}};
    const std::expected<Filleted, Misfit> closed =
        Fillet(loop, Line{.Start = {.X = -5, .Y = 1}, .End = {.X = 15, .Y = 1}},
               {.Radius = 0.2, .Picks = {At(loop, 0.3), Point{.X = -4, .Y = 1}}});
    REQUIRE(closed);
    CHECK(std::get<Spline>(closed->First) == loop);
    const std::expected<Filleted, Misfit> waves =
        Fillet(wave, arch, {.Radius = 0.2, .Picks = {At(wave, 0.4), At(arch, 0.1)}});
    CHECK(waves.has_value());
  }

  TEST_CASE("fillet: curves that never meet are bridged by the circle nearest the picks")
  {
    const Arc left  = Circle({.X = 0, .Y = 0}, 2);
    const Arc right = {
        .Center = {.X = 7, .Y = 0}, .Radius = 2, .Start = Degrees(90), .Sweep = Degrees(180)};
    const std::expected<Filleted, Misfit> top =
        Fillet(left, right, {.Radius = 2, .Picks = {Point{.X = 0, .Y = 2}, At(right, 0.2)}});
    REQUIRE(top);
    CHECK(top->Round->Center.Y > 0.0);
    const Line far = {.Start = {.X = -5, .Y = -4}, .End = {.X = 5, .Y = -4}};
    const std::expected<Filleted, Misfit> under = Fillet(
        far, left, {.Radius = 1, .Picks = {Point{.X = -4, .Y = -4}, Point{.X = -2, .Y = 0}}});
    REQUIRE(under);
    CHECK(std::abs(under->Round->Center.Y + 3) < 1e-14);
    CHECK(Why(Fillet(far, left,
                     {.Radius = 0.1, .Picks = {Point{.X = -4, .Y = -4}, Point{.X = -2, .Y = 0}}}))
          == Misfit::Apart);
    CHECK(Why(Fillet(far, left,
                     {.Radius = 0, .Picks = {Point{.X = -4, .Y = -4}, Point{.X = -2, .Y = 0}}}))
          == Misfit::Apart);
    const Ellipse oval = {.Map   = {.X = {.X = 4, .Y = 0}, .Y = {.X = 0, .Y = 2}},
                          .Start = Degrees(0),
                          .Sweep = Degrees(360)};
    const std::expected<Filleted, Misfit> smooth = Fillet(
        far, oval, {.Radius = 1, .Picks = {Point{.X = -4, .Y = -4}, Point{.X = -4, .Y = 0}}});
    REQUIRE(smooth);
    CHECK(std::abs(smooth->Round->Center.Y + 3) < 1e-12);
  }

  TEST_CASE("fillet: misfits for numbers not finite and curves of no way")
  {
    const Line                 across = {.Start = {.X = 0, .Y = 0}, .End = {.X = 10, .Y = 0}};
    const Line                 upward = {.Start = {.X = 12, .Y = -2}, .End = {.X = 12, .Y = 10}};
    const double               nan    = std::numeric_limits<double>::quiet_NaN();
    const std::array<Point, 2> picks  = {Point{.X = 5, .Y = 0}, Point{.X = 12, .Y = 5}};
    CHECK(Why(Fillet(across, upward, {.Radius = nan, .Picks = picks})) == Misfit::NotFinite);
    CHECK(Why(Fillet(across, upward, {.Radius = 1, .Picks = {Point{.X = nan, .Y = 0}, picks[1]}}))
          == Misfit::NotFinite);
    CHECK(Why(Fillet(across, upward, {.Radius = 1, .Picks = {picks[0], Point{.X = 0, .Y = nan}}}))
          == Misfit::NotFinite);
    CHECK(Why(Fillet(Line{.Start = {.X = nan, .Y = 0}}, upward, {.Radius = 1, .Picks = picks}))
          == Misfit::NotFinite);
    CHECK(Why(Fillet(across, Line{.End = {.X = nan, .Y = 0}}, {.Radius = 1, .Picks = picks}))
          == Misfit::NotFinite);
    CHECK(Why(Fillet(across, upward, {.Radius = -1, .Picks = picks})) == Misfit::Degenerate);
    const std::array<Curve, 9> flat = {
        Line{.Start = {.X = 1, .Y = 1}, .End = {.X = 1, .Y = 1}},
        Ray{.Start = {.X = 1, .Y = 1}},
        Xline{.Through = {.X = 1, .Y = 1}},
        Arc{.Center = {.X = 1, .Y = 1}, .Radius = 0, .Sweep = Degrees(90)},
        Arc{.Center = {.X = 1, .Y = 1}, .Radius = 1},
        Polyline{},
        Polyline{.Vertices = {{.Position = {.X = 1, .Y = 1}}, {.Position = {.X = 1, .Y = 1}}}},
        Ellipse{.Map = {.X = {.X = 1, .Y = 0}, .Y = {.X = 2, .Y = 0}}, .Sweep = Degrees(90)},
        Spline{.Degree = 3, .Knots = {0, 1}, .Poles = {{.X = 0, .Y = 0}}},
    };
    for (const Curve &curve : flat)
    {
      CHECK(Why(Fillet(curve, upward, {.Radius = 1, .Picks = picks})) == Misfit::Degenerate);
      CHECK(Why(Fillet(across, curve, {.Radius = 1, .Picks = picks})) == Misfit::Degenerate);
    }
    CHECK(Why(Fillet(Ellipse{.Map = {.X = {.X = 1, .Y = 0}, .Y = {.X = 0, .Y = 2}}}, upward,
                     {.Radius = 1, .Picks = picks}))
          == Misfit::Degenerate);
  }

  TEST_CASE("fillet: a polyline's corners rounded where the radius fits, the rest reported")
  {
    const Polyline                        box = {.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                                              {.Position = {.X = 10, .Y = 0}},
                                                              {.Position = {.X = 10, .Y = 4}},
                                                              {.Position = {.X = 0, .Y = 4}}},
                                                 .Closed   = true};
    const std::expected<Reshaped, Misfit> one = Fillet(box, 1, Exacting);
    REQUIRE(one);
    CHECK(one->Left.empty());
    REQUIRE(one->Shape.Vertices.size() == 8);
    CHECK(one->Shape.Closed);
    CHECK(one->Shape.Vertices[0].Position == Point{.X = 1, .Y = 0});
    CHECK(one->Shape.Vertices[1].Position == Point{.X = 9, .Y = 0});
    CHECK(one->Shape.Vertices[1].Bulge    == doctest::Approx(std::tan(Pi / 8)));
    CHECK(one->Shape.Vertices[2].Position == Point{.X = 10, .Y = 1});

    const std::expected<Reshaped, Misfit> two = Fillet(box, 2, Exacting);
    REQUIRE(two);
    CHECK(two->Left.empty());
    CHECK(two->Shape.Vertices.size() == 6);

    const std::expected<Reshaped, Misfit> three = Fillet(box, 3, Exacting);
    REQUIRE(three);
    REQUIRE(three->Left.size()  == 2);
    CHECK(three->Left[0].Vertex == 2);
    CHECK(three->Left[0].Why    == Misfit::TooLarge);
    CHECK(three->Left[1].Vertex == 3);

    const Polyline open = {.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                        {.Position = {.X = 10, .Y = 0}},
                                        {.Position = {.X = 10, .Y = 4}, .Bulge = 0.4},
                                        {.Position = {.X = 0, .Y = 4}},
                                        {.Position = {.X = -5, .Y = 4}},
                                        {.Position = {.X = 5, .Y = 4}},
                                        {.Position = {.X = 5, .Y = 4}},
                                        {.Position = {.X = 5, .Y = 9}}}};
    const std::expected<Reshaped, Misfit> mixed = Fillet(open, 1, Exacting);
    REQUIRE(mixed);
    REQUIRE(mixed->Left.size()  == 3);
    CHECK(mixed->Left[0].Vertex == 4);
    CHECK(mixed->Left[0].Why    == Misfit::Parallel);
    CHECK(mixed->Left[1].Why    == Misfit::Degenerate);
    CHECK(mixed->Left[2].Why    == Misfit::Degenerate);
    CHECK(mixed->Shape.Vertices.front().Position == Point{.X = 0, .Y = 0});
    CHECK(mixed->Shape.Vertices.back().Position  == Point{.X = 5, .Y = 9});

    CHECK(!Fillet(box, -1, Exacting).has_value());
    CHECK(Fillet(box, std::numeric_limits<double>::infinity(), Exacting).error()
          == Misfit::NotFinite);
    CHECK(Fillet(Polyline{}, 1, Exacting)->Shape.Vertices.empty());
  }

  TEST_CASE("fillet: nested roots' signs and the double nearest one, told exactly")
  {
    using Corner::Big;
    using Corner::Form;
    const auto held = [](double value) -> Big
    {
      return Corner::Held<Big>(value);
    };
    CHECK(Corner::Signed(held(0), held(1), held(2))  == 1);
    CHECK(Corner::Signed(held(1), held(1), held(2))  == 1);
    CHECK(Corner::Signed(held(3), held(-1), held(4)) == 1);
    CHECK(Corner::Signed(held(1), held(-1), held(4)) == -1);
    CHECK(Corner::Signed(held(2), held(-1), held(4)) == 0);
    CHECK(Corner::Signed(held(-2), held(0), held(4)) == -1);

    const Form<Big> zero  = {.Below = held(1)};
    const Form<Big> outer = {.Root = held(4), .Outer = held(1), .Base = held(9), .Below = held(1)};
    const Form<Big> mixed = {
        .Whole = held(5), .Root = held(4), .Outer = held(-1), .Base = held(9), .Below = held(1)};
    const Form<Big> equal = {
        .Whole = held(3), .Root = held(4), .Outer = held(-1), .Base = held(9), .Below = held(1)};
    const Form<Big> small = {
        .Whole = held(1), .Root = held(4), .Outer = held(-1), .Base = held(9), .Below = held(-1)};
    const Form<Big> same  = {.Whole = held(1),
                             .Inner = held(1),
                             .Root  = held(4),
                             .Outer = held(2),
                             .Base  = held(5),
                             .Bend  = held(2),
                             .Below = held(1)};
    const Form<Big> under = {
        .Whole = held(1), .Root = held(4), .Outer = held(1), .Base = held(-9), .Below = held(1)};
    CHECK(Corner::Signed(zero, Expansion::Of(0))  == 0);
    CHECK(Corner::Signed(outer, Expansion::Of(0)) == 1);
    CHECK(Corner::Signed(mixed, Expansion::Of(0)) == 1);
    CHECK(Corner::Signed(equal, Expansion::Of(0)) == 0);
    CHECK(Corner::Signed(small, Expansion::Of(0)) == 1);
    CHECK(Corner::Signed(same, Expansion::Of(0))  == 1);
    CHECK(Corner::Signed(under, Expansion::Of(0)) == 1);

    const double    max   = std::numeric_limits<double>::max();
    const Form<Big> third = {.Whole = held(-1), .Below = held(3)};
    CHECK(Corner::Pinpointed(third, -0.3333, 1e-3) == -1.0 / 3.0);
    CHECK(Corner::Pinpointed(third, std::numeric_limits<double>::quiet_NaN(), 0) == -1.0 / 3.0);
    CHECK(Corner::Pinpointed(third, 7, 1e-3)        == -1.0 / 3.0);
    CHECK(Corner::Pinpointed(third, -1e308, 1e308)  == -1.0 / 3.0);
    CHECK(Corner::Pinpointed(third, 1e308, 1e308)   == -1.0 / 3.0);
    CHECK(Corner::Pinpointed(third, 1.7e308, 4e307) == -1.0 / 3.0);
    CHECK(Corner::Pinpointed(third, -1, 1e-3)       == -1.0 / 3.0);
    CHECK(Corner::Pinpointed(zero, 1, 1)            == 0.0);
    const Form<Big> odd = {
        .Whole = Big{.Value = Expansion::Sum(Expansion::Of(1), Expansion::Of(0x3p-53))},
        .Below = held(1)};
    CHECK(Corner::Pinpointed(odd, 1, 1e-15) == 1 + 0x1p-51);
    const Form<Big> largest = {.Whole = held(max), .Below = held(1)};
    CHECK(Corner::Pinpointed(largest, max, 0) == max);
  }

  TEST_CASE("fillet: parts that would run past a turn, and places on a curve or off it")
  {
    const Curve bowed = Polyline{
        .Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 2}, {.Position = {.X = 10, .Y = 0}}}};
    const std::expected<Corner::Guide, Misfit> bow = Corner::Guided(bowed, {.X = 5, .Y = -5});
    REQUIRE(bow);
    CHECK(Corner::Cut(*bow, {.Along = 1.5}, -1).error() == Misfit::TooLarge);
    CHECK(Corner::Cut(*bow, {.Along = -0.5}, 1).error() == Misfit::TooLarge);
    const Curve arc =
        Arc{.Center = {.X = 0, .Y = 0}, .Radius = 1, .Start = Degrees(0), .Sweep = Degrees(30)};
    const std::expected<Corner::Guide, Misfit> bent = Corner::Guided(arc, {.X = 1, .Y = 0});
    REQUIRE(bent);
    CHECK(Corner::Cut(*bent, {.Along = -11.5}, 1).error() == Misfit::TooLarge);
    CHECK(Corner::Cut(*bent, {.Along = 12.5}, -1).error() == Misfit::TooLarge);
    CHECK(Corner::Within(*bent, 12.0 - 1e-15));
    CHECK(!Corner::Within(*bent, 6.0));
    CHECK(!Corner::Within(*bow, 1.2));

    const Curve                                one   = Bezier::Cubic{.Start = {.X = 0, .Y = 0},
                                                                     .Out   = {.X = 1, .Y = 0},
                                                                     .In    = {.X = 2, .Y = 0},
                                                                     .End   = {.X = 3, .Y = 0}};
    const Curve                                other = Bezier::Cubic{.Start = {.X = 0, .Y = 1},
                                                                     .Out   = {.X = 1, .Y = 1},
                                                                     .In    = {.X = 2, .Y = 1},
                                                                     .End   = {.X = 3, .Y = 1}};
    const std::expected<Corner::Guide, Misfit> near  = Corner::Guided(one, {});
    const std::expected<Corner::Guide, Misfit> far   = Corner::Guided(other, {});
    CHECK(!Corner::Polished({.Guides = {&*near, &*far}, .Distances = {0.1, 0.1}}, {0.5, 0.5}));
  }

  TEST_CASE("fillet: corners picked right on them, circles that do not fit, and touching circles")
  {
    const Line across = {.Start = {.X = 0, .Y = 0}, .End = {.X = 10, .Y = 0}};
    const Line upward = {.Start = {.X = 12, .Y = -2}, .End = {.X = 12, .Y = 10}};
    CHECK(Why(Fillet(across, upward,
                     {.Radius = 0, .Picks = {Point{.X = 12, .Y = 0}, Point{.X = 12, .Y = 5}}}))
          == Misfit::TooLarge);
    const std::expected<Filleted, Misfit> picked = Fillet(
        across, upward, {.Radius = 1, .Picks = {Point{.X = 12, .Y = 0}, Point{.X = 12, .Y = 5}}});
    REQUIRE(picked);
    CHECK(Distance(picked->Touches[0], Point{.X = 12, .Y = 0}) == 1);

    const Arc  ring = Circle({.X = 0, .Y = 0}, 3);
    const Line high = {.Start = {.X = -5, .Y = 2.9}, .End = {.X = 5, .Y = 2.9}};
    CHECK(Why(Fillet(high, ring,
                     {.Radius = 0.5, .Picks = {Point{.X = 0, .Y = 2.9}, Point{.X = 0.2, .Y = 3}}}))
          == Misfit::TooLarge);

    const std::expected<Filleted, Misfit> kissing =
        Fillet(Circle({.X = 0, .Y = 0}, 1), Circle({.X = 3, .Y = 0}, 2),
               {.Radius = 0, .Picks = {Point{.X = 0, .Y = 1}, Point{.X = 3, .Y = 2}}});
    REQUIRE(kissing);
    CHECK(kissing->Touches[0] == Point{.X = 1, .Y = 0});

    const Line low = {.Start = {.X = 1, .Y = 3}, .End = {.X = 8, .Y = 3}};
    const std::expected<Filleted, Misfit> span =
        Fillet(across, low, {.Radius = 1, .Picks = {Point{.X = 1, .Y = 0}, Point{.X = 5, .Y = 3}}});
    REQUIRE(span);
    CHECK(std::get<Line>(span->Second) == Line{.Start = {.X = 0, .Y = 3}, .End = {.X = 8, .Y = 3}});
  }

  TEST_CASE("fillet: bridges kept either way round an arc, and smooth curves at their edges")
  {
    const Arc left  = Circle({.X = 0, .Y = 0}, 2);
    const Arc right = {
        .Center = {.X = 7, .Y = 0}, .Radius = 2, .Start = Degrees(90), .Sweep = Degrees(180)};
    CHECK(Fillet(left, right, {.Radius = 2, .Picks = {Point{.X = 0, .Y = 2}, At(right, 0.6)}})
              .has_value());
    CHECK(Fillet(right, left, {.Radius = 2, .Picks = {At(right, 0.1), Point{.X = 0, .Y = 2}}})
              .has_value());
    const Line far = {.Start = {.X = -5, .Y = -4}, .End = {.X = 5, .Y = -4}};
    CHECK(
        Fillet(far, left, {.Radius = 1, .Picks = {Point{.X = 4, .Y = -4}, Point{.X = -2, .Y = 0}}})
            .has_value());

    const Ellipse back = {.Map   = {.X = {.X = 4, .Y = 0}, .Y = {.X = 0, .Y = 2}},
                          .Start = Degrees(0),
                          .Sweep = Degrees(-360)};
    const std::expected<Filleted, Misfit> reversed =
        Fillet(back, Line{.Start = {.X = -6, .Y = 1}, .End = {.X = 6, .Y = 1}},
               {.Radius = 0.3, .Picks = {Point{.X = 0, .Y = 2}, Point{.X = 5, .Y = 1}}});
    CHECK(reversed.has_value());
    const Arc ring = Circle({.X = 0, .Y = 5}, 1);
    CHECK(Fillet(back, ring, {.Radius = 1, .Picks = {Point{.X = 0, .Y = 2}, Point{.X = 1, .Y = 5}}})
              .has_value());

    const Bezier::Cubic cusp  = {.Start = {.X = 0, .Y = 0},
                                 .Out   = {.X = 0, .Y = 0},
                                 .In    = {.X = 4, .Y = 4},
                                 .End   = {.X = 4, .Y = 4}};
    const Bezier::Cubic other = {.Start = {.X = 0, .Y = 4},
                                 .Out   = {.X = 0, .Y = 4},
                                 .In    = {.X = 4, .Y = 0},
                                 .End   = {.X = 4, .Y = 0}};
    CHECK(Fillet(cusp, Line{.Start = {.X = 0, .Y = 2}, .End = {.X = 4, .Y = 2}},
                 {.Radius = 0.3, .Picks = {At(cusp, 0.9), Point{.X = 0.5, .Y = 2}}})
              .has_value());
    CHECK(
        Fillet(cusp, other, {.Radius = 0.3, .Picks = {At(cusp, 0.9), At(other, 0.9)}}).has_value());
    CHECK(Fillet(cusp, Ray{.Start = {.X = 2, .Y = -1}, .Way = {.X = 0, .Y = 1}},
                 {.Radius = 0, .Picks = {At(cusp, 0.9), Point{.X = 2, .Y = 0}}})
              .has_value());
    CHECK(Fillet(cusp, Ray{.Start = {.X = -1, .Y = 1}, .Way = {.X = 1, .Y = 0}},
                 {.Radius = 0, .Picks = {At(cusp, 0.9), Point{.X = 0, .Y = 1}}})
              .has_value());

    const Polyline dome = {.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                        {.Position = {.X = 10, .Y = 0}, .Bulge = 1},
                                        {.Position = {.X = 0, .Y = 0}}}};
    const std::expected<Reshaped, Misfit> wide = Fillet(dome, 3, Exacting);
    REQUIRE(wide);
    REQUIRE(wide->Left.size() == 1);
    CHECK(wide->Left[0].Why   == Misfit::TooLarge);
    CHECK(Fillet(
              Polyline{.Vertices = {{.Position = {.X = std::numeric_limits<double>::quiet_NaN()}}}},
              1, Exacting)
              .error()
          == Misfit::NotFinite);
  }

  TEST_CASE("fillet: a circle first, a Bézier curve kept toward its start, smooth joints passed by")
  {
    const Arc  above  = Circle({.X = 0, .Y = 2}, 1);
    const Line across = {.Start = {.X = -5, .Y = 0}, .End = {.X = 5, .Y = 0}};
    const std::expected<Filleted, Misfit> first = Fillet(
        above, across, {.Radius = 0.5, .Picks = {Point{.X = -1, .Y = 2}, Point{.X = -4, .Y = 0}}});
    REQUIRE(first);
    CHECK(first->Round->Center == Point{.X = 0, .Y = 0.5});

    const Bezier::Cubic                   arch = {.Start = {.X = 0, .Y = 0},
                                                  .Out   = {.X = 2, .Y = 4},
                                                  .In    = {.X = 6, .Y = 4},
                                                  .End   = {.X = 8, .Y = 0}};
    const std::expected<Filleted, Misfit> start =
        Fillet(arch, Line{.Start = {.X = -5, .Y = 1}, .End = {.X = 15, .Y = 1}},
               {.Radius = 0.1, .Picks = {At(arch, 0.02), Point{.X = 5, .Y = 1}}});
    REQUIRE(start);
    CHECK(std::get<Bezier::Cubic>(start->First).Start == arch.Start);

    const Bezier::Cubic low  = {.Start = {.X = 0, .Y = 0},
                                .Out   = {.X = 1, .Y = 0},
                                .In    = {.X = 2, .Y = 0},
                                .End   = {.X = 3, .Y = 0}};
    const Bezier::Cubic high = {.Start = {.X = 0, .Y = 1},
                                .Out   = {.X = 1, .Y = 1},
                                .In    = {.X = 2, .Y = 1},
                                .End   = {.X = 3, .Y = 1}};
    CHECK(Why(Fillet(low, high,
                     {.Radius = 0.1, .Picks = {Point{.X = 1, .Y = 0}, Point{.X = 1, .Y = 1}}}))
          == Misfit::Apart);

    const Ellipse oval = {.Map   = {.X = {.X = 4, .Y = 0}, .Y = {.X = 0, .Y = 2}},
                          .Start = Degrees(0),
                          .Sweep = Degrees(360)};
    CHECK(Fillet(Circle({.X = 0, .Y = 5}, 1), oval,
                 {.Radius = 1, .Picks = {Point{.X = 1, .Y = 5}, Point{.X = 0, .Y = 2}}})
              .has_value());
    CHECK(Fillet(oval, Circle({.X = 0, .Y = 5}, 2),
                 {.Radius = 0.75, .Picks = {Point{.X = -3, .Y = 1.3}, Point{.X = -2, .Y = 5}}})
              .has_value());

    const Polyline joint = {.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                         {.Position = {.X = 10, .Y = 0}, .Bulge = std::tan(Pi / 8)},
                                         {.Position = {.X = 15, .Y = 5}},
                                         {.Position = {.X = 15, .Y = 10}}}};
    const std::expected<Reshaped, Misfit> smooth = Fillet(joint, 1, Exacting);
    REQUIRE(smooth);
    CHECK(smooth->Left.empty());
    CHECK(smooth->Shape == joint);
    const Polyline dot = {.Vertices = {{.Position = {.X = 1, .Y = 1}}}, .Closed = true};
    CHECK(Fillet(dot, 1, Exacting)->Shape.Vertices.empty());
  }
} // namespace Geoxy
