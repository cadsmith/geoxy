/***************************************************************************************************
 *  @file      tests/geoxy/polish.cpp
 *  @brief     Tests for crossings polished in double-double: Bézier curves, splines, polylines,
 *             arcs and ellipses, each crossing the double nearest the exact point.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "crossing.hpp"
#include "crossings/bows.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/arc.hpp"
#include "geoxy/at.hpp"
#include "geoxy/bezier/cubic.hpp"
#include "geoxy/bezier/quadratic.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/spline.hpp"
#include "geoxy/xline.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy::Crossings
{
  namespace
  {
    /**
     *  Two curves crossing, and where, as an 80-digit oracle finds them and rounds them once.
     */
    struct Pair
    {
      Curve              First  = {}; /**< One curve.                                 */
      Curve              Second = {}; /**< The other.                                 */
      std::vector<Point> Points = {}; /**< Where they cross, each the double nearest. */
    };

    /**
     *  A double from its bits: a number of a case, written so where it lies near a named constant
     *  that a literal would be taken for.
     *
     *  @param  bits The bits.
     *  @return      The double.
     */
    constexpr auto Bits(std::uint64_t bits) -> double
    {
      return std::bit_cast<double>(bits);
    }

    /**
     *  Two curves met: every crossing exactly an expected point, every expected point found; and
     *  polished again from places 1e-9 off, the same point.
     *
     *  @param pair The curves and the points.
     */
    void Exact(const Pair &pair)
    {
      const Intersections meeting =
          Intersect(pair.First, pair.Second, {.Length = 1e-9, .Radians = 1e-9});
      std::vector<const Hit *> crossings;
      for (const Hit &hit : meeting.Hits)
      {
        if (hit.Kind == Contact::Cross)
        {
          crossings.push_back(&hit);
        }
      }
      CHECK(crossings.size() == pair.Points.size());
      for (const Hit *hit : crossings)
      {
        CHECK(std::ranges::find(pair.Points, hit->At) != pair.Points.end());
        const std::optional<Crossed> again =
            Polished({&pair.First, &pair.Second},
                     {hit->First.Along + 1e-9, hit->Second.Along - 1e-9}, hit->At);
        CHECK(again.value_or(Crossed{}).At == hit->At);
      }
    }

    const Bezier::Cubic Wave = {.Start = {.X = 0, .Y = 0},
                                .Out   = {.X = 1, .Y = 3},
                                .In    = {.X = 2, .Y = -3},
                                .End   = {.X = 3, .Y = 0}}; /**< A cubic crossed by the others. */

    const Bezier::Quadratic Hump = {.Start   = {.X = 0, .Y = 0},
                                    .Control = {.X = 1, .Y = 2},
                                    .End     = {.X = 2, .Y = 0}}; /**< Its top at (1, 1). */

    const Ellipse Leaning = {
        .Map = {.X = {.X = 2, .Y = 0.5}, .Y = {.X = -0.25, .Y = 1}, .Origin = {.X = 1.5, .Y = 0.1}},
        .Start = {},
        .Sweep = Turn}; /**< An ellipse sheared, its numbers exact. */

    /**
     *  Where the hump crosses height 1/2: 1 +- sqrt(1/2).
     *
     *  @return The points.
     */
    auto HumpAtHalf() -> std::vector<Point>
    {
      return {Point{.X = 0.2928932188134525, .Y = 0.5}, Point{.X = 1.7071067811865475, .Y = 0.5}};
    }
  } // namespace

  TEST_CASE("polish: Bézier curves, splines and straight pieces, to the double nearest")
  {
    Exact({.First  = Hump,
           .Second = Line{.Start = {.X = -1, .Y = 0.5}, .End = {.X = 3, .Y = 0.5}},
           .Points = HumpAtHalf()});
    Exact({.First  = Ray{.Start = {.X = -1, .Y = 0.5}, .Way = {.X = 1, .Y = 0}},
           .Second = Hump,
           .Points = HumpAtHalf()});
    Exact({.First  = Xline{.Through = {.X = 0, .Y = 0.5}, .Way = {.X = 1, .Y = 0}},
           .Second = Hump,
           .Points = HumpAtHalf()});
    Exact({.First  = Wave,
           .Second = Bezier::Cubic{.Start = {.X = 0, .Y = 1},
                                   .Out   = {.X = 1, .Y = -2},
                                   .In    = {.X = 2, .Y = 4},
                                   .End   = {.X = 3, .Y = -1}},
           .Points = {Point{.X = 0.2077068641349944, .Y = 0.4996681141146634},
                      Point{.X = 1.1965469021948743, .Y = 0.43655090732253543},
                      Point{.X = 2.8589041284069734, .Y = -0.36543610847612573}}});
    Exact({.First  = Spline{.Degree  = 2,
                            .Knots   = {0, 0, 0, 1, 1, 1},
                            .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 1}, {.X = 2, .Y = 0}},
                            .Weights = {1, 2, 1}},
           .Second = Line{.Start = {.X = -1, .Y = 0.25}, .End = {.X = 3, .Y = 0.25}},
           .Points = {Point{.X = 0.260490027112548, .Y = 0.25},
                      Point{.X = 1.739509972887452, .Y = 0.25}}});
    Exact({.First  = Polyline{.Vertices = {{.Position = {.X = -1, .Y = -1}, .Bulge = 0},
                                           {.Position = {.X = 1, .Y = 1.5}, .Bulge = 0},
                                           {.Position = {.X = 3, .Y = -0.5}, .Bulge = 0}},
                              .Closed   = false},
           .Second = Wave,
           .Points = {Point{.X = 0.22606562566852254, .Y = 0.5325820320856531},
                      Point{.X = 0.4317345322407574, .Y = 0.7896681653009469},
                      Point{.X = 2.8609461421185713, .Y = -0.3609461421185715}}});
  }

  TEST_CASE("polish: arcs and ellipses by their quadratics, to the double nearest")
  {
    const Arc ring = {.Center = {.X = 1.5, .Y = 0}, .Radius = 1.25, .Start = {}, .Sweep = Turn};
    const std::vector<Point> onRing  = {Point{.X = 0.5962112430534608, .Y = 0.8635194744862617},
                                        Point{.X = 2.4037887569465393, .Y = -0.8635194744862617}};
    const std::vector<Point> onSlant = {Point{.X = 0.4667355008944791, .Y = 0.8144634912703572},
                                        Point{.X = 0.5429219059327847, .Y =  0.8511624226808954},
                                        Point{.X = 2.1364986526622194, .Y = -0.7828379526076579},
                                        Point{.X = 2.8152357140700413, .Y = -0.45608430297837743}};
    Exact({.First = Wave, .Second = ring, .Points = onRing});
    Exact({.First = ring, .Second = Wave, .Points = onRing});
    Exact({.First = Wave, .Second = Leaning, .Points = onSlant});
    Exact({.First  = Leaning,
           .Second = Ellipse{.Map   = {.X      = {.X = 0.5, .Y = -0.125},
                                       .Y      = {.X = 0.75, .Y = 1.5},
                                       .Origin = {.X = 1.25, .Y = 0.3}},
                             .Start = {},
                             .Sweep = Turn},
           .Points = {Point{.X = 2.124692311258562, .Y = 1.2177074652334803},
                      Point{.X =  1.1276583371808278, .Y =  1.0674094620870362},
                      Point{.X = 0.35464684774436944, .Y = -0.9789798135082255},
                      Point{.X = 0.9035682926140716, .Y = -1.017188992227753}}});
  }

  TEST_CASE("polish: bulges by their exact circles, against arcs, ellipses, Béziers and bulges")
  {
    /**
     *  Cases of the curves audit, each crossing the double nearest its exact rational reference.
     */
    Exact(
        {.First = Arc{.Center = {.X = -0x1.f213579b34e08p-2, .Y = 0x1.0342ee172dd9ap-2},
                      .Radius = 0x1.3708d0b75d4e0p-1,
                      .Start  = Degrees(90.0),
                      .Sweep  = Degrees(-142.0)},
         .Second =
             Polyline{
                 .Vertices = {{.Position = {.X = 0x1.681e829bf57d4p-2, .Y = 0x1.8d434348f0338p-2},
                               .Bulge    = 0x1.6d488db8bc6bcp-1},
                              {.Position = {.X = -0x1.46dae7586f19ep-1, .Y = -0x1.c10c912c1cf54p-2},
                               .Bulge    = 0x1.e76aaa61ef998p-3}},
                 .Closed   = false},
         .Points = {Point{.X = 0x1.38e21442ed14fp-4, .Y = 0x1.ed7789179f702p-2}}});
    Exact(
        {.First = Arc{.Center = {.X = 0x1.01eca3a83f924p-2, .Y = 0x1.fd0b34cad85a0p-3},
                      .Radius = 0x1.a3c4c7c08ef46p-1,
                      .Start  = Degrees(30.0),
                      .Sweep  = Degrees(84.0)},
         .Second =
             Polyline{
                 .Vertices = {{.Position = {.X = -0x1.298fe5c7e51c6p-1, .Y = 0x1.ad9a2b6a9be2ep-1},
                               .Bulge    = -0x1.8500e56368550p-1},
                              {.Position = {.X = 0x1.95f340e794ad4p-1, .Y = 0x1.27333895a9cfcp-2},
                               .Bulge    = 0x1.7098a08ff367ap+0},
                              {.Position = {.X = -0x1.420b1fc3696f0p-1, .Y = 0x1.fd33c9e747b18p-2},
                               .Bulge    = 0x0.0p+0},
                              {.Position = {.X = 0x1.96548d2cde29ep-1, .Y = 0x1.38c97b187eee0p-5},
                               .Bulge    = 0x1.a1fd6825e6e00p-8}},
                 .Closed   = false},
         .Points = {Point{.X = 0x1.c1b1d601b1de6p-1, .Y = 0x1.8e11625b8989ep-1},
                    Point{.X = 0x1.7e82bd78eee3ep-2, .Y = 0x1.0f30ddd4b86acp+0}}});
    Exact(
        {.First =
             Ellipse{.Map   = {.X      = {.X = 0x1.53ca56a19472dp-1, .Y = 0x1.7d23f2b622e15p-3},
                               .Y      = {.X = -0x1.04d08456a91c1p-2, .Y = 0x1.d1098fb42ccc1p-1},
                               .Origin = {.X = -0x1.ae9fbaae2e8b4p-2, .Y = 0x1.f9ef5ec7e8cecp-2}},
                     .Start = Degrees(165.0),
                     .Sweep = Degrees(360.0)},
         .Second =
             Polyline{
                 .Vertices = {{.Position = {.X = -0x1.399b71f785d5cp-2, .Y = -0x1.749c5dbc27500p-3},
                               .Bulge    = -0x1.e3709baca3ee0p-3},
                              {.Position = {.X = -0x1.b87bbf5865012p-1, .Y = -0x1.82acb73d44bb8p-3},
                               .Bulge    = 0x0.0p+0},
                              {.Position = {.X = 0x1.37176b5e3af50p-3, .Y = -0x1.b7dc5358b9e80p-7},
                               .Bulge    = 0x0.0p+0},
                              {.Position = {.X = 0x1.92ed717b4e0e8p-2, .Y = 0x1.842506afb2df0p-3},
                               .Bulge    = 0x1.52c77996a0a40p+0}},
                 .Closed   = true},
         .Points = {Point{.X = -0x1.a1f957b688c72p-1, .Y = -0x1.731304078fbf8p-3},
                    Point{.X = -0x1.8d33faead07f1p-1, .Y = -0x1.c75f9d60a9072p-3},
                    Point{.X =  0x1.07e2bb133be4fp-2, .Y =  0x1.3511e061f210ap-4},
                    Point{.X = 0x1.28710f6fc2b88p-2, .Y = 0x1.8be7da001c0dfp-2}}});
    Exact(
        {.First =
             Ellipse{.Map   = {.X      = {.X = -0x1.8bac873432c6bp-2, .Y = 0x1.6b04f74364824p-1},
                               .Y      = {.X = -0x1.91bd475bbb5e1p-1, .Y = -0x1.b5e07ab7ab0efp-2},
                               .Origin = {.X = -0x1.6e8396bc12816p-2, .Y = 0x1.5afa6d550bc48p-3}},
                     .Start = Degrees(57.0),
                     .Sweep = Degrees(360.0)},
         .Second =
             Polyline{
                 .Vertices = {{.Position = {.X = -0x1.a024193cc034cp-2, .Y = -0x1.d1f6430c789b0p-2},
                               .Bulge    = 0x0.0p+0},
                              {.Position = {.X = 0x1.8be0a61e9af70p-1, .Y = -0x1.006032960e4a0p-5},
                               .Bulge    = 0x0.0p+0},
                              {.Position = {.X = 0x1.f49067452a894p-1, .Y = -0x1.12f2442eb7006p-1},
                               .Bulge    = 0x0.0p+0},
                              {.Position = {.X = -0x1.9550cae6bf622p-1, .Y = 0x1.04e43a6ce2670p-1},
                               .Bulge    = 0x0.0p+0}},
                 .Closed   = false},
         .Points = {Point{.X = 0x1.935270bb3a676p-2, .Y = -0x1.889148cc925cep-3},
                    Point{.X = 0x1.a52d75921cb1ep-2, .Y = -0x1.4a572fb19fd59p-3}}});
    Exact(
        {.First =
             Bezier::Quadratic{.Start   = {.X = 0x1.fca80b0458218p-2, .Y = -0x1.0b52fa020f8a0p-2},
                               .Control = {.X = -0x1.71f8bc7a293a2p-1, .Y = 0x1.2c7b66b66c0e4p-1},
                               .End     = {.X = -0x1.8d256d4582f54p-1, .Y = -0x1.6c8117ba902ecp-1}},
         .Second =
             Polyline{
                 .Vertices = {{.Position = {.X = -0x1.75759d1876920p-3, .Y = -0x1.84662f84958e0p-2},
                               .Bulge    = 0x1.1cadb9bff2ad0p-2},
                              {.Position = {.X = -0x1.ac1adf4ea037cp-1, .Y = 0x1.340ff51178048p-2},
                               .Bulge    = 0x0.0p+0}},
                 .Closed   = false},
         .Points = {Point{.X = -0x1.aba489917741fp-2, .Y = 0x1.bfd8b182e5b77p-5}}});
    Exact(
        {.First =
             Bezier::Quadratic{.Start   = {.X = -0x1.979dca0e2db60p-1, .Y = -0x1.1c4536aca73d0p-1},
                               .Control = {.X = -0x1.f4addcd2186f8p-3, .Y = 0x1.bed14176d2a8ap-1},
                               .End     = {.X = 0x1.16a1d77e7d4c4p-1, .Y = 0x1.b543a29ed97acp-1}},
         .Second =
             Polyline{
                 .Vertices = {{.Position = {.X = -0x1.d676b3b3c4398p-1, .Y = -0x1.c5344c0ba64e6p-1},
                               .Bulge    = 0x0.0p+0},
                              {.Position = {.X = -0x1.a720e15f4af08p-1, .Y = -0x1.507599ce93898p-3},
                               .Bulge    = 0x0.0p+0},
                              {.Position = {.X = 0x1.c07e810367ad4p-1, .Y = -0x1.740c820040ce4p-2},
                               .Bulge    = 0x0.0p+0},
                              {.Position = {.X = 0x1.94c9c09fe012cp-1, .Y = -0x1.faf66f910927cp-1},
                               .Bulge    = 0x0.0p+0}},
                 .Closed   = false},
         .Points = {Point{.X = -0x1.46c59e6f026e2p-1, .Y = -0x1.7d85ecff57f15p-3}}});
    Exact(
        {.First = Bezier::Cubic{.Start = {.X = 0x1.ba3ad3f37eae0p-3, .Y = 0x1.29f64338f96b0p-1},
                                .Out   = {.X = -0x1.16589dfc53632p-1, .Y = -0x1.2d6e4e9717560p-2},
                                .In    = {.X = 0x1.2a78cfe54397ep-1, .Y = 0x1.8968b6f6b4d88p-1},
                                .End   = {.X = -0x1.ae5f7c3db2e30p-1, .Y = -0x1.d73dd0c694d54p-2}},
         .Second =
             Polyline{.Vertices =
                          {{.Position = {.X = -0x1.c67ba0fddffc0p-5, .Y = 0x1.cc61ccbbe47c0p-1},
                            .Bulge    = 0x0.0p+0},
                           {.Position = {.X = -0x1.3150d101a7990p-1, .Y = -0x1.45218b8ca4ad4p-1},
                            .Bulge    = 0x0.0p+0},
                           {.Position = {.X = 0x1.b5d67f37a3b18p-3, .Y = -0x1.0fab472c91620p-4},
                            .Bulge    = 0x0.0p+0},
                           {.Position = {.X = 0x1.0e458a8fe46f8p-2, .Y = Bits(0x3fd4500c874478e0)},
                            .Bulge    = -0x1.3aba6f028fa6ap+0}},
                      .Closed = false},
         .Points = {Point{.X = -0x1.9c7f6a376c47cp-2, .Y = -0x1.60dcbde366d85p-4}}});
    Exact(
        {.First = Bezier::Cubic{.Start = {.X = 0x1.07c53f3d73f98p-3, .Y = -0x1.3b1b1401bd4c0p-1},
                                .Out   = {.X = 0x1.093005edbcc98p-3, .Y = 0x1.69dbfed3ba3aep-1},
                                .In    = {.X = -0x1.be285069e5598p-1, .Y = 0x1.4143c5f49ed62p-1},
                                .End   = {.X = -0x1.48e210c86d7a0p-4, .Y = -0x1.944edb3ba67aep-1}},
         .Second =
             Polyline{
                 .Vertices = {{.Position = {.X = -0x1.e5e2f09771ca8p-3, .Y = 0x1.6f846c7ec8832p-1},
                               .Bulge    = 0x1.626e1cc1f0770p-2},
                              {.Position = {.X = -0x1.cdb48e2af5880p-7, .Y = 0x1.1cbd21c945460p-4},
                               .Bulge    = -0x1.6ef7e1b658c05p+0},
                              {.Position = {.X = -0x1.3425b99fe3200p-9, .Y = 0x1.2a2241908e880p-5},
                               .Bulge    = 0x0.0p+0}},
                 .Closed   = false},
         .Points = {Point{.X = 0x1.013c57319e5fap-6, .Y = 0x1.94a4ed383e8e4p-5},
                    Point{.X = 0x1.7817bab8bc637p-8, .Y = 0x1.2756ec125bc3bp-4},
                    Point{.X = -0x1.cf6ca2de03680p-3, .Y = 0x1.4c27eaa0dc726p-2}}});
    Exact(
        {.First =
             Polyline{
                 .Vertices = {{.Position = {.X = 0x1.08b5970f92e1ap-1, .Y = -0x1.d6a8ca12e6590p-4},
                               .Bulge    = 0x1.7347bd4f14f10p+0},
                              {.Position = {.X = -0x1.adf9f09494980p-7, .Y = 0x1.8982a4787d370p-4},
                               .Bulge    = 0x0.0p+0},
                              {.Position = {.X = -0x1.9dddba22e5028p-1, .Y = 0x1.d034cf939f6a8p-3},
                               .Bulge    = 0x1.3ea926814df78p+0},
                              {.Position = {.X = 0x1.0b0e208426f44p-2, .Y = 0x1.ec59d57a335c2p-1},
                               .Bulge    = -0x1.15b075b318652p-1},
                              {.Position = {.X = -0x1.44108401801e8p-3, .Y = -0x1.5aad234895440p-4},
                               .Bulge    = 0x0.0p+0}},
                 .Closed   = false},
         .Second =
             Polyline{
                 .Vertices = {{.Position = {.X = -0x1.e466fac804db0p-2, .Y = 0x1.9c75ec5ce790cp-2},
                               .Bulge    = 0x0.0p+0},
                              {.Position = {.X = -0x1.fcdce189a3e30p-2, .Y = 0x1.2fbd73c209c5cp-2},
                               .Bulge    = 0x0.0p+0},
                              {.Position = {.X = 0x1.9ea5c9be4afb4p-2, .Y = 0x1.7d13635cfac64p-2},
                               .Bulge    = -0x1.003ec6afc728fp-1},
                              {.Position = {.X = 0x1.2643a6e4d30e0p-3, .Y = -0x1.ab718413a6c26p-1},
                               .Bulge    = 0x0.0p+0}},
                 .Closed   = false},
         .Points = {Point{.X = 0x1.2a44f68638ba8p-1, .Y = -0x1.2ebe0caa04b54p-9},
                    Point{.X = 0x1.06b9f19795319p-3, .Y = 0x1.655a6c2303bc8p-2},
                    Point{.X = 0x1.ce275222ec9e6p-2, .Y = 0x1.40a56dd9a0ba8p-2},
                    Point{.X = 0x1.66e00b7f66ff5p-2, .Y = 0x1.7867c1716def3p-2}}});
    Exact(
        {.First =
             Polyline{
                 .Vertices = {{.Position = {.X = -0x1.e9def9cc07b6ap-1, .Y = 0x1.ba365cc36ad14p-2},
                               .Bulge    = 0x1.185f5efa19ba0p-5},
                              {.Position = {.X = 0x1.f9edab7ab0fe8p-1, .Y = 0x1.77afe6fdbe620p-1},
                               .Bulge    = 0x1.a006bff5b0428p-2},
                              {.Position = {.X = 0x1.c573fd94cc3c0p-1, .Y = 0x1.5f0e2c7ab4b8ap-1},
                               .Bulge    = 0x0.0p+0},
                              {.Position = {.X = -0x1.994ef99fd5f1cp-1, .Y = -0x1.28511c6cf1facp-1},
                               .Bulge    = 0x1.a4459b574ef1cp-1},
                              {.Position = {.X = 0x1.bb0d4aa5d5e48p-1, .Y = 0x1.f3df37d7fffa6p-1},
                               .Bulge    = -0x1.7cc04a068b566p+0},
                              {.Position = {.X = 0x1.ec7c42248d594p-1, .Y = -0x1.8c1c8e9ee8b2cp-1},
                               .Bulge    = 0x0.0p+0}},
                 .Closed   = false},
         .Second =
             Polyline{
                 .Vertices = {{.Position = {.X = 0x1.1497bf8546f24p-2, .Y = -0x1.a371595014b90p-1},
                               .Bulge    = 0x0.0p+0},
                              {.Position = {.X = 0x1.f5e94ef78b460p-2, .Y = 0x1.bd9d9cb2fee20p-4},
                               .Bulge    = -0x1.30109e5753da0p-4},
                              {.Position = {.X = -0x1.6cb8eef96a910p-1, .Y = 0x1.d502a2623ca90p-3},
                               .Bulge    = 0x0.0p+0}},
                 .Closed   = false},
         .Points = {Point{.X = 0x1.d900d368dfd7fp-4, .Y = 0x1.b9604cd0244f9p-4},
                    Point{.X = 0x1.2c95c419dde69p-2, .Y = -0x1.70d9417a6b72bp-1}}});
  }

  TEST_CASE(
      "polish: a crossing's place on a conic is its point's, whatever start it is polished from")
  {
    const Curve ring =
        Arc{.Center = {.X = 1.5, .Y = 0}, .Radius = 1.25, .Start = {}, .Sweep = Turn};
    const Curve leaning = Leaning;
    const Curve wave    = Wave;
    const Curve bulge   = Polyline{.Vertices = {{.Position = {.X = -1, .Y = 0.5}, .Bulge = 0.5},
                                                {.Position = {.X = 4, .Y = 0.5}, .Bulge = 0}},
                                   .Closed   = false};
    for (const std::array<const Curve *, 2> &curves :
         {std::array<const Curve *, 2>{&ring, &leaning},
          std::array<const Curve *, 2>{&leaning, &ring}, std::array<const Curve *, 2>{&wave, &ring},
          std::array<const Curve *, 2>{&leaning, &wave},
          std::array<const Curve *, 2>{&bulge, &wave}, std::array<const Curve *, 2>{&wave, &bulge}})
    {
      const Intersections meeting =
          Intersect(*curves[0], *curves[1], {.Length = 1e-9, .Radians = 1e-9});
      CHECK(!meeting.Hits.empty());
      for (const Hit &hit : meeting.Hits)
      {
        const std::array<double, 2>  places = {hit.First.Along, hit.Second.Along};
        const std::optional<Crossed> near   = Polished(curves, places, hit.At);
        const std::optional<Crossed> off =
            Polished(curves, {places[0] + 1e-7, places[1] - 1e-7}, hit.At);
        REQUIRE(near.has_value());
        REQUIRE(off.has_value());
        CHECK(off->At == near->At);
        for (std::size_t side = 0; side < 2; ++side)
        {
          /**
           *  A conic's place is its point's, as Located places it; a polynomial's is where Newton's
           *  method settles, the same from either start.
           */
          const std::optional<Place> placed = Located(*curves.at(side), places.at(side), near->At);
          CHECK(off->Along.at(side) == doctest::Approx(near->Along.at(side)).epsilon(1e-15));
          if (!std::holds_alternative<Bezier::Cubic>(*curves.at(side)))
          {
            REQUIRE(placed.has_value());
            CHECK(std::bit_cast<std::uint64_t>(near->Along.at(side))
                  == std::bit_cast<std::uint64_t>(placed->Along));
            CHECK(std::bit_cast<std::uint64_t>(off->Along.at(side))
                  == std::bit_cast<std::uint64_t>(placed->Along));
          }
        }
      }
    }
  }

  TEST_CASE("polish: a crossing's points on its curves are At's at its places")
  {
    /**
     *  On a polyline the point Located worked out already, on an arc At's own: both At's at the
     *  crossing's places, to the bit, two bulges crossed and a bulge and an arc.
     */
    const Curve bulged =
        Polyline{.Vertices = {{.Position = {.X = -0x1.298fe5c7e51c6p-1, .Y = 0x1.ad9a2b6a9be2ep-1},
                               .Bulge    = -0x1.8500e56368550p-1},
                              {.Position = {.X = 0x1.95f340e794ad4p-1, .Y = 0x1.27333895a9cfcp-2},
                               .Bulge    = 0x1.7098a08ff367ap+0},
                              {.Position = {.X = -0x1.420b1fc3696f0p-1, .Y = 0x1.fd33c9e747b18p-2},
                               .Bulge    = 0x0.0p+0}},
                 .Closed   = false};
    const Curve across   = Polyline{.Vertices = {{.Position = {.X = -1, .Y = 0.3}, .Bulge = 0.4},
                                                 {.Position = {.X = 1, .Y = 0.8}, .Bulge = 0}},
                                    .Closed   = false};
    const Curve round    = Arc{.Center = {.X = 0x1.01eca3a83f924p-2, .Y = 0x1.fd0b34cad85a0p-3},
                               .Radius = 0x1.a3c4c7c08ef46p-1,
                               .Start  = Degrees(30.0),
                               .Sweep  = Degrees(84.0)};
    std::size_t polished = 0;
    for (const Curve *other : {&across, &round})
    {
      const Intersections meeting = Intersect(bulged, *other, {.Length = 1e-9, .Radians = 1e-9});
      for (const Hit &hit : meeting.Hits)
      {
        const std::optional<Crossed> crossed =
            Polished({&bulged, other}, {hit.First.Along, hit.Second.Along}, hit.At);
        polished            += crossed.has_value() ? 1U : 0U;
        const Crossed found  = crossed.value_or(Crossed{});
        CHECK(found.On[0] == At(bulged, found.Along[0]));
        CHECK(found.On[1] == At(*other, found.Along[1]));
      }
    }
    CHECK(polished >= 2);
  }

  TEST_CASE("polish: nothing for a tangency, a far start or a crossing past an end")
  {
    const Curve hump  = Hump;
    const Curve top   = Line{.Start = {.X = -1, .Y = 1}, .End = {.X = 3, .Y = 1}};
    const Curve under = Line{.Start = {.X = -1, .Y = 0}, .End = {.X = 3, .Y = 0}};
    const Curve unit  = Arc{.Center = {}, .Radius = 1, .Start = {}, .Sweep = Turn};
    const Curve next  = Arc{.Center = {.X = 2, .Y = 0}, .Radius = 1, .Start = {}, .Sweep = Turn};
    const Curve bulge = Polyline{.Vertices = {{.Position = {.X = -1, .Y = 0}, .Bulge = 1},
                                              {.Position = {.X = 1, .Y = 0}, .Bulge = 0}},
                                 .Closed   = false};
    const Curve none  = Polyline{};
    const Curve flat  = Ellipse{.Map = Scaling(0, 0), .Start = {}, .Sweep = Turn};

    CHECK(!Polished({&bulge, &hump}, {0.5, 0.5}, {}));
    CHECK(!Polished({&hump, &bulge}, {0.5, 0.5}, {}));
    CHECK(!Polished({&bulge, &unit}, {0.5, 0.25}, {}));
    CHECK(!Polished({&unit, &bulge}, {0.25, 0.5}, {}));
    CHECK(!Polished({&none, &hump}, {0.0, 0.5}, {}));
    CHECK(!Polished({&flat, &hump}, {0.5, 0.5}, {}));

    /**
     *  Parallel at the start, or tangent so that Newton's method creeps: nothing either way.
     */
    CHECK(!Polished({&top, &under}, {0.5, 0.5}, {}));
    CHECK(!Polished({&top, &unit}, {0.25, 0.25}, {.X = 0, .Y = 1}));
    CHECK(!Polished({&unit, &next}, {0.0, 0.5}, {.X = 1, .Y = 0}));
    CHECK(!Polished({&hump, &top}, {0.501, 0.5}, {.X = 1, .Y = 1}));
    CHECK(!Polished({&top, &unit}, {0.251, 0.25}, {.X = 0.004, .Y = 1}));
    CHECK(!Polished({&unit, &next}, {0.0, 0.5}, {.X = 1, .Y = 1e-3}));

    /**
     *  Started far from the crossing it settles at: nothing, the crossing is not this hit's.
     */
    const Curve half = Line{.Start = {.X = -1, .Y = 0.5}, .End = {.X = 3, .Y = 0.5}};
    CHECK(!Polished({&hump, &half}, {0.2, 0.3}, {}));
    CHECK(!Polished({&half, &unit}, {0.4, 0.3}, {}));
    const Curve over = Arc{.Center = {.X = 1, .Y = 0}, .Radius = 1, .Start = {}, .Sweep = Turn};
    CHECK(!Polished({&unit, &over}, {0.0, 0.5}, {.X = 0.6, .Y = 0.9}));
    CHECK(!Polished({&unit, &over}, {0.0, 0.5}, {.X = 0.5, .Y = 0.9}));
    CHECK(!Polished({&hump, &half}, {0.1464466094067262, 0.4}, {}));

    /**
     *  Where the crossing lies past a line's end, on either curve.
     */
    const Curve shy = Line{.Start = {.X = -1, .Y = 0.5}, .End = {.X = 0.25, .Y = 0.5}};
    CHECK(!Polished({&hump, &shy}, {0.1464466, 1.0343146}, {}));
    CHECK(!Polished({&shy, &hump}, {1.0343146, 0.1464466}, {}));
    const Curve late = Line{.Start = {.X = 0.5, .Y = 0.5}, .End = {.X = 3, .Y = 0.5}};
    CHECK(!Polished({&hump, &late}, {0.1464466, -0.0828427}, {}));

    /**
     *  A polyline's straight segment crossing toward a bulge beyond its vertex: Newton's method
     *  steps onto the bulge, which it does not work out as a polynomial, and gives up.
     */
    const Curve bend = Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0},
                                             {.Position = {.X = 1, .Y = 0}, .Bulge = 0.5},
                                             {.Position = {.X = 2, .Y = 0}, .Bulge = 0}},
                                .Closed   = false};
    const Curve post = Line{.Start = {.X = 1.5, .Y = -1}, .End = {.X = 1.5, .Y = 1}};
    CHECK(!Polished({&bend, &post}, {0.99, 0.5}, {.X = 1.5, .Y = -0.1}));
    CHECK(!Polished({&post, &bend}, {0.5, 0.99}, {.X = 1.5, .Y = -0.1}));

    /**
     *  A place not a number never settles.
     */
    const double nan   = std::numeric_limits<double>::quiet_NaN();
    const Curve  slant = Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 1, .Y = 1}};
    CHECK(!Polished({&half, &slant}, {nan, 0.5}, {}));
    CHECK(!Polished({&bend, &post}, {nan, 0.5}, {}));
  }

  TEST_CASE("polish: two bulges Newton's method leaves, crossed by their exact circles")
  {
    /**
     *  A bulge with half a circle of radius 1 and one of radius 1.5, both below their chords, their
     *  circles crossing at (1/3, +-2 sqrt(2) / 3): from a hit 0.5 off, farther than Newton's method
     *  may move one, the crossing of their exact circles nearer it, the double nearest, either
     *  first.
     */
    const Polyline small = {.Vertices = {{.Position = {.X = -1, .Y = 0}, .Bulge = 1},
                                         {.Position = {.X = 1, .Y = 0}, .Bulge = 0}},
                            .Closed   = false};
    const Polyline large = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 1},
                                         {.Position = {.X = 3, .Y = 0}, .Bulge = 0}},
                            .Closed   = false};
    const Curve    one   = small;
    const Curve    other = large;
    const Point    under = {.X = 0x1.5555555555555p-2, .Y = -0x1.e2b7dddfefa66p-1};
    const Point    off   = {.X = 0.1, .Y = -0.5};
    CHECK(Bowed({&small, &large}, {0.5, 0.5}, off) == under);
    CHECK(Bowed({&large, &small}, {0.5, 0.5}, off) == under);
    for (const auto &[first, second] : {std::pair{&one, &other}, std::pair{&other, &one}})
    {
      const std::optional<Crossed> crossed = Polished({first, second}, {0.5, 0.5}, off);
      REQUIRE(crossed);
      CHECK(crossed->At    == under);
      CHECK(crossed->On[0] == At(*first, crossed->Along[0]));
      CHECK(crossed->On[1] == At(*second, crossed->Along[1]));
    }

    /**
     *  Hit midway between the two crossings, neither nearer by more than rounding: the exact way's
     *  choice, one of the two.
     */
    const std::optional<Point> midway = Bowed({&small, &large}, {0.5, 0.5}, {.X = 1.0 / 3, .Y = 0});
    REQUIRE(midway);
    CHECK(midway->X           == under.X);
    CHECK(std::abs(midway->Y) == -under.Y);

    /**
     *  Half circles touching at (1, 0), their radicand nothing: no bound settles it, and the exact
     *  way gives the one point, its y the +0 Newton's method's sums leave. So does a crossing on
     *  the x axis, the other far: no bound settles a y of nothing.
     */
    const Polyline right = {.Vertices = {{.Position = {.X = 0, .Y = -1}, .Bulge = 1},
                                         {.Position = {.X = 0, .Y = 1}, .Bulge = 0}},
                            .Closed   = false};
    const Polyline left  = {.Vertices = {{.Position = {.X = 2, .Y = 1}, .Bulge = 1},
                                         {.Position = {.X = 2, .Y = -1}, .Bulge = 0}},
                            .Closed   = false};
    const Polyline cup   = {.Vertices = {{.Position = {.X = 2, .Y = 1}, .Bulge = -1},
                                         {.Position = {.X = 0, .Y = 1}, .Bulge = 0}},
                            .Closed   = false};
    for (const Polyline *across : {&left, &cup})
    {
      const std::optional<Point> touch = Bowed({&right, across}, {0.5, 0.5}, {.X = 0.9, .Y = 0.1});
      REQUIRE(touch);
      CHECK(*touch == Point{.X = 1, .Y = 0});
      CHECK_FALSE(std::signbit(touch->Y));
    }
    const Curve                  bowl    = right;
    const Curve                  hollow  = left;
    const std::optional<Crossed> touched = Polished({&bowl, &hollow}, {0.5, 0.5}, {.X = 1, .Y = 0});
    REQUIRE(touched);
    CHECK(touched->At == Point{.X = 1, .Y = 0});

    /**
     *  Circles apart, or of one center: nothing, the exact way's too.
     */
    const Polyline far   = {.Vertices = {{.Position = {.X = 4, .Y = 0}, .Bulge = 1},
                                         {.Position = {.X = 6, .Y = 0}, .Bulge = 0}},
                            .Closed   = false};
    const Polyline round = {.Vertices = {{.Position = {.X = -2, .Y = 0}, .Bulge = 1},
                                         {.Position = {.X = 2, .Y = 0}, .Bulge = 0}},
                            .Closed   = false};
    CHECK_FALSE(Bowed({&small, &far}, {0.5, 0.5}, {.X = 2.5, .Y = 0}));
    CHECK_FALSE(Bowed({&small, &round}, {0.5, 0.5}, {.X = 0, .Y = -1.5}));
    const Curve apart = far;
    CHECK_FALSE(Polished({&one, &apart}, {0.5, 0.5}, {.X = 2.5, .Y = 0}));
  }

  TEST_CASE("polish: two bulges by one step from a guess at their crossing, settled or left")
  {
    const Belly unit  = {.Start = {.X = -1, .Y = 0}, .End = {.X = 1, .Y = 0}, .Bulge = 1};
    const Belly broad = {.Start = {.X = 0, .Y = 0}, .End = {.X = 3, .Y = 0}, .Bulge = 1};
    const Point under = {.X = 0x1.5555555555555p-2, .Y = -0x1.e2b7dddfefa66p-1};
    const Point above = {.X = 0x1.5555555555555p-2, .Y = 0x1.e2b7dddfefa66p-1};

    /**
     *  From a point near either crossing, the guess and its one step settle it, the double nearest;
     *  from one midway, the guess is not taken; from one far off, no other crossing is told farther
     *  than this one; a crossing on an axis no bound settles.
     */
    CHECK(Hopped({unit, broad}, {.X = 0.3333, .Y = -0.9428}) == under);
    CHECK(Hopped({unit, broad}, {.X = 0.3333, .Y = 0.9428})  == above);
    CHECK_FALSE(Hopped({unit, broad}, {.X = 1.0 / 3, .Y = 0}));
    CHECK_FALSE(Hopped({unit, broad}, {.X = 100, .Y = -50}));
    const Belly right = {.Start = {.X = 0, .Y = -1}, .End = {.X = 0, .Y = 1}, .Bulge = 1};
    const Belly cup   = {.Start = {.X = 2, .Y = 1}, .End = {.X = 0, .Y = 1}, .Bulge = -1};
    CHECK_FALSE(Hopped({right, cup}, {.X = 0.9, .Y = 0.1}));

    /**
     *  A bulge nearly a whole circle of radius 2.5e5 across one of radius 70: the step's error, the
     *  larger bulge times the step squared, is not shown under its reach, and Bowed settles the
     *  crossing the other ways.
     */
    const Belly whole = {.Start = {.X = -0x1.12fa232b85db2p-1, .Y = 0x1.50dce1a98e708p-1},
                         .End   = {.X = 0x1.4c425c7c2abbcp-2, .Y = -0x1.b7e08504bc80dp-1},
                         .Bulge = 0x1.496f06be6633bp+19};
    const Belly flat  = {.Start = {.X = 0x1.7338a9e7bcbep-2, .Y = 0x1.3c9a86bc9898p-3},
                         .End   = {.X = 0x1.727b4ed6beb5p-2, .Y = -0x1.db05de3ecc996p-1},
                         .Bulge = 0x1.fd2489827376ap-9};
    const Point near  = {.X = -0x1.066d5248fd009p-1, .Y = 0x1.a78a34e8ecf2ep-1};
    CHECK_FALSE(Hopped({whole, flat}, near));
    const auto polyline = [](const Belly &belly) -> Polyline
    {
      return {.Vertices = {{.Position = belly.Start, .Bulge = belly.Bulge},
                           {.Position = belly.End, .Bulge = 0}},
              .Closed   = false};
    };
    const Polyline round = polyline(whole);
    const Polyline line  = polyline(flat);
    CHECK(Bowed({&round, &line}, {0.5, 0.5}, near));

    /**
     *  Both bulges moved 1e8 along x: their crossing's y under 2^-8 of its x, which Newton's method
     *  may leave short of its last bit, is crossed exactly first, the double nearest (1e8 + 1/3, -2
     *  sqrt(2) / 3).
     */
    const Curve far   = Polyline{.Vertices = {{.Position = {.X = 1e8 - 1, .Y = 0}, .Bulge = 1},
                                              {.Position = {.X = 1e8 + 1, .Y = 0}, .Bulge = 0}},
                                 .Closed   = false};
    const Curve wider = Polyline{.Vertices = {{.Position = {.X = 1e8, .Y = 0}, .Bulge = 1},
                                              {.Position = {.X = 1e8 + 3, .Y = 0}, .Bulge = 0}},
                                 .Closed   = false};
    const std::optional<Crossed> moved =
        Polished({&far, &wider}, {0.5, 0.5}, {.X = 1e8 + 0.3333, .Y = -0.9428});
    REQUIRE(moved);
    CHECK(moved->At == Point{.X = 0x1.7d78401555555p+26, .Y = -0x1.e2b7dddfefa66p-1});
  }
} // namespace Geoxy::Crossings
