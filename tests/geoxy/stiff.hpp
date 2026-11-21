/***************************************************************************************************
 *  @file      tests/geoxy/stiff.hpp
 *  @brief     Splines whose weights lie far apart, and what a point and a box ask of them, answered
 *             in exact rational arithmetic by tools/reference/stiff.py.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TESTS_STIFF_HPP
#define GEOXY_TESTS_STIFF_HPP

#include "geoxy/box.hpp"
#include "geoxy/point.hpp"
#include "geoxy/spline.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace Geoxy::Tests
{
  /**
   *  A spline whose weights lie a power of two apart: heavy ones for a power above 0, light ones
   *  for one below.
   *
   *  @param  shape Which: "conic", a quadratic of one span, its middle weight the power; "bow", a
   *                cubic of one span, its two middle weights the power; "paired", a cubic of two
   *                spans, one weight the power and one a power as far below 1; any other, a
   *                quadratic of three spans, every second weight the power.
   *  @param  power The power its weights lie apart by.
   *  @return       The spline.
   */
  inline auto Stiff(std::string_view shape, int power) -> Spline
  {
    const double weight = std::ldexp(1.0, power);
    if (shape == "conic")
    {
      return {.Degree  = 2,
              .Knots   = {0, 0, 0, 1, 1, 1},
              .Poles   = {{.X = 0, .Y = 0}, {.X = 1, .Y = 1}, {.X = 2, .Y = 0}},
              .Weights = {1, weight, 1}};
    }
    if (shape == "bow")
    {
      return {
          .Degree = 3,
          .Knots  = {0, 0, 0, 0, 1, 1, 1, 1},
          .Poles = {{.X = 0, .Y = 0}, {.X = 0.5, .Y = 1.5}, {.X = 1.5, .Y = 1.5}, {.X = 2, .Y = 0}},
          .Weights = {1, weight, weight, 1}};
    }
    if (shape == "paired")
    {
      return {.Degree  = 3,
              .Knots   = {0, 0, 0, 0, 1, 2, 2, 2, 2},
              .Poles   = {{.X = 0, .Y = 0},
                          {.X = 0.5, .Y =  1.5},
                          {.X = 1.5, .Y =  1.5},
                          {.X = 2.5, .Y = 0.25},
                          {.X = 3, .Y = 1}},
              .Weights = {1, weight, 1, 1 / weight, 1}};
    }
    return {.Degree  = 2,
            .Knots   = {0, 0, 0, 1, 2, 3, 3, 3},
            .Poles   = {{.X = 0, .Y = 0},
                        {.X = 1, .Y = 1},
                        {.X = 2, .Y = 0},
                        {.X = 3, .Y = 1},
                        {.X = 4, .Y = 0}},
            .Weights = {1, weight, 1, weight, 1}};
  }

  /**
   *  A place of a spline as exact arithmetic tells it: its point, its parameter, each the double
   *  nearest, and how fast the point moves by the parameter there.
   */
  struct Told
  {
    Point  At    = {};  /**< The point.                          */
    double Along = 0.0; /**< The parameter.                      */
    double Speed = 0.0; /**< The point's speed by the parameter. */
  };

  /**
   *  What a point asks of a stiff spline, and the answer: the places, in order.
   */
  struct Asked
  {
    Spline            Curve  = {}; /**< The spline.           */
    Point             From   = {}; /**< The point asked from. */
    std::vector<Told> Places = {}; /**< The places, in order. */
  };

  /**
   *  A stiff spline cut to a box, and the answer: the stretches within, in order, each as the place
   *  it starts at and the place it ends at, neither with a point.
   */
  struct Boxed
  {
    Spline            Curve  = {}; /**< The spline.                         */
    Box               Within = {}; /**< The box.                            */
    std::vector<Told> Ends   = {}; /**< Each stretch's start, then its end. */
  };

  /**
   *  A step of a double at the splines' coordinates.
   */
  inline constexpr double Step = 0x1p-52;

  /**
   *  Whether a point is within a few steps of the one exact arithmetic gives.
   *
   *  @param  point The point.
   *  @param  told  The place.
   *  @return       Whether it is within sixteen steps.
   */
  inline auto Lands(Point point, const Told &told) -> bool
  {
    return std::max(std::abs(point.X - told.At.X), std::abs(point.Y - told.At.Y)) <= 16.0 * Step;
  }

  /**
   *  Whether a parameter tells a place as well as a double can: within two doubles of the one exact
   *  arithmetic gives, or near enough that the curve, at its speed there, moves less than eight
   *  steps. Beside a light weight a spline hardly moves, and its parameter need not be told to its
   *  last digit to tell the place.
   *
   *  @param  along The parameter.
   *  @param  told  The place.
   *  @return       Whether it does.
   */
  inline auto Tells(double along, const Told &told) -> bool
  {
    const double most = std::max(std::abs(along), std::abs(told.Along));
    const double gap  = std::abs(along - told.Along);
    return gap <= 2.0 * (most - std::nextafter(most, 0.0)) || gap * told.Speed <= 8.0 * Step;
  }

  /**
   *  The next number of a line, written as C writes a double in sixteens.
   *
   *  @param  words The line.
   *  @return       The number.
   */
  inline auto Number(std::istringstream &words) -> double
  {
    std::string word;
    words >> word;
    return std::strtod(word.c_str(), nullptr);
  }

  /**
   *  What points ask of stiff splines, read from a sheet: a line a question, the spline's shape and
   *  power, the point, then each place of the answer as its point, its parameter and its speed.
   *
   *  @param  sheet The sheet.
   *  @return       The questions, in the sheet's order.
   */
  inline auto Asks(std::string_view sheet) -> std::vector<Asked>
  {
    std::vector<Asked> asks;
    std::istringstream lines{std::string{sheet}};
    for (std::string line; std::getline(lines, line);)
    {
      std::istringstream words{line};
      std::string        shape;
      int                power = 0;
      if (!(words >> shape >> power))
      {
        continue;
      }
      Asked asked  = {.Curve = Stiff(shape, power), .From = {}, .Places = {}};
      asked.From.X = Number(words);
      asked.From.Y = Number(words);
      while (words.peek() != std::istringstream::traits_type::eof())
      {
        Told told  = {};
        told.At.X  = Number(words);
        told.At.Y  = Number(words);
        told.Along = Number(words);
        told.Speed = Number(words);
        asked.Places.push_back(told);
      }
      asks.push_back(asked);
    }
    return asks;
  }

  /**
   *  Stiff splines cut to boxes, read from a sheet: a line a box, the spline's shape and power, the
   *  box's least and greatest corners, then each stretch of the answer as the parameters it starts
   *  and ends at and the speeds there.
   *
   *  @param  sheet The sheet.
   *  @return       The boxes, in the sheet's order.
   */
  inline auto Boxes(std::string_view sheet) -> std::vector<Boxed>
  {
    std::vector<Boxed> boxes;
    std::istringstream lines{std::string{sheet}};
    for (std::string line; std::getline(lines, line);)
    {
      std::istringstream words{line};
      std::string        shape;
      int                power = 0;
      if (!(words >> shape >> power))
      {
        continue;
      }
      Boxed boxed        = {.Curve = Stiff(shape, power), .Within = {}, .Ends = {}};
      boxed.Within.Min.X = Number(words);
      boxed.Within.Min.Y = Number(words);
      boxed.Within.Max.X = Number(words);
      boxed.Within.Max.Y = Number(words);
      while (words.peek() != std::istringstream::traits_type::eof())
      {
        Told start  = {};
        Told end    = {};
        start.Along = Number(words);
        start.Speed = Number(words);
        end.Along   = Number(words);
        end.Speed   = Number(words);
        boxed.Ends.push_back(start);
        boxed.Ends.push_back(end);
      }
      boxes.push_back(boxed);
    }
    return boxes;
  }

  /**
   *  The place nearest a point on a stiff spline, each answered in exact rational arithmetic and
   *  moved less than two steps of a double by a step of the question.
   */
  inline constexpr std::string_view Nearests = R"(
conic 3 -0x1.bb6ccb3e71d6cp-1 0x1.0fca8dad01e7ep+1 0x1.3453da650066dp-1 0x1.30cb731218f36p-1 0x1.5b942ab73da79p-4 0x1.19addd6cfb451p+2
conic 3 0x1.baaf774ef11fcp+0 0x1.2b9f415b2c3e4p-2 0x1.b7ca5b666437ap+0 0x1.1ff5200c803ebp-2 0x1.f3c5c251a2938p-1 0x1.895d57b974187p+3
conic 8 0x1.985dc8e46fb40p-3 -0x1.3f2130d259140p-5 0x1.48959a113b003p-4 0x1.48958bbcf95f1p-4 0x1.652e7c1bd6cb5p-13 0x1.3263155b4e88dp+9
conic 8 0x1.10f2cdc1e1de0p+0 -0x1.fea6c7fdb2800p-10 0x1.88b8d6ad354bcp+0 0x1.dd1bd52d88117p-2 0x1.ff210b639112fp-1 0x1.9e79106334cc0p+7
conic 60 -0x1.23e936390e9ccp-1 0x1.2e797d29c9262p+0 0x1.3909c41a83af8p-2 0x1.3909c41a83af8p-2 0x1.c2de9c20a654ep-63 0x1.5d0abe0d6828ap+60
conic 60 0x1.5817168fcea57p+1 -0x1.058cb6740da46p-1 0x1.0000000000000p+1 0.0 0x1.0000000000000p+0 0x1.6a09e667f3bcdp+61
conic -60 0x1.802b989824150p-3 0x1.3fa25be0b3292p-1 0x1.802b989824150p-3 0x1.2a8913cf4c782p-61 0x1.f27ebaec7aad1p-3 0x1.d89ac991b6302p+0
conic -60 -0x1.de22b3f334950p-1 0x1.1775b91c42758p-2 0.0 0.0 0.0 0x1.6a09e667f3bcdp-59
conic 1000 0x1.eb51806f95680p-4 0x1.fbb9a023a95d8p-2 0x1.3b47001fc75bcp-2 0x1.3b47001fc75bcp-2 0x1.c787866dbf03fp-1003 0x1.5ad8c70b13376p+1000
conic 1000 0x1.dfd9d9fb1a06ep+0 0x1.ac1adc5fa5050p+0 0x1.19df7ecdba80fp+0 0x1.cc4102648afe2p-1 0x1.0000000000000p+0 0x1.d958cea2edf4dp+994
conic -1000 0x1.5a6cd73214e5fp+1 -0x1.215626306eb88p-1 0x1.0000000000000p+1 0.0 0x1.0000000000000p+0 0x1.6a09e667f3bcdp-999
conic -1000 0x1.803610e248086p+0 0x1.2ae6107b0283fp+1 0x1.803610e248086p+0 0x1.bb292edb1bf89p-1001 0x1.44b9c9b3148f4p-1 0x1.9d5e7b2025607p+1
bow 3 0x1.473f1b8a2d980p-5 0x1.1f13235083abep+0 0x1.746e37f3d3e78p-2 0x1.ed9e4f65ae118p-1 0x1.0c3aba1b9b4d7p-4 0x1.aa49ae97b2eafp+2
bow 3 0x1.56ac41fa101b8p+0 0x1.33c686f286518p+0 0x1.6b35490c5e697p+0 0x1.49a464cd71791p+0 0x1.a6e8ad35cb6e6p-1 0x1.0fcd971d83e70p+1
bow 8 0x1.3de676df87887p+1 0x1.14d07e6677fbep-1 0x1.e2dd5be4cf431p+0 0x1.5d5b2d254e143p-2 0x1.ffcdcb048a267p-1 0x1.6afee30d6719bp+9
bow 8 0x1.a7d32d76296a4p+0 0x1.8488d001b4f1cp-1 0x1.bcb608e041ac9p+0 0x1.9294f3e58c21ap-1 0x1.ff4483844d4e1p-1 0x1.14bb9dec8c485p+8
bow 60 0x1.41cdcc90bfa04p-1 -0x1.fe44d05652b02p-1 0.0 0.0 0.0 0x1.2f9422c23c47ep+62
bow 60 -0x1.0592688481f64p-1 -0x1.6a73e556bec5cp-1 0.0 0.0 0.0 0x1.2f9422c23c47ep+62
bow -60 0x1.b1eb2e4f0e708p-2 0x1.4e0859a267938p-1 0x1.b1eb2e4f0e708p-2 0x1.e213c14603ee4p-59 0x1.91a87ef3c378bp-2 0x1.0cfac8aff937ep+2
bow -60 0x1.ff76b047fdf60p-2 0x1.a85acf9eaa152p+0 0x1.ff76b047fdf60p-2 0x1.fb0f66b003495p-59 0x1.a32b955cce216p-2 0x1.2994be1f84653p+2
bow 1000 0x1.5b80df81bb312p+0 -0x1.bf89a474f1b3ap-2 0x1.0000000000000p+1 0.0 0x1.0000000000000p+0 0x1.2f9422c23c47ep+1002
bow 1000 -0x1.baea9a3b0e260p-4 0x1.c68fcfe1b0e56p+0 0x1.0000000000000p-1 0x1.8000000000000p+0 0x1.f576fe5fbc261p-503 0x1.25005a429bad9p+3
bow -1000 0x1.6dcc4561a3850p+0 0x1.2a9939a14d61ap+1 0x1.6dcc4561a3850p+0 0x1.0735eb5fd6c33p-998 0x1.26d3cc1858f33p-1 0x1.40bd6e320b190p+2
bow -1000 0x1.17bab31544faap+1 0x1.5aec5aa504a00p-9 0x1.0000000000000p+1 0.0 0x1.0000000000000p+0 0x1.2f9422c23c47ep-998
paired 3 -0x1.34eeee81ec373p-1 0x1.b78e96c922048p-3 0x1.09123b7d16110p-8 0x1.8d9706cc11716p-7 0x1.641b3676b7839p-12 0x1.2b01cad753799p+5
paired 3 0x1.0ee1acd0a1140p-6 0x1.0f0de1a4966f8p+0 0x1.4885013d8e51ep-2 0x1.e88e15a86edfcp-1 0x1.0c24bec8068a0p-4 0x1.8efcf44084fd8p+2
paired 8 0x1.436205d70bb98p+0 -0x1.5e1f63ad65ba0p-2 0x1.85102b86b620ep-6 0x1.23cc1bcb8b4dep-4 0x1.104814f38609ap-14 0x1.137b3beadc826p+10
paired 8 0x1.38c5aa8ecc58cp+0 -0x1.78ad117c7dc9cp-3 0x1.126e47b1a2e69p-4 0x1.9ba5564fd8c40p-3 0x1.a665d1a917eacp-13 0x1.c79e2131b2907p+9
paired 60 0x1.635c5d18cc74ep+1 0x1.11e0607a6bdc6p+0 0x1.62be503d4dec6p+1 0x1.0bb3dfe77a945p+0 0x1.fffff42fc035fp+0 0x1.b894531fb3364p+19
paired 60 0x1.4180c7fc1658ap+1 -0x1.d73db7ee02e68p-3 0x1.622fa594a3534p+1 0x1.0becf0f7bed20p+0 0x1.fffff41ab8df0p+0 0x1.bd044c9e172abp+19
paired -60 0x1.33ce671ade369p+1 -0x1.7bc5192d40910p-5 0x1.302bd10e1930dp+1 0x1.e6ac81b036cabp-3 0x1.101dc3fe8dc7bp-18 0x1.5568fa1b0807ap+16
paired -60 0x1.0dcaf47c7a703p+1 -0x1.e281e14fbda38p-1 0x1.fe5b4d923adc7p+0 0x1.98490adb76e7bp-3 0x1.40d38f8410cfbp-19 0x1.f16c95a37bb81p+18
paired 1000 0x1.ae68037a5ce5cp-1 -0x1.cc60cbeb7436cp-1 0.0 0.0 0.0 0x1.2f9422c23c47ep+1002
paired 1000 -0x1.1c1667d6b6420p-4 -0x1.b22ef23925cd0p-1 0.0 0.0 0.0 0x1.2f9422c23c47ep+1002
paired -1000 0x1.76fea120ed88ap+0 0x1.697118b1b956ep+0 0x1.971171d1e1a6bp+0 0x1.45a78e4181522p-3 0x1.848126170a008p-333 0x1.2655c221d9e4cp+333
paired -1000 0x1.ffb2bc1634e90p-4 0x1.08d26f8c7d2dcp+1 0x1.506b009330e1bp-2 0x1.0d2266dc271b0p-5 0x1.57bbcc8668a75p-334 0x1.48111bc707543p+333
waves 3 0x1.d293284200844p+1 0x1.10d29ec6625c4p-1 0x1.c72578e66e51dp+1 0x1.c57419334936fp-2 0x1.7a01820ad5651p+1 0x1.fb488d8865281p+2
waves 3 0x1.d38e09c8a6ad6p+1 -0x1.937fa69e0f706p-1 0x1.0000000000000p+2 0.0 0x1.8000000000000p+1 0x1.6a09e667f3bcdp+4
waves 8 0x1.152752a8bb214p+1 -0x1.896086de09e1bp-1 0x1.14dd03b37300cp+1 0x1.fa1f7c7b522fcp-1 0x1.8aa021fae2347p+0 0x1.efe28aeceec06p+1
waves 8 0x1.cb7dc8c482b4cp+1 0x1.909a631df5e60p-3 0x1.d93a0d5db5767p+1 0x1.362f73636f9aep-2 0x1.7fe437f07c063p+1 0x1.60988f4e91dc0p+8
waves 60 -0x1.bf1b6e495e0e2p-1 0x1.f2f260119a284p+0 0x1.1364a8eceb213p-1 0x1.1364a8eceb213p-1 0x1.29f72a39a639bp-61 0x1.35439a15bdeb4p+59
waves 60 0x1.f42c62681fb3ap+1 0x1.caef58877afd8p-3 0x1.ebbeb66fd401ep+1 0x1.44149902bfe1cp-3 0x1.8000000000000p+1 0x1.00864d1e46e2ap+61
waves -60 0x1.56fcd8ac5fd94p+0 0x1.22f691de567bap+1 0x1.56fcd8ac5fd94p+0 0x1.fffdaf649040fp-60 0x1.5629597cd51ebp-1 0x1.febd5a7b0c1e5p+1
waves -60 -0x1.e7bb81d88cad8p-3 -0x1.8407c9b79b7fap-1 0.0 0.0 0.0 0x1.6a09e667f3bcdp-59
waves 1000 0x1.9c21efb87d8c4p+0 0x1.292a2f5659b0ap+0 0x1.9c21efb87d8c4p+0 0x1.0000000000000p+0 0x1.6600d68afcac7p+0 0x1.c4c20adda1509p+1
waves 1000 -0x1.daf065d69de74p-1 -0x1.5676859bd3115p-1 0.0 0.0 0.0 0x1.6a09e667f3bcdp+1001
waves -1000 0x1.05c19dc9ef4d4p-1 -0x1.75233bd584e2ep-2 0x1.05c19dc9ef4d4p-1 0x1.7d49ef067b70fp-1000 0x1.d00cf844010b5p-2 0x1.8923bf25809c2p+1
waves -1000 0x1.627380b1364c0p+1 0x1.7e10c9f21a3a0p-4 0x1.627380b1364c0p+1 0x1.fdcefefd0a124p-1000 0x1.2de4fbbf6cdd3p+1 0x1.0770e2fb27b51p+2
)";

  /**
   *  The feet of the perpendiculars from a point to a stiff spline, each answered in exact rational
   *  arithmetic and moved less than two steps of a double by a step of the question.
   */
  inline constexpr std::string_view Perpendiculars = R"(
conic 3 -0x1.bb6ccb3e71d6cp-1 0x1.0fca8dad01e7ep+1 0x1.3453da650066dp-1 0x1.30cb731218f36p-1 0x1.5b942ab73da79p-4 0x1.19addd6cfb451p+2
conic 3 0x1.baaf774ef11fcp+0 0x1.2b9f415b2c3e4p-2 0x1.b7ca5b666437ap+0 0x1.1ff5200c803ebp-2 0x1.f3c5c251a2938p-1 0x1.895d57b974187p+3
conic 8 0x1.985dc8e46fb40p-3 -0x1.3f2130d259140p-5 0x1.48959a113b003p-4 0x1.48958bbcf95f1p-4 0x1.652e7c1bd6cb5p-13 0x1.3263155b4e88dp+9 0x1.013f06c92a65ep+0 0x1.fccfe92f50c7bp-1 0x1.7b492d10df091p-1 0x1.4f95c3bb90df2p-5 0x1.1e70b96676f42p+0 0x1.c318065d5368bp-1 0x1.f8b29e46b8309p-1 0x1.5150f8cd704ddp+3
conic 8 0x1.10f2cdc1e1de0p+0 -0x1.fea6c7fdb2800p-10 0x1.1074657765c24p-1 0x1.1073ca86fa335p-1 0x1.22854605930d4p-9 0x1.3e68cd281eec2p+7 0x1.ffddd3c3470bdp-1 0x1.fe00da6318d32p-1 0x1.eeddacaa187c7p-2 0x1.006c19f016724p-6 0x1.88b8d6ad354bcp+0 0x1.dd1bd52d88117p-2 0x1.ff210b639112fp-1 0x1.9e79106334cc0p+7
conic 60 -0x1.23e936390e9ccp-1 0x1.2e797d29c9262p+0 0x1.3909c41a83af8p-2 0x1.3909c41a83af8p-2 0x1.c2de9c20a654ep-63 0x1.5d0abe0d6828ap+60
conic 60 0x1.3d035f8e124acp+0 0x1.5fe7553b75678p+0 0x1.0000000000000p+0 0x1.0000000000000p+0 0x1.5bef4e6226f6bp-1 0x1.c37e198a85c13p-58
conic -60 -0x1.5e6c5d609c5a0p-4 0x1.fa4f852110de0p-4 0x1.166d383d21084p-121 0x1.c787e5510a86bp-122 0x1.c787e5510a86bp-63 0x1.c1d51a77c100bp-59
conic -60 0x1.62e3c73ad1884p+1 0x1.7b3f6c9f6598ap+0 0x1.0000000000000p+1 0x1.d5c638f0f8a38p-121 0x1.0000000000000p+0 0x1.14d0a4707d926p-58
conic 1000 0x1.eb51806f95680p-4 0x1.fbb9a023a95d8p-2 0x1.3b47001fc75bcp-2 0x1.3b47001fc75bcp-2 0x1.c787866dbf03fp-1003 0x1.5ad8c70b13376p+1000
conic 1000 0x1.dfd9d9fb1a06ep+0 0x1.ac1adc5fa5050p+0 0x1.19df7ecdba80fp+0 0x1.cc4102648afe2p-1 0x1.0000000000000p+0 0x1.d958cea2edf4dp+994
conic -1000 0x1.6cd6a07a5f20fp+1 0x1.f11685d3bb078p-1 0x1.0000000000000p+1 0.0 0x1.0000000000000p+0 0x1.848e7b22a15e6p-999
conic -1000 0x1.803610e248086p+0 0x1.2ae6107b0283fp+1 0x1.803610e248086p+0 0x1.bb292edb1bf89p-1001 0x1.44b9c9b3148f4p-1 0x1.9d5e7b2025607p+1
bow 3 0x1.473f1b8a2d980p-5 0x1.1f13235083abep+0 0x1.746e37f3d3e78p-2 0x1.ed9e4f65ae118p-1 0x1.0c3aba1b9b4d7p-4 0x1.aa49ae97b2eafp+2
bow 3 0x1.56ac41fa101b8p+0 0x1.33c686f286518p+0 0x1.6b35490c5e697p+0 0x1.49a464cd71791p+0 0x1.a6e8ad35cb6e6p-1 0x1.0fcd971d83e70p+1
bow 8 0x1.3de676df87887p+1 0x1.14d07e6677fbep-1 0x1.e2dd5be4cf431p+0 0x1.5d5b2d254e143p-2 0x1.ffcdcb048a267p-1 0x1.6afee30d6719bp+9
bow 8 0x1.a7d32d76296a4p+0 0x1.8488d001b4f1cp-1 0x1.9e6df2c50cd9bp-2 0x1.33a0e824ec6bep+0 0x1.542b23daea5dfp-8 0x1.8ad4623ba5a69p+5 0x1.f476474fe56ffp-2 0x1.69191cc5626bdp+0 0x1.433a5662c3936p-6 0x1.452963f2f0a02p+2 0x1.bcb608e041ac9p+0 0x1.9294f3e58c21ap-1 0x1.ff4483844d4e1p-1 0x1.14bb9dec8c485p+8
bow 60 0x1.41cdcc90bfa04p-1 -0x1.fe44d05652b02p-1 0x1.000000184929ep-1 0x1.7fffffff5a495p+0 0x1.8b7a3adc374acp-29 0x1.04d080ffe7e30p+0 0x1.41cdcc90bfa03p-1 0x1.8000000000000p+0 0x1.07373242fe80ep-3 0x1.0000000000000p+0 0x1.7ffffffc17a4fp+0 0x1.7ffffffe38e8ep+0 0x1.fffffff6ffe52p-1 0x1.32d9f89e53867p+0
bow 60 -0x1.0592688481f64p-1 -0x1.6a73e556bec5cp-1 0x1.7ffffffe9ad9cp+0 0x1.7ffffffcc3f4ap+0 0x1.fffffffb0dac0p-1 0x1.f1350f4fa6a5ap+0 0x1.f616dd2cc9516p+0 0x1.dbb6879a40bf3p-4 0x1.0000000000000p+0 0x1.02635974dbb3ap+62
bow -60 0x1.b1eb2e4f0e708p-2 0x1.4e0859a267938p-1 0x1.b1eb2e4f0e708p-2 0x1.e213c14603ee4p-59 0x1.91a87ef3c378bp-2 0x1.0cfac8aff937ep+2
bow -60 0x1.ff76b047fdf60p-2 0x1.a85acf9eaa152p+0 0x1.ff76b047fdf60p-2 0x1.fb0f66b003495p-59 0x1.a32b955cce216p-2 0x1.2994be1f84653p+2
bow 1000 0x1.5b80df81bb312p+0 -0x1.bf89a474f1b3ap-2 0x1.2f5d3766fdd73p-8 0x1.c70bd31a7cc2dp-7 0x1.9843e3e62be67p-1009 0x1.29fbd07b41ff1p+1002 0x1.0000000000000p-1 0x1.8000000000000p+0 0x1.f6692b96e2fa0p-501 0x1.486ac3388d667p+0 0x1.5b80df81bb312p+0 0x1.8000000000000p+0 0x1.b701bf0376624p-1 0x1.0000000000000p+0 0x1.8000000000000p+0 0x1.8000000000000p+0 0x1.0000000000000p+0 0x1.0725f25e94ba6p+0
bow 1000 -0x1.baea9a3b0e260p-4 0x1.c68fcfe1b0e56p+0 0x1.0000000000000p-1 0x1.8000000000000p+0 0x1.f576fe5fbc261p-503 0x1.25005a429bad9p+3
bow -1000 0x1.6dcc4561a3850p+0 0x1.2a9939a14d61ap+1 0x1.6dcc4561a3850p+0 0x1.0735eb5fd6c33p-998 0x1.26d3cc1858f33p-1 0x1.40bd6e320b190p+2
bow -1000 0x1.4c169154293e0p+1 0x1.a1a2cde6379fcp+0 0x1.0000000000000p+1 0.0 0x1.0000000000000p+0 0x1.a49d46cfe0f0ap-997
paired 3 -0x1.34eeee81ec373p-1 0x1.b78e96c922048p-3 0x1.09123b7d16110p-8 0x1.8d9706cc11716p-7 0x1.641b3676b7839p-12 0x1.2b01cad753799p+5
paired 3 0x1.0ee1acd0a1140p-6 0x1.0f0de1a4966f8p+0 0x1.4885013d8e51ep-2 0x1.e88e15a86edfcp-1 0x1.0c24bec8068a0p-4 0x1.8efcf44084fd8p+2
paired 8 0x1.436205d70bb98p+0 -0x1.5e1f63ad65ba0p-2 0x1.85102b86b620ep-6 0x1.23cc1bcb8b4dep-4 0x1.104814f38609ap-14 0x1.137b3beadc826p+10 0x1.01334acfa88e9p-1 0x1.7fdd1f97dcbd3p+0 0x1.44cace976f59ap-1 0x1.41db469efffc9p-7 0x1.9804f66db7891p+0 0x1.4a3bc72e69a13p+0 0x1.c95db177c6ca5p+0 0x1.55357952d94f4p+3
paired 8 0x1.38c5aa8ecc58cp+0 -0x1.78ad117c7dc9cp-3 0x1.126e47b1a2e69p-4 0x1.9ba5564fd8c40p-3 0x1.a665d1a917eacp-13 0x1.c79e2131b2907p+9 0x1.012d2b91b086cp-1 0x1.7fdbd6442a57fp+0 0x1.422547a09d0c6p-1 0x1.40844ec480853p-7 0x1.85d336bfd5dd1p+0 0x1.4dea96c945615p+0 0x1.c7967ac602a7dp+0 0x1.46bb2569bb584p+3
paired 60 0x1.635c5d18cc74ep+1 0x1.11e0607a6bdc6p+0 0x1.62be503d4dec6p+1 0x1.0bb3dfe77a945p+0 0x1.fffff42fc035fp+0 0x1.b894531fb3364p+19
paired 60 0x1.4180c7fc1658ap+1 -0x1.d73db7ee02e68p-3 0x1.75085598bc7bep-3 0x1.17c640328d5cep-1 0x1.873294e593f97p-63 0x1.eabcf97c649e9p+60 0x1.0000000000000p-1 0x1.8000000000000p+0 0x1.c9609edb80db6p-2 0x1.882e9e90fb32bp-59 0x1.622fa594a3534p+1 0x1.0becf0f7bed20p+0 0x1.fffff41ab8df0p+0 0x1.bd044c9e172abp+19
paired -60 -0x1.37722a6dd6faep-2 0x1.54f24634e2584p+0 0x1.6be6472b90c46p-109 0x1.a7cee17fae255p-111 0x1.f138b8e7b8eeap-57 0x1.0ccd7962d2e0fp-51
paired -60 0x1.c4fb3b52b3baap+1 -0x1.6e4ed96cf649ep-1 0x1.4000000000000p+1 0x1.0000000000000p-2 0x1.cbc324cc614c6p+0 0x1.b89e61deafcc4p-58
paired 1000 0x1.7bac26a43998cp+1 0x1.901101733f7b0p-3 0x1.6bbf78bb11331p-2 0x1.10cf9a8c4ce65p+0 0x1.a2be7057d4e06p-1001 0x1.973e01cc889c8p+998 0x1.0000000000000p-1 0x1.8000000000000p+0 0x1.491ccdff1c22dp-2 0x1.541ac5f7b4bb1p-998 0x1.8000000000000p+1 0x1.0000000000000p+0 0x1.0000000000000p+1 0x1.608e626cfb535p-997
paired 1000 0x1.ae68037a5ce5cp-1 -0x1.cc60cbeb7436cp-1 0x1.0000000000000p-1 0x1.8000000000000p+0 0x1.8ba7498d604abp-1 0x1.9e1a786bd57c6p-999 0x1.49f6640c8b877p+0 0x1.579b8597174b5p+0 0x1.0000000000000p+1 0x1.034eeb6fa0e3dp+333
paired -1000 0x1.76fea120ed88ap+0 0x1.697118b1b956ep+0 0x1.971171d1e1a6bp+0 0x1.45a78e4181522p-3 0x1.848126170a008p-333 0x1.2655c221d9e4cp+333 0x1.4000000000000p+1 0x1.0000000000000p-2 0x1.c174b4b062051p+0 0x1.2edf47dfe092dp-998 0x1.5bdecfb7ecda5p+1 0x1.2738de4f8d1ddp-1 0x1.0000000000000p+1 0x1.068f46807f8e1p+999
paired -1000 0x1.ffb2bc1634e90p-4 0x1.08d26f8c7d2dcp+1 0x1.506b009330e1bp-2 0x1.0d2266dc271b0p-5 0x1.57bbcc8668a75p-334 0x1.48111bc707543p+333 0x1.4000000000000p+1 0x1.0000000000000p-2 0x1.dc0edd4387729p+0 0x1.dcdbbbff510d1p-997 0x1.4dea3a37355c4p+1 0x1.a6faba9680536p-2 0x1.0000000000000p+1 0x1.05cc4887c2bc9p+997
waves 3 0x1.d293284200844p+1 0x1.10d29ec6625c4p-1 0x1.c72578e66e51dp+1 0x1.c57419334936fp-2 0x1.7a01820ad5651p+1 0x1.fb488d8865281p+2
waves 3 0x1.d38e09c8a6ad6p+1 -0x1.937fa69e0f706p-1 0x1.6a4525bbb31e7p+1 0x1.b2cc3d223e87bp-1 0x1.e7669ae3630fcp+0 0x1.e6989f2a1ad38p-1 0x1.79b98ea494081p+1 0x1.decfad4433413p-1 0x1.1aaf9b0453294p+1 0x1.0ca21a10f38e7p-2
waves 8 0x1.152752a8bb214p+1 -0x1.896086de09e1bp-1 0x1.65a2704bf8197p-1 0x1.65a1a1a0a8ee6p-1 0x1.26e0b47109c03p-8 0x1.0a440d93391efp+6 0x1.ff3094f3d9bb6p-1 0x1.fe8fd3a59cb53p-1 0x1.b5faa4357c034p-2 0x1.ea337230f050ap-7 0x1.14dd03b37300cp+1 0x1.fa1f7c7b522fcp-1 0x1.8aa021fae2347p+0 0x1.efe28aeceec06p+1 0x1.80163afd327d1p+1 0x1.fed3f455d19aap-1 0x1.3eafc72e05df9p+1 0x1.60495350ae9a8p-7 0x1.bbbf65d789c58p+1 0x1.11021ac09678fp-1 0x1.7fb71803a3ff4p+1 0x1.3d4aca3840d1fp+7
waves 8 0x1.cb7dc8c482b4cp+1 0x1.909a631df5e60p-3 0x1.7f7fb08d2bde5p+1 0x1.fdff3cfcc010cp-1 0x1.ffa74cff8fc3bp+0 0x1.ab89dea72542cp-6 0x1.7fac68cd731b3p+1 0x1.fe9d947c83ec1p-1 0x1.0dc83a2039592p+1 0x1.a03ec5015e7a2p-7 0x1.d93a0d5db5767p+1 0x1.362f73636f9aep-2 0x1.7fe437f07c063p+1 0x1.60988f4e91dc0p+8
waves 60 -0x1.bf1b6e495e0e2p-1 0x1.f2f260119a284p+0 0x1.1364a8eceb213p-1 0x1.1364a8eceb213p-1 0x1.29f72a39a639bp-61 0x1.35439a15bdeb4p+59
waves 60 0x1.f42c62681fb3ap+1 0x1.caef58877afd8p-3 0x1.ebbeb66fd401ep+1 0x1.44149902bfe1cp-3 0x1.8000000000000p+1 0x1.00864d1e46e2ap+61
waves -60 0x1.9fbb25999d430p+1 -0x1.2836d703c2130p-1 0x1.9fbb25999d430p+1 0x1.bf01c65ed269ap-60 0x1.3d036a83eaf83p+1 0x1.e1959d421765fp+1
waves -60 0x1.f0c699209950cp+1 0x1.41814cc5389c8p-2 0x1.f0c699209950cp+1 0x1.74f0772386497p-61 0x1.5e6c1c6e8e964p+1 0x1.27fae888a3633p+0
waves 1000 0x1.9c21efb87d8c4p+0 0x1.292a2f5659b0ap+0 0x1.9c21efb87d8c4p+0 0x1.0000000000000p+0 0x1.6600d68afcac7p+0 0x1.c4c20adda1509p+1
waves 1000 -0x1.e805c40f75decp-1 0x1.c6dab12d8fc2cp+0 0x1.a5af9e4ba9a6cp-2 0x1.a5af9e4ba9a6cp-2 0x1.6674cf3732645p-1002 0x1.f506db5e4f1bbp+999
waves -1000 0x1.05c19dc9ef4d4p-1 -0x1.75233bd584e2ep-2 0x1.05c19dc9ef4d4p-1 0x1.7d49ef067b70fp-1000 0x1.d00cf844010b5p-2 0x1.8923bf25809c2p+1
waves -1000 0x1.627380b1364c0p+1 0x1.7e10c9f21a3a0p-4 0x1.627380b1364c0p+1 0x1.fdcefefd0a124p-1000 0x1.2de4fbbf6cdd3p+1 0x1.0770e2fb27b51p+2
)";

  /**
   *  The places of a stiff spline whose tangent passes through a point, each answered in exact
   *  rational arithmetic and moved less than two steps of a double by a step of the question.
   */
  inline constexpr std::string_view Touches = R"(
conic 3 -0x1.bb6ccb3e71d6cp-1 0x1.0fca8dad01e7ep+1 0x1.1ae99ca232379p+0 0x1.b37333910aac6p-1 0x1.6933318e641d0p-1 0x1.9621cfbbda43ap-1
conic 3 0x1.1726a2efb427ap+1 0x1.8b53ec111bf1cp+0 0x1.d73e8606770d8p-1 0x1.bb38b0a1b5683p-1 0x1.57f6b6508170bp-2 0x1.4da63a4b583c0p-1
conic 8 0x1.7eeb5943ee1edp+1 -0x1.d59f819ae28b4p-1 0x1.0387c73bfc61bp+0 0x1.f8ab504d04df3p-1 0x1.c250c9920d699p-1 0x1.7aec371ee02b5p-3
conic 8 -0x1.dafa2b444dc78p-1 -0x1.cf8a80a7e1330p-3 0x1.fe5b92851bd5dp-1 0x1.fd6b7e7d27ef2p-1 0x1.48019dd26857cp-2 0x1.c03a5767de7cbp-6
conic 60 -0x1.23e936390e9ccp-1 0x1.2e797d29c9262p+0 0x1.0000000000000p+0 0x1.0000000000000p+0 0x1.0ed96c0785384p-1 0x1.045242d410603p-58
conic 60 0x1.5817168fcea57p+1 -0x1.058cb6740da46p-1 0x1.0000000000000p+0 0x1.0000000000000p+0 0x1.9e6b5232bdaa2p-1 0x1.37f2cdfd4d131p-56
conic -60 -0x1.5e6c5d609c5a0p-4 0x1.fa4f852110de0p-4 0x1.0000000000000p+1 0x1.fbe69f3be533cp-117 0x1.0000000000000p+0 0x1.0e6c96ff416ebp-55
conic -60 0x1.62e3c73ad1884p+1 0x1.7b3f6c9f6598ap+0 0x1.4057eed3f8571p-120 0x1.be3b3ea6fda69p-121 0x1.be3b3ea6fda69p-62 0x1.0f9c25dde77f5p-58
conic 1000 0x1.eb51806f95680p-4 0x1.fbb9a023a95d8p-2 0x1.0000000000000p+0 0x1.0000000000000p+0 0x1.5ecd2dd93d173p-2 0x1.8f9161a5efe32p-998
conic 1000 0x1.dfd9d9fb1a06ep+0 0x1.ac1adc5fa5050p+0 0x1.0000000000000p+0 0x1.0000000000000p+0 0x1.0fe43be880b39p-2 0x1.43af075c472b6p-997
conic -1000 0x1.5a6cd73214e5fp+1 -0x1.215626306eb88p-1 0x1.0000000000000p+1 0.0 0x1.0000000000000p+0 0x1.99d21745c6b58p-999
conic -1000 0x1.6cd6a07a5f20fp+1 0x1.f11685d3bb078p-1 0.0 0.0 0x1.ef90f052e5361p-1001 0x1.8cfc1b64ac89cp-998
bow 3 0x1.473f1b8a2d980p-5 0x1.1f13235083abep+0 0x1.82774f9f81579p-1 0x1.65b1671796955p+0 0x1.32d69fcabed8bp-2 0x1.5edde88c5dcdap+0
bow 3 0x1.46213b60defefp+1 -0x1.49f101fb36dcdp-1 0x1.94c3eadff8e6ep+0 0x1.11e8257b26f3dp+0 0x1.d3ad73b2bd191p-1 0x1.34fc169f3801ap+2
bow 8 0x1.3de676df87887p+1 0x1.14d07e6677fbep-1 0x1.7ab5522fe16c5p+0 0x1.740369cd91ee4p+0 0x1.ecd43371c0877p-1 0x1.f0f6ef8f75471p+0
bow 8 0x1.44aba44a434cap+0 0x1.d22749c5568b8p+0 0x1.17580a46f2c23p-1 0x1.78502bb56b5c5p+0 0x1.d0ba14759c8f4p-5 0x1.54533fc1573cap+0 0x1.81f4ec638e935p+0 0x1.6ab63e6af2954p+0 0x1.f51e65aa292cfp-1 0x1.1eddb74b4f734p+2
bow 60 -0x1.0592688481f64p-1 -0x1.6a73e556bec5cp-1 0x1.fffffff94dfa9p-2 0x1.7ffffff7fbfc8p+0 0x1.feff9efd53117p-33 0x1.1aa923116253bp+3
bow 60 -0x1.904366a82a508p-2 0x1.0621104b08ea5p+1 0x1.7ffffffb8a241p+0 0x1.7ffffffe65ef4p+0 0x1.fffffff602e7bp-1 0x1.270631fec183fp+0
bow -60 0x1.b1eb2e4f0e708p-2 0x1.4e0859a267938p-1 0x1.ec3d500dbeaabp-91 0x1.187e2978c830bp-89 0x1.f2a7662477f00p-32 0x1.576ae42a5eabcp-58 0x1.0000000000000p+1 0x1.67f26fa21d6e3p-88 0x1.fffffff600607p-1 0x1.7888f0f9e739fp-57
bow -60 0x1.ff76b047fdf60p-2 0x1.a85acf9eaa152p+0 0x1.0000000000000p+1 0x1.793b3d123d038p-89 0x1.fffffffac2bbap-1 0x1.8478b6e396061p-58
bow 1000 -0x1.baea9a3b0e260p-4 0x1.c68fcfe1b0e56p+0 0x1.8000000000000p+0 0x1.8000000000000p+0 0x1.0000000000000p+0 0x1.1378b97254dadp+0
bow 1000 0x1.5e13da6cf8c1ep+1 0x1.a9699a2488360p-4 0x1.8000000000000p+0 0x1.8000000000000p+0 0x1.0000000000000p+0 0x1.3603edd7463fdp+1
bow -1000 0x1.6dcc4561a3850p+0 0x1.2a9939a14d61ap+1 0.0 0.0 0x1.d493df7c95d80p-502 0x1.51bb9adff0586p-998
bow -1000 0x1.17bab31544faap+1 0x1.5aec5aa504a00p-9 0.0 0.0 0x1.8e13439921390p-496 0x1.d070040d51eb1p-989
paired 3 -0x1.34eeee81ec373p-1 0x1.b78e96c922048p-3 0x1.0d0ec176d0263p-1 0x1.7347c3f19fed8p+0 0x1.ca07c5092f10dp-2 0x1.6b3a230ccc77ap-2 0x1.7bc320af96c0fp+1 0x1.f700693925960p-1 0x1.ec85a39d35761p+0 0x1.885c77ccf17eep-1
paired 3 0x1.0ee1acd0a1140p-6 0x1.0f0de1a4966f8p+0 0x1.159739e12a1c8p-1 0x1.773ee9ca3c59bp+0 0x1.08b662d1670c1p-1 0x1.386bbc086d0fbp-2 0x1.72badaf9b7cc8p+1 0x1.f4d097d906533p-1 0x1.dc6a8fb4cf215p+0 0x1.8c443c47bfdd5p+0
paired 8 0x1.436205d70bb98p+0 -0x1.5e1f63ad65ba0p-2 0x1.7fffcc21b4eb6p+1 0x1.ffff255fe9b34p-1 0x1.ffcbcdb8706b7p+0 0x1.94aebdeb84718p-7
paired 8 0x1.38c5aa8ecc58cp+0 -0x1.78ad117c7dc9cp-3 0x1.7fffbaf5df24cp+1 0x1.fffef42583636p-1 0x1.ffbec9a85d6f2p+0 0x1.ab036c90421cfp-7
paired 60 0x1.635c5d18cc74ep+1 0x1.11e0607a6bdc6p+0 0x1.0000000000001p-1 0x1.8000000000000p+0 0x1.a63db81960f99p+0 0x1.c2c8ba1a540fcp-51
paired 60 0x1.1d0bb5fdb6cf0p+0 -0x1.04411a256d048p-2 0x1.8000000000000p+1 0x1.0000000000000p+0 0x1.0000000000000p+1 0x1.a9fd9b556693ep-59
paired -60 0x1.d2ed8a328acfcp+1 0x1.7f116c06814bap+0 0x1.8fefa8d92ee40p-113 0x1.bd8823e734412p-114 0x1.6ea8df0faba6ep-58 0x1.8004b4e45b024p-54 0x1.4000000000000p+1 0x1.0000000000000p-2 0x1.cc5940652a204p+0 0x1.c3192009f85b0p-58
paired -60 0x1.c3d89ca5c8b74p+0 0x1.2aa80d3246c1fp+1 0x1.69df1d4e67350p-119 0x1.39708dfd3fd9ap-118 0x1.8ea9929c09581p-61 0x1.45b6504735618p-57
paired 1000 0x1.ae68037a5ce5cp-1 -0x1.cc60cbeb7436cp-1 0x1.8000000000000p+1 0x1.0000000000000p+0 0x1.0000000000000p+1 0x1.836b038b12707p-999
paired 1000 -0x1.1c1667d6b6420p-4 -0x1.b22ef23925cd0p-1 0x1.8000000000000p+1 0x1.0000000000000p+0 0x1.0000000000000p+1 0x1.ba99d516f9524p-999
paired -1000 0x1.83869d11da320p+0 0x1.2e32d6a5c7c46p+1 0.0 0.0 0x1.14b914eb8151cp-1001 0x1.099792c7d0567p-997
paired -1000 0x1.c0f8918fde7b0p-1 -0x1.6a033a723ec37p-1 0x1.4000000000000p+1 0x1.0000000000000p-2 0x1.b2b2682b95b16p+0 0x1.8bcbbdf9d45fep-999
waves 3 0x1.d293284200844p+1 0x1.10d29ec6625c4p-1 0x1.07e810741778cp+0 0x1.e16f6f1166738p-1 0x1.6e0765ead6aabp-1 0x1.e395b01f64fc7p-3 0x1.9e1585883131cp+0 0x1.7e2de3de110ccp-1 0x1.5cb6afb4cb9bap+0 0x1.443348bf6c797p+1 0x1.8657ad7e71dc8p+1 0x1.d4b35461bd288p-1 0x1.4a4e0be3ed2cfp+1 0x1.bd15430b04581p-2
waves 3 0x1.d38e09c8a6ad6p+1 -0x1.937fa69e0f706p-1 0x1.11e229262fc3fp+0 0x1.d90f0d743d146p-1 0x1.beddd34f79087p-1 0x1.506a43d75e611p-2 0x1.23740c13ed22ep+0 0x1.bbf21d5ed2464p-1 0x1.0dc0a3ed827cbp+0 0x1.832597090fa54p-1
waves 8 0x1.cb7dc8c482b4cp+1 0x1.909a631df5e60p-3 0x1.0058ec5d0367ap+0 0x1.fef21dd535dc9p-1 0x1.83b42befab998p-1 0x1.113269f9fa2dbp-7 0x1.0111a62076e9ep+0 0x1.fdf05c9edb9b0p-1 0x1.0233d10889f88p+0 0x1.b4c5a1fd68258p-5
waves 8 0x1.097de1f989a56p+1 0x1.17fbd3e570ca8p+0 0x1.001cd22665f5dp+0 0x1.feff43d4a222ap-1 0x1.47a292dc421bap-1 0x1.0458540afe176p-7 0x1.7ff2b887f48b7p+1 0x1.fefed45c124ebp-1 0x1.2ea612c864678p+1 0x1.05903f39fb341p-7
waves -60 0x1.56fcd8ac5fd94p+0 0x1.22f691de567bap+1 0x1.0000000000000p+2 0x1.5cb7c6acf5430p-122 0x1.8000000000000p+1 0x1.8a118b669342dp-59
waves -60 0x1.07926653242e4p+1 0x1.df53005565042p-1 0x1.eb41354f4bb79p-119 0x1.331453e507c3ep-119 0x1.331453e507c3ep-60 0x1.354565d6684aep-58 0x1.0000000000000p+2 0x1.12ba29803ca46p-119 0x1.8000000000000p+1 0x1.269f1bd3b364dp-58
waves 1000 0x1.9c21efb87d8c4p+0 0x1.292a2f5659b0ap+0 0x1.0000000000000p+0 0x1.0000000000000p+0 0x1.2b82b45f13f77p-1 0x1.1d55f4f441ffdp-999 0x1.8000000000000p+1 0x1.0000000000000p+0 0x1.2f1f7ae114409p+1 0x1.073540694839bp-999
waves 1000 0x1.b179679bc3af4p-1 0x1.61f3ac823cdf0p+0 0x1.8000000000000p+1 0x1.0000000000000p+0 0x1.319625e84d2aap+1 0x1.0e89861f01ec3p-999
waves -1000 0x1.17e747b4c8cc5p+2 0x1.3feb1567373c2p+1 0.0 0.0 0x1.7feb74fcf50bap-1001 0x1.01f98e72dd817p-998
waves -1000 0x1.b80549e8a1a4ap+0 0x1.e0881ba2559ccp+0 0x1.0000000000000p+2 0.0 0x1.8000000000000p+1 0x1.92e5a92b4edfdp-999
waves 60 0x1.045eb1c9f9070p+2 0x1.22de9e52d3ee8p+1 0x1.0000000000000p+0 0x1.0000000000000p+0 0x1.10b7bd1559bbfp-1 0x1.4842e15adf335p-59
waves 60 0x1.fd58c322e34c6p+1 0x1.2ae1b5712e334p+1 0x1.0000000000000p+0 0x1.0000000000000p+0 0x1.0a56ea4b52a22p-1 0x1.55aa9552f25d4p-59
)";

  /**
   *  Stiff splines cut to a box, each answered in exact rational arithmetic and moved less than two
   *  steps of a double by a step of the question.
   */
  inline constexpr std::string_view Clips = R"(
conic 3 -0x1.205a02e23b202p-2 0x1.8e53441ddf256p-2 0x1.24e8a8552996ap+0 0x1.6d17f29554f02p+0 0x1.39f62c3e5ea07p-5 0x1.24d45f6c2a1bfp+3 0x1.83d4f754c3885p-1 0x1.117a7ae83e91bp+0
conic 3 0x1.489996873cedep-1 -0x1.737d9d9e1f72cp-2 0x1.0870f57f2e91dp+1 0x1.4798eccbbff44p-2 0x1.f15d492b5fa08p-1 0x1.6395ca84bc157p+3 0x1.0000000000000p+0 0x1.6a09e667f3bcdp+4
conic 8 0x1.a9a47b8490dd0p-1 0x1.7cca3b93b35a0p-2 0x1.5142946610e02p+1 0x1.84bc445ab9f16p-1 0x1.fcdd92caf15f9p-1 0x1.53e0c209285c6p+5 0x1.ff689e76e4651p-1 0x1.1e594440d1f19p+8
conic 8 -0x1.80e5bb42edcf8p-3 0x1.f3e6e102251f0p-4 0x1.49f3142345eacp-1 0x1.6104b792817a7p-1 0x1.1c9eb3ff7ab32p-12 0x1.1736502725f02p+9 0x1.ce53d025e8ff8p-9 0x1.70c9476cab1a9p+6
conic 60 0x1.54fadb181c3c0p-1 0x1.a91a460d9b208p-2 0x1.a03f6e8319377p+0 0x1.131460e1773f1p+0 0x1.fe69b15300479p-61 0x1.432579f81a94cp+58 0x1.0000000000000p+0 0x1.ef5bd7f405862p+59
conic 60 -0x1.02b77509a9150p-4 0x1.d21ee98aac7a4p-2 0x1.5a6206fa89a83p+0 0x1.b8a0b0a9755d2p+0 0x1.abc9a019da334p-62 0x1.add428aa6380fp+59 0x1.0000000000000p+0 0x1.69067ff40ec54p+58
conic -60 0x1.01966663bec9fp+0 -0x1.554cd898e061fp-2 0x1.746774ea6b948p+0 0x1.1cd354463b924p-1 0x1.00cb33b1e6249p-1 0x1.fffc3841d0c2ep+1 0x1.3d918c31b758ep-1 0x1.af133f0800f45p+1
conic -60 0x1.54356ce59eb90p-3 -0x1.bb9092d61aea6p-3 0x1.097237f8aeb72p+1 0x1.f6d971724fad2p-1 0x1.d9c9f8f388d53p-3 0x1.b691ec92d753bp+0 0x1.0000000000000p+0 0x1.6a09e667f3bcdp-59
conic 1000 0x1.19f198993f530p-3 0x1.2a72c68c01f28p-3 0x1.63991f4f6b8e3p+0 0x1.174c070cd0934p-1 0x1.5d5bfd68e79b8p-1004 0x1.0835a88e45b08p+1001 0x1.33423a6c0f0c1p-1001 0x1.2b24a5750b44fp+999
conic 1000 0x1.a05b7b3838c50p-5 -0x1.bdce0792ab3eep-3 0x1.932a9fdd00e44p+0 0x1.33807e3c006b5p+0 0x1.b6a6da51f2fa1p-1006 0x1.462c38b76d913p+1001 0x1.0000000000000p+0 0x1.de9408313a73bp+999
conic -1000 -0x1.f949402cd7720p-2 -0x1.a14a9baaf5219p-2 0x1.8097371226880p-1 0x1.09aee4bd952eap+0 0.0 0x1.6a09e667f3bcdp-999 0x1.bf47007e6d772p-2 0x1.e81826429944ep+1
conic -1000 -0x1.eed1f9d6f31b4p-3 -0x1.4dda7f6654bedp-2 0x1.8d021818270acp+0 0x1.c4795150c9702p-1 0.0 0x1.6a09e667f3bcdp-999 0x1.4cdbf80462fcbp-1 0x1.87fee8e2c832ap+1
bow 3 0x1.05e8709e35280p-6 0x1.c8040b060bc64p-1 0x1.919bb76dc30f2p+0 0x1.1ffa4348adf1ap+0 0x1.be232d245a170p-5 0x1.03eb4dba3f394p+3 0x1.9e27d0b707cbdp-4 0x1.0064dfa96581ep+2 0x1.cc3b05e91f068p-1 0x1.0064dfa96581ep+2 0x1.d0f5251e8ea52p-1 0x1.1fc46444e3b27p+2
bow 3 0x1.643ba8f94911ep-1 0x1.64371b951e3c8p-3 0x1.3bf8fd5047848p+1 0x1.d4d1d2b2fcd80p-1 0x1.e2537d5610389p-1 0x1.e6dffe1b5a4e6p+2 0x1.fd3b63cc95ef5p-1 0x1.e6065c087765ap+4
bow 8 0x1.14c11de6d02d4p+0 0x1.2da198f987c60p-4 0x1.296474cb1db8cp+1 0x1.4bbb6431318c8p+0 0x1.fbd65ca18d3c3p-1 0x1.7928ce960371dp+4 0x1.fff730a17d3cep-1 0x1.12958c1ccd8eap+10
bow 8 0x1.66f6f80723580p-4 0x1.367379cd19824p-1 0x1.2a453f6d2c21bp+0 0x1.e2077f236d75cp+0 0x1.ce60fcfa0515cp-11 0x1.b0bc0c5633366p+8 0x1.53f166f1b2e37p-1 0x1.02338a997f254p+0
bow 60 -0x1.cb98280c2b891p-2 -0x1.faaaa4fd25903p-2 0x1.1e15561b50ff2p+0 0x1.a87d8fd513d6ep-1 0.0 0x1.2f9422c23c47ep+62 0x1.a5cd1f325a68fp-62 0x1.e5dda62d006dep+59
bow 60 -0x1.7542a70f80420p-5 -0x1.e38745543ef18p-2 0x1.1d48ef7b2acf2p-1 0x1.41f953b654d54p-1 0.0 0x1.2f9422c23c47ep+62 0x1.eccc5b24e879ep-63 0x1.999231276ef52p+60
bow -60 0x1.1ebfbaace9260p-3 -0x1.9110cc94e8440p-5 0x1.2a281cd523af5p-1 0x1.b3a0f0577e3aep-3 0x1.3001a2f534445p-2 0x1.df1293644362ap+0 0x1.b4a1186805ca7p-2 0x1.44094381eade4p+2
bow -60 0x1.516fad235e056p-1 -0x1.548045e67c550p-5 0x1.4817e6b60e5fcp+0 0x1.877efe9dc725dp-1 0x1.c3aadb6a28ab4p-2 0x1.582421850a2e2p+2 0x1.189f425571d82p-1 0x1.64d8e3f1bc051p+2
bow 1000 -0x1.c8911e409df58p-5 -0x1.cd25696e57000p-2 0x1.a0399de8729e7p+0 0x1.9a17ec97ed144p-1 0.0 0x1.2f9422c23c47ep+1002 0x1.871a8a766761fp-1002 0x1.07b918ed77d50p+1000
bow 1000 -0x1.7b3fb21bd7c26p-3 0x1.8a498718c0580p-4 0x1.4a4c339efdc10p-1 0x1.9658066525d47p-1 0x1.7682f5d610441p-1006 0x1.09dd6ebb30535p+1002 0x1.7f82861c710d1p-1002 0x1.0d471bd502954p+1000
bow -1000 0x1.bd1b789ee535ap-1 -0x1.7ef15f6558f8ap-3 0x1.8b32637b02504p+0 0x1.16b2aaac97a7bp+0 0x1.e996ad599744cp-2 0x1.7a2b81544def7p+2 0x1.334d6545bb8f8p-1 0x1.19c967af5db37p+2
bow -1000 0x1.30452091d9971p+0 -0x1.5e80f30e00725p-2 0x1.4c73d2400d4e7p+1 0x1.551cf647be5b4p-1 0x1.10434fb44047ap-1 0x1.73d92d441c61ep+2 0x1.0000000000000p+0 0x1.2f9422c23c47ep-998
paired 3 -0x1.66d747e4ee6e0p-2 0x1.174397c918f3cp-3 0x1.aba255ceac374p-2 0x1.8616f731dc20ap-2 0x1.0f4aca3f740ccp-8 0x1.fc63fa4eaefd8p+4 0x1.c6ba0aa4a5ed9p-7 0x1.612b7c58b4e0dp+4
paired 8 -0x1.3d56eaab0023fp-2 -0x1.c28241cfec090p-4 0x1.7124d0d2a6b75p+0 0x1.11696bc1b5058p-1 0.0 0x1.2f9422c23c47ep+10 0x1.78f9b1c51f9a5p-11 0x1.f8b3b5544bdc4p+8
paired 8 0x1.f6bb21d99a698p+0 0x1.3581c748f3e94p-1 0x1.e3e39f27886d9p+1 0x1.435a439a67e7ap+0 0x1.d1c92f162097fp+0 0x1.80f17b7901803p+3 0x1.0000000000000p+1 0x1.5a22073490377p-7
paired 60 0x1.9ad9cede9e6b0p+0 0x1.20e770f3dc996p-1 0x1.79f4085d92ddap+1 0x1.662a0297bbb9bp+0 0x1.ffffe48c5a330p+0 0x1.197599edcd64cp+20 0x1.fffff93120e21p+0 0x1.553d599ee890cp+18
paired 60 -0x1.26ceee82f26a0p-4 -0x1.d59c1190a145ep-3 0x1.906dec872d929p+0 0x1.3e6746eb29c13p+0 0.0 0x1.2f9422c23c47ep+62 0x1.9e34779367580p-60 0x1.1b7aa7bf49dbdp+57
paired -60 0x1.c51af54c34636p-1 -0x1.3c478a9829dfcp-4 0x1.3d870087c84a4p+1 0x1.15bd1735c682cp-1 0x1.4c89f2cd2c84dp-20 0x1.53b1f982783f2p+20 0x1.00466743874f6p-17 0x1.d8ecd74cf469fp+12
paired 1000 -0x1.6b1cd0a431536p-3 0x1.c67b9f34dc95cp-3 0x1.0256e435686cep-1 0x1.7b7205b5d737ep-1 0x1.da2113851e77cp-1005 0x1.b8cbd5157b51ap+1001 0x1.4d54c826d5a9bp-1002 0x1.36d2a44220dd5p+1000
paired 1000 0x1.09a678e781c80p-7 0x1.38d94b8b2c904p-2 0x1.df8aefe3f6690p-1 0x1.a81017b19c010p+0 0x1.5d371306a1b88p-1004 0x1.81042a46f817ap+1001 0x1.0000000000000p+1 0x1.0af8f9e5e03ddp+332
paired -1000 0x1.14ffae0fb3513p+0 0x1.5ce87cd7cdd20p-5 0x1.d22e0f35c426ap+1 0x1.4fa8c829d646ap+0 0x1.26bd935c7a83cp-333 0x1.9b6c26b8a1da8p+333 0x1.0000000000000p+1 0x1.5a22073490377p+1001
paired -1000 0x1.476f6e39db55ep+0 -0x1.a8e4c10ce9400p-6 0x1.9f3685fe77dbep+1 0x1.1cfae7d7c5130p+0 0x1.47936707bd9c6p-333 0x1.78c96513a61e8p+333 0x1.0000000000000p+1 0x1.5a22073490377p+1001
paired -60 0x1.a87794767478cp-1 -0x1.a79adce0c161ep-2 0x1.104821a9057d2p+1 0x1.60c67181e1346p-2 0x1.41b540b99061cp-20 0x1.5454b32de5550p+20 0x1.6b1572a447833p-19 0x1.593f15723b4dcp+18
paired 3 -0x1.2e199830be7e8p-5 0x1.909b6548165f8p-1 0x1.0bc5d003ef48ep-1 0x1.72ecf3fa6ecbcp+0 0x1.5ba8974b461b2p-5 0x1.3ee4bbb5c1c92p+3 0x1.bf4f75c540defp-2 0x1.767cfbf387a8ap-2
waves 3 0x1.11a9438247c77p+1 -0x1.55ff6e4c1b0b7p-2 0x1.dcad7aaa40163p+1 0x1.27690d9e3f468p-1 0x1.761aab6a5605ap+1 0x1.3c873833ecbd1p+2 0x1.7d0b45b4dda8dp+1 0x1.93438b0d9f9c6p+3
waves 3 -0x1.d6fb1f9471b9dp-2 0x1.4caca28dac574p-3 0x1.b6e5c4e274897p+1 0x1.3dd454d7a691ap+0 0x1.874c1d6467b8bp-7 0x1.05e13b4ad76efp+4 0x1.766b4115f0e06p+1 0x1.478af387d3f54p+2
waves 8 0x1.64ce60d8e1ba4p+1 -0x1.27c6953e973fbp-2 0x1.8e00f8b36b7f1p+2 0x1.e099c0612c4c4p-1 0x1.7c4f3fbe44e2ep+1 0x1.76c1dfd781219p+1 0x1.8000000000000p+1 0x1.6a09e667f3bcdp+9
waves 8 0x1.f655567708fd2p+0 -0x1.b6bf2a9587afep-3 0x1.2392a4b171ac0p+2 0x1.3165bbab05a64p+0 0x1.7d8ddc046a558p+0 0x1.f9862969f3075p+1 0x1.8000000000000p+1 0x1.6a09e667f3bcdp+9
waves 60 0x1.7bb211beb447cp+1 0x1.94fd07ceeae08p-2 0x1.9b4d4fd5b9e8ep+2 0x1.6d7c6ee120993p+0 0x1.e264cb42fca7ep+0 0x1.4b088858ec463p-1 0x1.8000000000000p+1 0x1.089879540c21dp+60
waves 60 0x1.fd2a73ec08170p+0 -0x1.a58233a7f5766p-3 0x1.4091d81cc45d0p+2 0x1.2c3b47f8129d1p+0 0x1.7f4a9b8ebcd81p+0 0x1.fff3f37bde918p+1 0x1.8000000000000p+1 0x1.6a09e667f3bcdp+61
waves -60 0x1.09c489afb8516p+1 -0x1.51a088ba9cd4cp-4 0x1.48afd017b215ap+2 0x1.048072e0d19bap+0 0x1.0fcd2ae2e8e9dp+1 0x1.5b4ab504ca31ap+0 0x1.8000000000000p+1 0x1.6a09e667f3bcdp-59
waves -60 0x1.1882caf81ddfcp-1 -0x1.6c654908f2afbp-2 0x1.a3f47160458bdp+1 0x1.7fc9096116ea6p-3 0x1.dc03ee6c953eep-2 0x1.995c0de73e3bfp+1 0x1.3e25abe6b4e76p+1 0x1.d7ff42de9777dp+1
waves 1000 0x1.43dbc92c7ea46p-1 0x1.331047aa78f06p-1 0x1.62228f9ca1b6cp+1 0x1.ad8cd81b70586p+0 0x1.b8aad2c261340p-1001 0x1.87166d7c37b5ep+998 0x1.bbc3b2eee291ap+0 0x1.0de225761666dp+1
waves 1000 0x1.05f4ea9b7dcb2p-1 0x1.89b866c803ffap-2 0x1.4f76c5ce580bep+1 0x1.581bf009a5d9ap+0 0x1.0c327cf469610p-1001 0x1.5962f6f0c43b6p+999 0x1.ac8b28c72fde2p+0 0x1.66077aa9c552ap+1
waves -1000 0x1.9fc129b17daa8p-1 -0x1.16c9475cc2f15p-2 0x1.b21f963d8dc91p+0 0x1.eaa6e370b8c34p-4 0x1.13f8464e0163ep-1 0x1.f0ee853904219p+1 0x1.8a0034ae8ad68p-1 0x1.745169264f8acp+1
waves -1000 0x1.516aafc1d51a8p+1 -0x1.9c3bd0ef6e524p-3 0x1.818dcada914e0p+2 0x1.0b742cd61d272p+0 0x1.29ae7b8823b7ap+1 0x1.f9aef197f11f6p+1 0x1.8000000000000p+1 0x1.6a09e667f3bcdp-999
)";

} // namespace Geoxy::Tests

#endif /* GEOXY_TESTS_STIFF_HPP */
