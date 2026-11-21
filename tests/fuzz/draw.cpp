/***************************************************************************************************
 *  @file      tests/fuzz/draw.cpp
 *  @brief     The fuzzer's bytes read as doubles and points.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "draw.hpp"

#include "exact/triangle.hpp"

#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace Geoxy::Fuzz
{
  namespace
  {
    /**
     *  Where survey and CAD coordinates sit: the origin, a unit square, a building, UTM eastings
     *  and northings, and far out.
     */
    constexpr std::array Offsets{0.0, 0.5, 1e3, 5e5, 5e6, -5e6, 1e8, -1e8};

    /**
     *  The edges of Shewchuk's limits on both sides, and the doubles at the ends of the range.
     */
    constexpr std::array Specials{
        0.0,
        -0.0,
        1.0,
        -1.0,
        0x1.0p-400,
        0x1.0p400,
        0x1.0p-401,
        0x1.0p401,
        0x1.0p-150,
        0x1.0p150,
        0x1.0p-151,
        0x1.0p151,
        std::numeric_limits<double>::min(),
        std::numeric_limits<double>::denorm_min(),
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::quiet_NaN(),
        0x1.fffffffffffffp-1,
        0x1.0000000000001p0,
    };

    /**
     *  A byte as a signed number.
     *
     *  @param  byte The byte.
     *  @return      The number, from -128 to 127.
     */
    auto Signed(std::uint8_t byte) -> int
    {
      return static_cast<int>(static_cast<std::int8_t>(byte));
    }
  } // namespace

  auto Draw::Byte() -> std::uint8_t
  {
    if (m_Bytes.empty())
    {
      return 0;
    }
    const std::uint8_t byte = m_Bytes[m_At];
    m_At                    = (m_At + 1) % m_Bytes.size();
    return byte;
  }

  auto Draw::Short() -> int
  {
    const auto low  = static_cast<unsigned>(Byte());
    const auto high = static_cast<unsigned>(Byte());
    return static_cast<int>(
        static_cast<std::int16_t>(static_cast<std::uint16_t>(low | (high << 8U))));
  }

  auto Draw::Word() -> std::uint64_t
  {
    std::uint64_t word = 0;
    for (unsigned shift = 0; shift < 64; shift += 8)
    {
      word |= static_cast<std::uint64_t>(Byte()) << shift;
    }
    return word;
  }

  auto Draw::Earlier() -> double
  {
    const std::uint8_t which = Byte();
    if (m_Count == 0)
    {
      return 0.0;
    }
    return m_Drawn.at(which % std::min(m_Count, m_Drawn.size()));
  }

  /**
   *  Eight ways, a byte choosing: any bits; a small whole number; a survey offset plus a small
   *  binary fraction; a double of any exponent near the range's ends; a coordinate drawn before,
   *  nudged; one drawn before plus a far smaller step, so the two differ in their low bits; an edge
   *  of Shewchuk's limits; a decimal with three places, as a drawing holds.
   */
  auto Draw::Coordinate() -> double
  {
    double value = 0.0;
    switch (Byte() % 8U)
    {
    case 0:
      value = std::bit_cast<double>(Word());
      break;
    case 1:
      value = static_cast<double>(Short());
      break;
    case 2:
    {
      const double offset = Offsets.at(Byte() % Offsets.size());
      const auto   step   = static_cast<double>(Short());
      value               = offset + std::ldexp(step, -(Byte() % 48));
      break;
    }
    case 3:
    {
      const std::uint64_t word     = Word();
      const double        fraction = std::ldexp(static_cast<double>(word >> 12U), -52);
      const double        sign     = (word & 1U) == 0 ? 1.0 : -1.0;
      value                        = sign * std::ldexp(1.0 + fraction, Signed(Byte()) * 5);
      break;
    }
    case 4:
      value = Nudged(Earlier());
      break;
    case 5:
    {
      const double earlier = Earlier();
      const auto   step    = static_cast<double>(Short());
      value                = earlier + std::ldexp(step, -(Byte() % 64));
      break;
    }
    case 6:
      value = Specials.at(Byte() % Specials.size());
      break;
    default:
    {
      const double whole = static_cast<double>(Short()) * 1000.0;
      value              = whole + (static_cast<double>(Short()) / 1000.0);
      break;
    }
    }
    m_Drawn.at(m_Count % m_Drawn.size()) = value;
    ++m_Count;
    return value;
  }

  auto Draw::Nudged(double value) -> double
  {
    int steps = static_cast<int>(Byte() % 17U) - 8;
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

  auto Draw::Free() -> Point
  {
    const double x = Coordinate();
    return {.X = x, .Y = Coordinate()};
  }

  /**
   *  The fraction along the line is a whole number, a binary fraction, a half or any coordinate:
   *  the first three land on the line exactly where the arithmetic is exact.
   */
  auto Draw::Near(Point first, Point second) -> Point
  {
    double fraction = 0.5;
    switch (Byte() % 4U)
    {
    case 0:
      fraction = static_cast<double>(Signed(Byte()));
      break;
    case 1:
    {
      const auto numerator = static_cast<double>(Signed(Byte()));
      fraction             = std::ldexp(numerator, -(Byte() % 16));
      break;
    }
    case 2:
      fraction = Coordinate();
      break;
    default:
      break;
    }
    const Point  online = first + ((second - first) * fraction);
    const double x      = Nudged(online.X);
    return {.X = x, .Y = Nudged(online.Y)};
  }

  auto Draw::Moved(Point point, std::uint8_t how) -> Point
  {
    const int step = static_cast<int>(how >> 3U);
    switch (how % 5U)
    {
    case 1:
      return {.X = point.Y, .Y = point.X};
    case 2:
      return {.X = -point.Y, .Y = point.X};
    case 3:
      return {.X = std::ldexp(point.X, (step - 16) * 16),
              .Y = std::ldexp(point.Y, (step - 16) * 16)};
    case 4:
    {
      const double offset = Offsets.at(static_cast<std::size_t>(step) % Offsets.size());
      return point + Vector{.X = offset, .Y = -offset};
    }
    default:
      return point;
    }
  }

  auto Draw::Corners() -> Exact::Triangle
  {
    Exact::Triangle corners;
    corners.First  = Free();
    corners.Second = Free();
    switch (Byte() % 4U)
    {
    case 0:
      corners.Third = Free();
      break;
    case 1:
      corners.Third = Near(corners.First, corners.Second);
      break;
    case 2:
    {
      const double x = Nudged(corners.First.X);
      corners.Third  = {.X = x, .Y = Nudged(corners.First.Y)};
      break;
    }
    default:
      corners.Second = Near(corners.First, corners.Second);
      corners.Third  = Near(corners.First, corners.Second);
      break;
    }
    const std::uint8_t how = Byte();
    return {.First  = Moved(corners.First, how),
            .Second = Moved(corners.Second, how),
            .Third  = Moved(corners.Third, how)};
  }

  /**
   *  Cocircular points come three ways. The eight images of one offset about a center, (u, v) and
   *  (v, u) with either sign each, lie on one circle exactly where the sums are exact, and each is
   *  nudged or not. Points of a circle by Pythagoras' rational parametrization lie on it but for
   *  rounding. Three nearly collinear corners have no circle at all.
   */
  auto Draw::Cocircular() -> Quad
  {
    std::array<Point, 4> points = {};
    switch (Byte() % 4U)
    {
    case 0:
      for (Point &point : points)
      {
        point = Free();
      }
      break;
    case 1:
    {
      const Point  center = Free();
      const double u      = Coordinate();
      const double v      = Coordinate();
      for (Point &point : points)
      {
        const std::uint8_t which = Byte();
        const double       signX = (which & 1U) == 0 ? 1.0 : -1.0;
        const double       signY = (which & 2U) == 0 ? 1.0 : -1.0;
        const bool         turn  = (which & 4U) != 0;
        const Vector       away  = {.X = signX * (turn ? v : u), .Y = signY * (turn ? u : v)};
        point                    = center + away;
        if ((which & 8U) != 0)
        {
          const double x = Nudged(point.X);
          point          = {.X = x, .Y = Nudged(point.Y)};
        }
      }
      break;
    }
    case 2:
    {
      const Point  center = Free();
      const double radius = Coordinate();
      for (Point &point : points)
      {
        const auto   numerator = static_cast<double>(Signed(Byte()));
        const double t         = std::ldexp(numerator, -(Byte() % 8));
        const double ends      = 1.0 + (t * t);
        const Vector away   = {.X = radius * (1.0 - (t * t)) / ends, .Y = radius * 2.0 * t / ends};
        const Point  online = center + away;
        const double x      = Nudged(online.X);
        point               = {.X = x, .Y = Nudged(online.Y)};
      }
      break;
    }
    default:
      points[0] = Free();
      points[1] = Free();
      points[2] = Near(points[0], points[1]);
      points[3] = (Byte() & 1U) == 0 ? Free() : Near(points[0], points[1]);
      break;
    }
    const std::uint8_t how = Byte();
    for (Point &point : points)
    {
      point = Moved(point, how);
    }
    return {.Circle = {.First = points[0], .Second = points[1], .Third = points[2]},
            .Asked  = points[3]};
  }
} // namespace Geoxy::Fuzz
