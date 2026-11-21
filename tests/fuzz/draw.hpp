/***************************************************************************************************
 *  @file      tests/fuzz/draw.hpp
 *  @brief     The fuzzer's bytes read as the doubles and points the exact predicates find hardest.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TESTS_FUZZ_DRAW_HPP
#define GEOXY_TESTS_FUZZ_DRAW_HPP

#include "exact/triangle.hpp"

#include "geoxy/point.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>

namespace Geoxy::Fuzz
{
  /**
   *  Four points: the three corners of a circle, then the point asked about.
   */
  struct Quad
  {
    Exact::Triangle Circle; /**< The corners the circle passes through. */
    Point           Asked;  /**< The point asked about.                 */
  };

  /**
   *  The fuzzer's bytes read as doubles and points, a byte choosing how each is made.
   *
   *  Random doubles almost never land where an exact predicate works: its first stage decides them.
   *  So most draws are built to reach the later stages: points on a line or a circle nudged by a
   *  few steps of a double, coordinates differing in their low bits, mixed magnitudes, survey-scale
   *  offsets, and the edges of Shewchuk's limits on both sides. A few are any bits at all, NaN and
   *  infinity among them.
   *
   *  The bytes are read round and round, so every input draws as much as it is asked, and an empty
   *  one draws zeros.
   */
  class Draw
  {
  public:
    /**
     *  A draw from the fuzzer's bytes.
     *
     *  @param bytes The bytes, outliving the draw.
     */
    explicit Draw(std::span<const std::uint8_t> bytes) : m_Bytes(bytes) {}

    /**
     *  The next byte.
     *
     *  @return The byte, or 0 when there are none.
     */
    auto Byte() -> std::uint8_t;

    /**
     *  A coordinate, made the way the next byte says.
     *
     *  @return The coordinate, possibly NaN, infinite or outside Shewchuk's limits.
     */
    auto Coordinate() -> double;

    /**
     *  A double moved up to eight steps either way, as many as the next byte says.
     *
     *  @param  value The double.
     *  @return       The double the steps away, the double itself for none.
     */
    auto Nudged(double value) -> double;

    /**
     *  A point of two coordinates.
     *
     *  @return The point.
     */
    auto Free() -> Point;

    /**
     *  A point on the line through two points, rounded, then nudged.
     *
     *  @param  first  One point of the line.
     *  @param  second The other.
     *  @return        The point.
     */
    auto Near(Point first, Point second) -> Point;

    /**
     *  Three corners, free, nearly collinear or nearly coincident, then moved as a whole.
     *
     *  @return The corners.
     */
    auto Corners() -> Exact::Triangle;

    /**
     *  Four points, free, nearly cocircular or with collinear corners, then moved as a whole.
     *
     *  @return The points.
     */
    auto Cocircular() -> Quad;

    /**
     *  A point moved as the whole of its figure is: swapped, turned a quarter, scaled by a power of
     *  two, offset, or left.
     *
     *  @param  point The point.
     *  @param  how   The byte the figure drew for its move.
     *  @return       The point moved.
     */
    static auto Moved(Point point, std::uint8_t how) -> Point;

  private:
    /**
     *  A whole number of two bytes.
     *
     *  @return The number, from -32768 to 32767.
     */
    auto Short() -> int;

    /**
     *  Eight bytes.
     *
     *  @return The bytes as a word, the first the lowest.
     */
    auto Word() -> std::uint64_t;

    /**
     *  A coordinate drawn before, the next byte saying which.
     *
     *  @return The coordinate, or 0 before any.
     */
    auto Earlier() -> double;

    std::span<const std::uint8_t> m_Bytes;
    std::size_t                   m_At    = 0;
    std::array<double, 16>        m_Drawn = {};
    std::size_t                   m_Count = 0;
  };
} // namespace Geoxy::Fuzz

#endif /* GEOXY_TESTS_FUZZ_DRAW_HPP */
