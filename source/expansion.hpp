/***************************************************************************************************
 *  @file      source/expansion.hpp
 *  @brief     Numbers held exactly as whole numbers of 64-bit limbs times a power of two, so that
 *             sums and products of doubles never round; and the double nearest a quotient that may
 *             hold a square root.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_EXPANSION_HPP
#define GEOXY_EXPANSION_HPP

#include "scratch.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>

/**
 *  Exact numbers: every double is one, and so is every sum, difference and product of them; nothing
 *  rounds until a double is asked for. A number is a whole number in two's complement, in limbs of
 *  64 bits from the least, times 2^(64 Low): the places of its limbs are whole multiples of 64, so
 *  that a sum only lines limbs up, never shifts bits, and a product only multiplies them. Where
 *  Shewchuk's expansions add each double by a chain of rounded sums that wait on one another, here
 *  a limb waits on a carry alone: a crossing's numbers, a few limbs each, are summed and multiplied
 *  some ten times sooner.
 */
namespace Geoxy::Expansion
{
  /**
   *  A number's limbs: up to 16 held in place, as nearly every number a crossing works out has, the
   *  rest, should a number have more, in the thread's pool; so that working a crossing out asks the
   *  system for no memory. Copies copy the limbs in use, not the room.
   */
  class Parts
  {
  public:
    /**
     *  How many limbs are held in place.
     */
    static constexpr std::size_t Held = 16;

    Parts() = default;

    /**
     *  A copy.
     *
     *  @param other The limbs.
     */
    Parts(const Parts &other);

    /**
     *  A copy, the other's room in the pool taken over.
     *
     *  @param other The limbs.
     */
    Parts(Parts &&other) noexcept;

    /**
     *  A copy.
     *
     *  @param  other The limbs.
     *  @return       Them.
     */
    auto operator=(const Parts &other) -> Parts &;

    /**
     *  A copy, the other's room in the pool taken over.
     *
     *  @param  other The limbs.
     *  @return       Them.
     */
    auto operator=(Parts &&other) noexcept -> Parts &;

    ~Parts() = default;

    /**
     *  How many limbs there are: none for zero.
     *
     *  @return The count.
     */
    [[nodiscard]] auto Count() const -> std::size_t { return m_Count; }

    /**
     *  Whether the number is zero.
     *
     *  @return Whether it has no limbs.
     */
    [[nodiscard]] auto Empty() const -> bool { return m_Count == 0; }

    /**
     *  The place of the least limb: it counts 2^(64 Low).
     *
     *  @return The place.
     */
    [[nodiscard]] auto Low() const -> std::int32_t { return m_Low; }

    /**
     *  The limbs.
     *
     *  @return Them, from the least, the largest's top bit the sign.
     */
    [[nodiscard]] auto Limbs() const -> std::span<const std::uint64_t>
    {
      return std::span{m_Count <= Held ? m_Held.data() : m_Heap.data(), m_Count};
    }

    /**
     *  Room for some limbs, what was held dropped: Settle tells what was put there.
     *
     *  @param  count How many.
     *  @return       The room.
     */
    auto Room(std::size_t count) -> std::span<std::uint64_t>;

    /**
     *  The limbs put in the room made the number's, at a place: the least that are zero and the
     *  largest that only repeat the sign left out, so that no number has more limbs than it needs,
     *  and zero none.
     *
     *  @param low The place of the room's least limb.
     */
    void Settle(std::int32_t low);

  private:
    std::int32_t                    m_Low   = 0;  /**< The place of the least limb. */
    std::uint32_t                   m_Count = 0;  /**< How many limbs there are.    */
    std::array<std::uint64_t, Held> m_Held  = {}; /**< The limbs, while they fit.   */
    List<std::uint64_t>             m_Heap  = {}; /**< The limbs, once they do not. */
  };

  /**
   *  A number held exactly.
   */
  struct Number
  {
    Geoxy::Expansion::Parts Parts = {}; /**< Its limbs. */
  };

  /**
   *  A double as a number.
   *
   *  @param  value The double, finite.
   *  @return       The number.
   */
  auto Of(double value) -> Number;

  /**
   *  The sum of two numbers, exactly.
   *
   *  @param  first  One.
   *  @param  second The other.
   *  @return        The sum.
   */
  auto Sum(const Number &first, const Number &second) -> Number;

  /**
   *  The difference of two numbers, exactly.
   *
   *  @param  first  The number.
   *  @param  second The number taken from it.
   *  @return        The difference.
   */
  auto Difference(const Number &first, const Number &second) -> Number;

  /**
   *  The product of a number and a double, exactly.
   *
   *  @param  number The number.
   *  @param  factor The double, finite.
   *  @return        The product.
   */
  auto Scaled(const Number &number, double factor) -> Number;

  /**
   *  The product of two numbers, exactly.
   *
   *  @param  first  One.
   *  @param  second The other.
   *  @return        The product.
   */
  auto Product(const Number &first, const Number &second) -> Number;

  /**
   *  The sign of a number.
   *
   *  @param  number The number.
   *  @return        1, -1 or 0.
   */
  auto Sign(const Number &number) -> int;

  /**
   *  The double nearest a number, the nearest even where two are as near: below the normal doubles
   *  the nearest subnormal, past the largest double's halfway point an infinity.
   *
   *  @param  number The number.
   *  @return        The double.
   */
  auto Approximate(const Number &number) -> double;

  /**
   *  The parts of a quotient that may hold a square root: (top + bottom sqrt(root)) / below.
   */
  struct Quotient
  {
    Number Top    = {}; /**< What is added to the root's share.     */
    Number Bottom = {}; /**< What the square root is multiplied by. */
    Number Root   = {}; /**< Under the square root, not negative.   */
    Number Below  = {}; /**< What all is divided by, not zero.      */
  };

  /**
   *  The double nearest (top + side bottom sqrt(root)) / below, the nearest even where two are as
   *  near: the place where two lines or two conics cross, rounded once. Double-double arithmetic
   *  settles most; else each candidate within its error bound is checked by the exact sign of how
   *  far the number lies from the midpoints either side of it. Where the double-double estimate
   *  leaves the doubles' range, its own value.
   *
   *  @param  parts top, bottom, root and below, root not negative and below not zero.
   *  @param  side  1 or -1: which root.
   *  @return       The double nearest.
   */
  auto Rounded(const Quotient &parts, int side) -> double;

} // namespace Geoxy::Expansion

#endif /* GEOXY_EXPANSION_HPP */
