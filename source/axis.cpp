/***************************************************************************************************
 *  @file      source/axis.cpp
 *  @brief     Whether a point lies exactly on a line's carrier or an arc's circle, by the sign of
 *             an exact sum of products.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "axis.hpp"

#include "real.hpp"
#include "twofold.hpp"

#include "geoxy/arc.hpp"
#include "geoxy/line.hpp"
#include "geoxy/point.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>

namespace Geoxy::Axis
{
  namespace
  {
    /**
     *  A sum of doubles held exactly: Shewchuk's expansion, its parts none overlapping another,
     *  none zero, the least first, grown one double at a time.
     */
    class Sum
    {
    public:
      /**
       *  A double added, exactly: Shewchuk's grow-expansion, zeros left out; a zero adds nothing,
       *  and is not added at all.
       *
       *  @param value The double.
       */
      void Add(double value)
      {
        if (Zero(value))
        {
          return;
        }
        double      carried = value;
        std::size_t kept    = 0;
        for (std::size_t index = 0; index < m_Count; ++index)
        {
          const Twofold::Two sum = Twofold::TwoSum(carried, m_Parts.at(index));
          if (!Zero(sum.Low))
          {
            m_Parts.at(kept) = sum.Low;
            ++kept;
          }
          carried = sum.High;
        }
        if (!Zero(carried))
        {
          m_Parts.at(kept) = carried;
          ++kept;
        }
        m_Count = kept;
      }

      /**
       *  A product of two doubles added, exactly.
       *
       *  @param first  One.
       *  @param second The other.
       */
      void Times(double first, double second)
      {
        const Twofold::Two product = Twofold::TwoProduct(first, second);
        Add(product.Low);
        Add(product.High);
      }

      /**
       *  Whether the sum is zero.
       *
       *  @return Whether it is.
       */
      [[nodiscard]] auto Nothing() const -> bool { return m_Count == 0; }

    private:
      std::array<double, 24> m_Parts = {}; /**< The parts, the least first. */
      std::size_t            m_Count = 0;  /**< How many there are.         */
    };

    /**
     *  Whether every number may be asked of.
     *
     *  @param  values The numbers.
     *  @return        Whether they may.
     */
    auto Heeded(std::initializer_list<double> values) -> bool
    {
      return std::ranges::all_of(values, Held);
    }

    /**
     *  A difference of two doubles, exactly, in two doubles.
     *
     *  @param  head The double.
     *  @param  tail The double taken from it.
     *  @return      The difference.
     */
    auto Severed(double head, double tail) -> Twofold::Two
    {
      return Twofold::TwoSum(head, -tail);
    }
  } // namespace

  auto Held(double value) -> bool
  {
    const double size = std::abs(value);
    return std::isfinite(size) && (Zero(size) || (size >= 0x1p-100 && size <= 0x1p100));
  }

  auto Fits(double size) -> bool
  {
    const double magnitude = std::abs(size);
    return magnitude > 0x1p-280 && magnitude < 0x1p280;
  }

  auto On(const Line &line, Point point) -> bool
  {
    if (!Heeded({line.Start.X, line.Start.Y, line.End.X, line.End.Y, point.X, point.Y}))
    {
      return false;
    }
    const Twofold::Two wayX = Severed(line.End.X, line.Start.X);
    const Twofold::Two wayY = Severed(line.End.Y, line.Start.Y);
    const Twofold::Two gapX = Severed(point.X, line.Start.X);
    const Twofold::Two gapY = Severed(point.Y, line.Start.Y);
    if (Zero(wayX.High) && Zero(wayY.High))
    {
      return false;
    }
    Sum cross;
    for (const double way : {wayX.High, wayX.Low})
    {
      for (const double gap : {gapY.High, gapY.Low})
      {
        cross.Times(way, gap);
      }
    }
    for (const double way : {wayY.High, wayY.Low})
    {
      for (const double gap : {gapX.High, gapX.Low})
      {
        cross.Times(-way, gap);
      }
    }
    return cross.Nothing();
  }

  auto On(const Arc &arc, Point point) -> bool
  {
    if (!Heeded({arc.Center.X, arc.Center.Y, arc.Radius, point.X, point.Y}))
    {
      return false;
    }
    Sum gap;
    for (const Twofold::Two apart :
         {Severed(point.X, arc.Center.X), Severed(point.Y, arc.Center.Y)})
    {
      gap.Times(apart.High, apart.High);
      gap.Times(2.0 * apart.High, apart.Low);
      gap.Times(apart.Low, apart.Low);
    }
    gap.Times(-arc.Radius, arc.Radius);
    return gap.Nothing();
  }
} // namespace Geoxy::Axis
