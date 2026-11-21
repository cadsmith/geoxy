/***************************************************************************************************
 *  @file      geoxy/cloud.hpp
 *  @brief     Points kept in a tree, put in and taken out one at a time: those in a box and those
 *             nearest a point, found without asking each.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_CLOUD_HPP
#define GEOXY_CLOUD_HPP

#include "geoxy/api.hpp"
#include "geoxy/box.hpp"
#include "geoxy/point.hpp"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <span>
#include <vector>

namespace Geoxy
{
  struct Boughs;

  /**
   *  Points in an R-tree: the ends, middles and centres a drawing snaps to, or any others, so that
   *  a box or a point is asked of the few points it comes near, not of every one.
   *
   *  Every answer is the one asking every point gives: the tree only spares the asking. A point of
   *  numbers that are not finite is kept apart and asked every time.
   *
   *  A point is told by its index: its place among the points the cloud was made of, or what Insert
   *  gave; an index erased is given again to a later point.
   */
  class Cloud
  {
  public:
    /**
     *  A point near another: which, and how far.
     */
    struct Close
    {
      std::size_t Index    = 0;   /**< The point's index. */
      double      Distance = 0.0; /**< How far it is.     */

      /**
       *  Whether two are equal exactly.
       *
       *  @param  lhs One point.
       *  @param  rhs The other.
       *  @return     Whether they are equal.
       */
      [[nodiscard]] friend constexpr auto operator==(const Close &lhs, const Close &rhs) -> bool
      {
        return lhs.Index == rhs.Index && std::equal_to<double>{}(lhs.Distance, rhs.Distance);
      }
    };

    /**
     *  How far a point may be asked for at most: any distance.
     */
    static constexpr double Anywhere = std::numeric_limits<double>::infinity();

    /**
     *  A cloud of no point.
     */
    GEOXY_API Cloud();

    /**
     *  A cloud of points, packed all at once, which is faster and packs closer than putting them in
     *  one at a time.
     *
     *  @param points The points, in the order their indices give.
     */
    GEOXY_API explicit Cloud(std::span<const Point> points);

    /**
     *  A cloud of the same points at the same indices.
     *
     *  @param other The cloud.
     */
    GEOXY_API Cloud(const Cloud &other);

    /**
     *  A cloud taken over.
     *
     *  @param other The cloud, left to be assigned to or to go.
     */
    GEOXY_API Cloud(Cloud &&other) noexcept;

    /**
     *  The same points as another cloud, at the same indices.
     *
     *  @param  other The cloud.
     *  @return       This cloud.
     */
    GEOXY_API auto operator=(const Cloud &other) -> Cloud &;

    /**
     *  Another cloud taken over.
     *
     *  @param  other The cloud, left to be assigned to or to go.
     *  @return       This cloud.
     */
    GEOXY_API auto operator=(Cloud &&other) noexcept -> Cloud &;

    /**
     *  The cloud goes, and its points with it.
     */
    GEOXY_API ~Cloud();

    /**
     *  How many points it holds.
     *
     *  @return The count.
     */
    [[nodiscard]] auto Size() const -> std::size_t { return m_Count; }

    /**
     *  Whether an index is a point's.
     *
     *  @param  index The index.
     *  @return       Whether the cloud holds a point there.
     */
    [[nodiscard]] auto Holds(std::size_t index) const -> bool
    {
      return index < m_Held.size() && m_Held[index] != 0;
    }

    /**
     *  A point.
     *
     *  @pre          The cloud holds a point at the index.
     *  @param  index Its index.
     *  @return       The point.
     */
    [[nodiscard]] auto operator[](std::size_t index) const -> Point { return m_Points[index]; }

    /**
     *  A point more.
     *
     *  @param  point The point.
     *  @return       Its index: the last one erased, or one past the greatest where none is spare.
     */
    GEOXY_API auto Insert(Point point) -> std::size_t;

    /**
     *  A point no more. An index that is no point's changes nothing.
     *
     *  @param index The point's index.
     */
    GEOXY_API void Erase(std::size_t index);

    /**
     *  The points a box holds, its edges included.
     *
     *  @param  box The box.
     *  @return     Their indices, in order.
     */
    [[nodiscard]] GEOXY_API auto Within(const Box &box) const -> std::vector<std::size_t>;

    /**
     *  The points nearest a point, the nearest first: the same points and distances as asking every
     *  point gives.
     *
     *  @param  point  The point.
     *  @param  count  How many at most.
     *  @param  within How far at most; Anywhere, for any distance.
     *  @return        The nearest, by distance; of two as near, the lesser index first.
     */
    [[nodiscard]] GEOXY_API auto Nearest(Point point, std::size_t count, double within) const
        -> std::vector<Close>;

  private:
    std::vector<Point>        m_Points = {}; /**< The points, and what erased ones left. */
    std::vector<std::uint8_t> m_Held   = {}; /**< Whether each index is a point's.       */
    std::vector<std::size_t>  m_Spare  = {}; /**< The indices erased, the last one last. */
    std::size_t               m_Count  = 0;  /**< How many points are held.              */
    std::unique_ptr<Boughs>   m_Boughs;      /**< Their boxes, in a tree.                */
  };
} // namespace Geoxy

#endif /* GEOXY_CLOUD_HPP */
