/***************************************************************************************************
 *  @file      geoxy/tree.hpp
 *  @brief     Curves kept in a tree of their boxes, put in and taken out one at a time: those near
 *             a box, those nearest a point and those a ray strikes, found without asking each.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_TREE_HPP
#define GEOXY_TREE_HPP

#include "geoxy/api.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"
#include "geoxy/ray.hpp"
#include "geoxy/tolerance.hpp"

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
   *  Curves in an R-tree of their boxes: each curve's box in a leaf, each leaf's box in a limb, up
   *  to one root, so that a box, a point or a ray is asked of the few limbs it comes near and of
   *  the curves in them, not of every curve.
   *
   *  Where a Grid is made once for curves of like size and met all at once, a Tree takes curves of
   *  any sizes, one more or one fewer at a time. Every answer is the one asking every curve gives:
   *  the tree only spares the asking. A curve with no box, a ray, an endless line, one of numbers
   *  that are not finite, is kept apart and asked every time.
   *
   *  The tree holds its own curves. A curve is told by its index: its place among the curves the
   *  tree was made of, or what Insert gave; an index erased is given again to a later curve.
   */
  class Tree
  {
  public:
    /**
     *  A curve near a point: which, where on it, and how far.
     */
    struct Close
    {
      std::size_t Index    = 0;   /**< The curve's index.                    */
      Place       Nearest  = {};  /**< The place on it nearest the point.    */
      double      Distance = 0.0; /**< How far that place is from the point. */

      /**
       *  Whether two are equal exactly.
       *
       *  @param  lhs One curve.
       *  @param  rhs The other.
       *  @return     Whether they are equal.
       */
      [[nodiscard]] friend constexpr auto operator==(const Close &lhs, const Close &rhs) -> bool
      {
        return lhs.Index == rhs.Index
            && lhs.Nearest == rhs.Nearest
            && std::equal_to<double>{}(lhs.Distance, rhs.Distance);
      }
    };

    /**
     *  Where a ray strikes a curve.
     */
    struct Strike
    {
      std::size_t Index = 0;   /**< The curve's index.         */
      Point       At    = {};  /**< Where.                     */
      double      Along = 0.0; /**< The ray's parameter there. */

      /**
       *  Whether two are equal exactly.
       *
       *  @param  lhs One strike.
       *  @param  rhs The other.
       *  @return     Whether they are equal.
       */
      [[nodiscard]] friend constexpr auto operator==(const Strike &lhs, const Strike &rhs) -> bool
      {
        return lhs.Index == rhs.Index
            && lhs.At == rhs.At
            && std::equal_to<double>{}(lhs.Along, rhs.Along);
      }
    };

    /**
     *  How far a point may be asked for at most: any distance.
     */
    static constexpr double Anywhere = std::numeric_limits<double>::infinity();

    /**
     *  A tree of no curve.
     */
    GEOXY_API Tree();

    /**
     *  A tree of curves, copied and packed all at once, which is faster and packs closer than
     *  putting them in one at a time.
     *
     *  @param curves The curves, in the order their indices give.
     */
    GEOXY_API explicit Tree(std::span<const Curve> curves);

    /**
     *  A tree of the same curves at the same indices.
     *
     *  @param other The tree.
     */
    GEOXY_API Tree(const Tree &other);

    /**
     *  A tree taken over.
     *
     *  @param other The tree, left to be assigned to or to go.
     */
    GEOXY_API Tree(Tree &&other) noexcept;

    /**
     *  The same curves as another tree, at the same indices.
     *
     *  @param  other The tree.
     *  @return       This tree.
     */
    GEOXY_API auto operator=(const Tree &other) -> Tree &;

    /**
     *  Another tree taken over.
     *
     *  @param  other The tree, left to be assigned to or to go.
     *  @return       This tree.
     */
    GEOXY_API auto operator=(Tree &&other) noexcept -> Tree &;

    /**
     *  The tree goes, and its curves with it.
     */
    GEOXY_API ~Tree();

    /**
     *  How many curves it holds.
     *
     *  @return The count.
     */
    [[nodiscard]] auto Size() const -> std::size_t { return m_Count; }

    /**
     *  Whether an index is a curve's.
     *
     *  @param  index The index.
     *  @return       Whether the tree holds a curve there.
     */
    [[nodiscard]] auto Holds(std::size_t index) const -> bool
    {
      return index < m_Held.size() && m_Held[index] != 0;
    }

    /**
     *  A curve.
     *
     *  @pre          The tree holds a curve at the index.
     *  @param  index Its index.
     *  @return       The curve.
     */
    [[nodiscard]] auto operator[](std::size_t index) const -> const Curve &
    {
      return m_Curves[index];
    }

    /**
     *  A curve more.
     *
     *  @param  curve The curve.
     *  @return       Its index: the last one erased, or one past the greatest where none is spare.
     */
    GEOXY_API auto Insert(Curve curve) -> std::size_t;

    /**
     *  A curve no more. An index that is no curve's changes nothing.
     *
     *  @param index The curve's index.
     */
    GEOXY_API void Erase(std::size_t index);

    /**
     *  The curves that may lie in a box or cross it: each whose box, a few steps of a double
     *  larger, overlaps it, and each with no box. No curve that reaches the box is left out, and
     *  few others are in; what Window and Crossing choose from.
     *
     *  @param  box The box.
     *  @return     Their indices, in order.
     */
    [[nodiscard]] GEOXY_API auto Near(const Box &box) const -> std::vector<std::size_t>;

    /**
     *  The curves lying wholly within a box, as Inside tells it.
     *
     *  @param  box The box.
     *  @return     Their indices, in order.
     */
    [[nodiscard]] GEOXY_API auto Window(const Box &box) const -> std::vector<std::size_t>;

    /**
     *  The curves that meet a box, within it or across its edges, as Intersects tells it.
     *
     *  @param  box The box.
     *  @return     Their indices, in order.
     */
    [[nodiscard]] GEOXY_API auto Crossing(const Box &box) const -> std::vector<std::size_t>;

    /**
     *  The curves nearest a point, the nearest first, each with its place nearest the point as
     *  Nearest gives it: the same curves, places and distances as asking every curve gives.
     *
     *  @param  point  The point.
     *  @param  count  How many at most.
     *  @param  within How far at most; Anywhere, for any distance.
     *  @return        The nearest, by distance; of two as near, the lesser index first.
     */
    [[nodiscard]] GEOXY_API auto Nearest(Point point, std::size_t count, double within) const
        -> std::vector<Close>;

    /**
     *  Where a ray strikes the curves, the nearest its start first: every hit Intersect finds
     *  between the ray and a curve, and both ends of every stretch they share.
     *
     *  @param  ray       The ray.
     *  @param  tolerance The tolerance, as Intersect takes it.
     *  @return           The strikes, by the ray's parameter; of two at one, the lesser index
     *                    first.
     */
    [[nodiscard]] GEOXY_API auto Cast(const Ray &ray, const Tolerance &tolerance) const
        -> std::vector<Strike>;

  private:
    std::vector<Curve>        m_Curves = {}; /**< The curves, and what erased ones left. */
    std::vector<std::uint8_t> m_Held   = {}; /**< Whether each index is a curve's.       */
    std::vector<std::size_t>  m_Spare  = {}; /**< The indices erased, the last one last. */
    std::size_t               m_Count  = 0;  /**< How many curves are held.              */
    double m_Size = 0.0; /**< The largest size of the numbers of a curve it has held. */
    std::unique_ptr<Boughs> m_Boughs; /**< Their boxes, in a tree.                */
  };
} // namespace Geoxy

#endif /* GEOXY_TREE_HPP */
