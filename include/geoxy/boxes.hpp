/***************************************************************************************************
 *  @file      geoxy/boxes.hpp
 *  @brief     Boxes kept in a tree by the caller's tags, inserted, moved and erased one at a time:
 *             those overlapping a box, those within it, those nearest a point, and the box of them
 *             all.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_BOXES_HPP
#define GEOXY_BOXES_HPP

#include "geoxy/api.hpp"
#include "geoxy/box.hpp"
#include "geoxy/point.hpp"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <unordered_map>
#include <vector>

namespace Geoxy
{
  struct Boughs;

  /**
   *  Boxes in an R-tree, each told by a tag the caller chooses: the boxes of what a caller keeps
   *  itself, so that a box or a point is asked of the few it comes near, not of every one.
   *
   *  Only boxes are kept, never what they bound: a caller answers the exact questions itself, of
   *  the few tags a query gives. Every answer is exact, the one asking every box gives: the tree
   *  only spares the asking. A box that is empty or has a side that is not a number is kept, and
   *  overlaps nothing.
   *
   *  Made from many boxes at once, the tree is packed in one pass; after that each Insert and Erase
   *  changes one leaf and the boxes above it. Extents is the exact box of them all, read from the
   *  root. Each query has two forms: one that gives the tags in a vector, in order, and one that
   *  tells a visitor of each as it is found, allocating nothing, until the visitor says stop.
   */
  class Boxes
  {
  public:
    /**
     *  A box and its tag, as a tree is made of many at once.
     */
    struct Entry
    {
      std::uint64_t Tag    = 0;  /**< The tag. */
      Box           Bounds = {}; /**< The box. */
    };

    /**
     *  A tag near a point: which, and how far, as the caller measured it.
     */
    struct Close
    {
      std::uint64_t Tag      = 0;   /**< The tag.       */
      double        Distance = 0.0; /**< How far it is. */

      /**
       *  Whether two are equal exactly.
       *
       *  @param  lhs One.
       *  @param  rhs The other.
       *  @return     Whether they are equal.
       */
      [[nodiscard]] friend constexpr auto operator==(const Close &lhs, const Close &rhs) -> bool
      {
        return lhs.Tag == rhs.Tag && std::equal_to<double>{}(lhs.Distance, rhs.Distance);
      }
    };

    /**
     *  How far what a tag stands for lies from a point, as the caller measures it: no less than the
     *  distance from the point to the tag's box; a tag it gives no number for is passed by.
     */
    using Measure = std::function<double(std::uint64_t tag)>;

    /**
     *  What is told of each tag a box query finds: true to go on, false to stop.
     */
    using Visit = std::function<bool(std::uint64_t tag)>;

    /**
     *  What is told of each tag in order of distance: true to go on to the next, false to stop.
     */
    using Walk = std::function<bool(const Close &close)>;

    /**
     *  No limit on how far Nearest looks.
     */
    static constexpr double Anywhere = std::numeric_limits<double>::infinity();

    /**
     *  What Nearest looks for: how many tags at most, and how far at most.
     */
    struct Seeking
    {
      std::size_t Count  = 1;        /**< How many at most.                       */
      double      Within = Anywhere; /**< How far at most, Anywhere for no limit. */
    };

    /**
     *  No boxes.
     */
    GEOXY_API Boxes();

    /**
     *  Many boxes at once, packed in one pass.
     *
     *  @param entries The boxes and their tags; where a tag comes more than once, its last box.
     */
    GEOXY_API explicit Boxes(std::span<const Entry> entries);

    /**
     *  The same boxes, in a tree of their own.
     *
     *  @param other The boxes.
     */
    GEOXY_API Boxes(const Boxes &other);

    /**
     *  The boxes taken over.
     *
     *  @param other The boxes, left with none.
     */
    GEOXY_API Boxes(Boxes &&other) noexcept;

    /**
     *  The same boxes as another, in a tree of their own.
     *
     *  @param  other The boxes.
     *  @return       This.
     */
    GEOXY_API auto operator=(const Boxes &other) -> Boxes &;

    /**
     *  Another's boxes taken over.
     *
     *  @param  other The boxes, left with none.
     *  @return       This.
     */
    GEOXY_API auto operator=(Boxes &&other) noexcept -> Boxes &;

    /**
     *  The boxes let go.
     */
    GEOXY_API ~Boxes();

    /**
     *  How many boxes are kept.
     *
     *  @return The count.
     */
    [[nodiscard]] auto Size() const -> std::size_t { return m_Marks.size(); }

    /**
     *  Whether a tag's box is kept.
     *
     *  @param  tag The tag.
     *  @return     Whether it is.
     */
    [[nodiscard]] auto Holds(std::uint64_t tag) const -> bool { return m_Marks.contains(tag); }

    /**
     *  A tag's box, as it was put.
     *
     *  @param  tag The tag.
     *  @return     The box; none where the tag is not kept.
     */
    [[nodiscard]] GEOXY_API auto Bounds(std::uint64_t tag) const -> std::optional<Box>;

    /**
     *  A box inserted under a tag, or the tag's box moved to it where the tag is kept already.
     *
     *  @param tag The tag.
     *  @param box The box.
     */
    GEOXY_API void Insert(std::uint64_t tag, const Box &box);

    /**
     *  A tag's box taken out, the boxes above it made again on the way to the root.
     *
     *  @param tag The tag; nothing is done for one that is not kept.
     */
    GEOXY_API void Erase(std::uint64_t tag);

    /**
     *  The exact box of every box kept with finite sides, read from the root.
     *
     *  @return The box; empty where none is kept.
     */
    [[nodiscard]] GEOXY_API auto Extents() const -> Box;

    /**
     *  The tags whose boxes overlap a box, edges and corners included: what a crossing window picks
     *  from.
     *
     *  @param  box The box.
     *  @return     Their tags, in order.
     */
    [[nodiscard]] GEOXY_API auto Overlapping(const Box &box) const -> std::vector<std::uint64_t>;

    /**
     *  The tags whose boxes overlap a box, each told as it is found, in no given order.
     *
     *  @param box   The box.
     *  @param visit Told each tag: false to stop.
     */
    GEOXY_API void Overlapping(const Box &box, const Visit &visit) const;

    /**
     *  The tags whose boxes lie wholly within a box, exactly: what a window picks, for what lies
     *  within its own box.
     *
     *  @param  box The box, its edges included.
     *  @return     Their tags, in order.
     */
    [[nodiscard]] GEOXY_API auto Within(const Box &box) const -> std::vector<std::uint64_t>;

    /**
     *  The tags whose boxes lie wholly within a box, each told as it is found, in no given order.
     *
     *  @param box   The box, its edges included.
     *  @param visit Told each tag: false to stop.
     */
    GEOXY_API void Within(const Box &box, const Visit &visit) const;

    /**
     *  The tags nearest a point, as the caller measures them, nearest first.
     *
     *  The boxes tell which are near enough to measure, nearest first, so that only those that
     *  could be among the nearest are measured; a tag whose box is not finite is measured always.
     *
     *  @param  point    The point.
     *  @param  seeking  How many at most, and how far at most.
     *  @param  distance How far a tag lies from the point: never less than its box does.
     *  @return          The nearest, by distance, the lower tag first of those as near.
     */
    [[nodiscard]] GEOXY_API auto Nearest(Point point, const Seeking &seeking,
                                         const Measure &distance) const -> std::vector<Close>;

    /**
     *  Every tag in order of distance from a point, as the caller measures it, the lower tag first
     *  of those as near, told one at a time until the walk says stop: a pick that goes on until a
     *  tag passes an exact test, or until the distance passes an aperture, measures no more than it
     *  must.
     *
     *  @param point    The point.
     *  @param distance How far a tag lies from the point: never less than its box does.
     *  @param walk     Told each tag and its distance, nearest first: false to stop.
     */
    GEOXY_API void Nearest(Point point, const Measure &distance, const Walk &walk) const;

  private:
    std::unordered_map<std::uint64_t, std::size_t> m_Marks = {}; /**< Each tag's mark.  */
    std::vector<std::uint64_t> m_Tags  = {}; /**< Each mark's tag.                       */
    std::vector<Box>           m_Kept  = {}; /**< Each mark's box, exact.                */
    std::vector<std::size_t>   m_Spare = {}; /**< The marks given back, the last last.   */
    std::unique_ptr<Boughs>    m_Boughs;     /**< The boxes, in a tree.                  */
  };
} // namespace Geoxy

#endif /* GEOXY_BOXES_HPP */
