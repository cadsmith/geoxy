/***************************************************************************************************
 *  @file      geoxy/prepared.hpp
 *  @brief     Curves prepared to be met many times, and what two of them meet in, with each hit's
 *             places on them or with its point alone.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_PREPARED_HPP
#define GEOXY_PREPARED_HPP

#include "geoxy/api.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/point.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/vector.hpp"

#include <array>
#include <functional>
#include <memory>
#include <vector>

namespace Geoxy
{
  /**
   *  What a polyline or an arc prepared keeps to be met piece by piece: its legs and their grid.
   */
  struct Legwork;

  /**
   *  What stands for a polyline whose slightest arcs lie too far from their centers to be worked
   *  out through them: the polyline with those arcs made chords, itself prepared.
   */
  struct Standby;

  /**
   *  A curve prepared to be met many times: what belongs to it alone worked out once and kept
   *  beside it, its box, the size of its numbers, a line's exact way and its square, an arc's ends
   *  and its radius squared, a polyline's legs, each segment's box and circle, and their grid.
   *
   *  Intersect works these out for both curves each time it is asked. A drawing meets one curve
   *  with many, so a curve prepared once is met at a fraction of that: most pairs lie far apart,
   *  and are told so from two boxes already there.
   *
   *  The curve is not copied: it must stay where it is, unchanged, for as long as this is used.
   */
  class Prepared
  {
  public:
    /**
     *  Two prepared curves to meet, the first first.
     */
    using Pair = std::array<std::reference_wrapper<const Prepared>, 2>;

    /**
     *  A curve prepared.
     *
     *  @param curve The curve, which must outlive this.
     */
    GEOXY_API explicit Prepared(const Curve &curve);

    /**
     *  No curve that is about to go: a prepared curve holds none of its own.
     */
    explicit Prepared(const Curve &&) = delete;

    /**
     *  The curve.
     *
     *  @return It, where it was when prepared.
     */
    [[nodiscard]] auto Shape() const -> const Curve & { return *m_Shape; }

  private:
    friend struct Abutting;
    friend struct Glancing;
    friend struct Preparing;

    Box                  m_Box     = {};      /**< The box it cannot leave, where it has one.     */
    double               m_Reach   = 0.0;     /**< The largest magnitude of that box's edges.     */
    double               m_Size    = 0.0;     /**< The size of its numbers, its grain's measure.  */
    bool                 m_Boxed   = false;   /**< Whether it has a box, of finite edges.         */
    bool                 m_Plain   = false;   /**< Whether it has, and is a line or an arc.       */
    bool                 m_Round   = false;   /**< Whether a plain arc of a radius, up to a turn. */
    const Curve         *m_Shape   = nullptr; /**< The curve.                                     */
    std::array<Point, 2> m_Ends    = {};      /**< A line's or such an arc's start and end.       */
    double               m_Sweep   = 0.0;     /**< Such an arc's sweep, in radians.               */
    Vector               m_Slip    = {};      /**< What rounding left of a line's way.            */
    double               m_Long    = 0.0;     /**< That way's length, rounded.                    */
    double               m_Square  = 0.0;     /**< Its square, or an arc's radius squared.        */
    double               m_Under   = 0.0;     /**< What that square leaves: the two are it.       */
    double               m_Inverse = 0.0;     /**< 1 / a line's square, rounded.                  */
    std::shared_ptr<const Legwork> m_Legwork = nullptr; /**< A polyline's or an arc's legwork,
                                                             shared by its copies.             */
    std::shared_ptr<const Standby> m_Standby = nullptr; /**< What stands for a polyline with arcs
                                                             too slight to work out through their
                                                             centers, prepared too.          */
  };

  /**
   *  A place where two curves meet, and how: a hit without its places on the curves.
   */
  struct Mark
  {
    Point   At   = {};             /**< Where, as the hit's point is. */
    Contact Kind = Contact::Cross; /**< How they meet there.          */

    /**
     *  Whether two marks are equal exactly, point and kind.
     *
     *  @param  lhs One mark.
     *  @param  rhs The other.
     *  @return     Whether they are equal.
     */
    [[nodiscard]] friend constexpr auto operator==(const Mark &lhs, const Mark &rhs) -> bool
    {
      return lhs.At == rhs.At && lhs.Kind == rhs.Kind;
    }
  };

  /**
   *  All that intersecting two curves found, each hit kept as a mark.
   */
  struct Marking
  {
    std::vector<Mark>    Marks    = {}; /**< The hits, in order along the first curve.    */
    std::vector<Overlap> Overlaps = {}; /**< The stretches shared, along the first curve. */
    bool                 Spent = false; /**< Whether a hit may be missing, as Intersections::Spent
                                             tells it: never silent.                          */
  };

  /**
   *  Where two prepared curves meet, as Intersect answers for their curves, to the bit: the same
   *  hits, the same stretches, written into intersections kept by the caller.
   *
   *  Two that lie far apart are told so from their boxes, kept; all else is worked out as Intersect
   *  does.
   *
   *  @param curves    The two prepared curves, the first first.
   *  @param tolerance The tolerance, as Intersect takes it.
   *  @param into      The intersections, whatever they held replaced.
   */
  GEOXY_API void Intersect(const Prepared::Pair &curves, const Tolerance &tolerance,
                           Intersections &into);

  /**
   *  Where two prepared curves meet, each hit's point and kind alone: what Intersect answers for
   *  their curves, to the bit and in its order, with no place on either curve worked out.
   *
   *  Two lines or arcs that cross plainly, or plainly miss, are told from what is kept: every
   *  crossing far from each end and from where the curves would touch, its point the double nearest
   *  the true one. Each place Intersect would work out besides, a foot on a line, an angle on an
   *  arc and the point there, is most of its work, and an arc's ends by sine and cosine most of the
   *  rest: neither is done. All else is met as Intersect does, and its hits kept as marks.
   *
   *  @param curves    The two prepared curves, the first first.
   *  @param tolerance The tolerance, as Intersect takes it.
   *  @param into      The marking, whatever it held replaced.
   */
  GEOXY_API void Intersect(const Prepared::Pair &curves, const Tolerance &tolerance, Marking &into);
} // namespace Geoxy

#endif /* GEOXY_PREPARED_HPP */
