/***************************************************************************************************
 *  @file      source/belt/cut.hpp
 *  @brief     The curves of a band's edge cut wherever two of them meet, each point numbered once.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_BELT_CUT_HPP
#define GEOXY_BELT_CUT_HPP

#include "belt/carrier.hpp"
#include "belt/work.hpp"
#include "scratch.hpp"

#include "geoxy/box.hpp"
#include "geoxy/point.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>

namespace Geoxy::Belt
{
  /**
   *  A point's number, where a set of them is kept.
   */
  struct Marks
  {
    Work  *Curves = nullptr; /**< The curves and their cuts. */
    Nodes *Points = nullptr; /**< The points' numbers.       */
  };

  /**
   *  How an end of a curve is made: as told when it was drawn; a moved piece's, the shape's vertex
   *  there moved the distance square to the piece; any other's, its doubles themselves.
   *
   *  @param  work  The curves.
   *  @param  curve The curve.
   *  @param  ends  Whether its end, not its start.
   *  @return       The making.
   */
  auto Ending(const Work &work, std::size_t curve, bool ends) -> Recipe;

  /**
   *  A point worked out again by its making, in two doubles a number.
   *
   *  @param  work The curves.
   *  @param  by   Its making, known.
   *  @param  near The point as doubles made it.
   *  @return      The doubles nearest the point.
   */
  auto Remade(Work &work, const Recipe &by, Point near) -> Point;

  /**
   *  An end of a curve, as where a point comes from.
   *
   *  @param  curve The curve.
   *  @param  ends  Whether its end, not its start.
   *  @return       Where the point comes from.
   */
  inline auto Sourced(std::size_t curve, bool ends) -> Birth
  {
    return {.Curve = static_cast<std::uint32_t>(curve), .Other = Starts + (ends ? 1U : 0U)};
  }

  /**
   *  The number of a point, and where it comes from kept with it.
   *
   *  @param  marks The curves and the points.
   *  @param  point The point.
   *  @param  from  Where it comes from.
   *  @return       Its number.
   */
  inline auto Numbered(const Marks &marks, Point point, const Birth &from) -> std::size_t
  {
    Work             &work = *marks.Curves;
    const std::size_t node = marks.Points->Of(point);
    if (node == work.By.size())
    {
      work.By.push_back(from);
    }
    return node;
  }

  /**
   *  The points two lines or circles that come near touching share, as two doubles find them: but
   *  where those find the two to miss, and doubles found them within their numbers' rounding of
   *  touching, the one point they come nearest at.
   *
   *  @param  coarse The points as doubles found them.
   *  @param  fine   The points as two doubles find them.
   *  @return        The points.
   */
  inline auto Finer(const Hits &coarse, const Hits &fine) -> Hits
  {
    if (coarse.Graze && fine.Count == 0)
    {
      return {.At = {coarse.Touch}, .Count = 1};
    }
    return fine;
  }

  /**
   *  Whether two boxes lie apart by more than a reach.
   *
   *  @param  first  One box.
   *  @param  second The other.
   *  @param  reach  The reach.
   *  @return        Whether they do.
   */
  inline auto Parted(const Box &first, const Box &second, double reach) -> bool
  {
    /**
     *  The widest of the four gaps between their sides, told at once: most pairs tried are apart,
     *  and which side tells it is as good as chance.
     */
    return std::max({second.Min.X - first.Max.X, first.Min.X - second.Max.X,
                     second.Min.Y - first.Max.Y, first.Min.Y - second.Max.Y})
         > reach;
  }

  /**
   *  Every two curves whose boxes come together cut where they meet: each tried against every other
   *  where they are few, else each two that share a cell of a grid, once, in the cell that holds
   *  the low corner of what their boxes share.
   *
   *  @param marks The curves and the points.
   */
  void Met(const Marks &marks);
} // namespace Geoxy::Belt

#endif /* GEOXY_BELT_CUT_HPP */
