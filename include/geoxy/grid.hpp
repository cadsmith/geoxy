/***************************************************************************************************
 *  @file      geoxy/grid.hpp
 *  @brief     Curves prepared once and put in a uniform grid by their boxes, met all at once or
 *             with one curve more, only where their boxes come near.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_GRID_HPP
#define GEOXY_GRID_HPP

#include "geoxy/api.hpp"
#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/point.hpp"
#include "geoxy/prepared.hpp"
#include "geoxy/region.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/workers.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <span>
#include <utility>
#include <vector>

namespace Geoxy
{
  /**
   *  Curves prepared once, their boxes kept side by side, and a uniform grid of cells over them:
   *  what intersecting every curve with every other, or one more with all of them, asks for.
   *
   *  Intersecting every pair of n curves is n² / 2 intersections, nearly all of curves far apart.
   *  Here each curve's box, grown by as far as Intersect looks past it, goes in each cell of a
   *  uniform tiling it covers, and two curves are met only where their grown boxes overlap: every
   *  pair Intersect finds anything for, and few more. A curve with no box, a ray, an endless line,
   *  one of numbers that are not finite, is met with every other, as Intersect would meet it.
   *
   *  The grid holds its own curves, so that none moves under what is prepared from it.
   */
  class Grid
  {
  public:
    /**
     *  Two curves of a grid, by their indices, the first first.
     */
    using Pair = std::array<std::size_t, 2>;

    /**
     *  Where two curves meet, told as it is found: the pair, and its marking, which lasts only as
     *  long as the call.
     */
    using Sink = std::function<void(const Pair &curves, const Marking &marking)>;

    /**
     *  A curve prepared and a grid to meet it with, the curve first.
     */
    using With =
        std::pair<std::reference_wrapper<const Prepared>, std::reference_wrapper<const Grid>>;

    /**
     *  A grid and the pool to meet it on.
     */
    using By = std::pair<std::reference_wrapper<const Grid>, std::reference_wrapper<Workers>>;

    /**
     *  The index a curve that is not in the grid stands at, met with it.
     */
    static constexpr std::size_t Outside = std::numeric_limits<std::size_t>::max();

    /**
     *  A grid of curves, copied and prepared, and its tiling made for a tolerance.
     *
     *  @param curves    The curves, in the order their indices give.
     *  @param tolerance The tolerance it is met within the most: a finer one finds its pairs in the
     *                   same tiling, a coarser one in a tiling made for it.
     */
    GEOXY_API Grid(std::span<const Curve> curves, const Tolerance &tolerance);

    /**
     *  A grid of the same curves, prepared again.
     *
     *  @param other The grid.
     */
    GEOXY_API Grid(const Grid &other);

    /**
     *  A grid taken over, its curves where they were.
     *
     *  @param other The grid, left with no curve.
     */
    Grid(Grid &&other) noexcept = default;

    /**
     *  The same curves as another grid, prepared again.
     *
     *  @param  other The grid.
     *  @return       This.
     */
    GEOXY_API auto operator=(const Grid &other) -> Grid &;

    /**
     *  Another grid taken over, its curves where they were.
     *
     *  @param  other The grid, left with no curve.
     *  @return       This.
     */
    auto operator=(Grid &&other) noexcept -> Grid & = default;

    /**
     *  The curves let go.
     */
    ~Grid() = default;

    /**
     *  How many curves it holds.
     *
     *  @return The count.
     */
    [[nodiscard]] auto Size() const -> std::size_t { return m_Curves.size(); }

    /**
     *  A curve, prepared.
     *
     *  @pre          The index is below Size.
     *  @param  index Its index.
     *  @return       It, whose Shape is the curve.
     */
    [[nodiscard]] auto operator[](std::size_t index) const -> const Prepared &
    {
      return m_Prepared[index];
    }

    /**
     *  One curve moved, or changed: put in place of the one at an index, prepared, and boxed and
     *  put in the tiling again, the others left as they are.
     *
     *  @param index Its index, below Size: std::out_of_range is thrown past it.
     *  @param curve The curve now there.
     */
    GEOXY_API void Moved(std::size_t index, Curve curve);

    /**
     *  The curves that may lie in a box or cross it: each whose box, grown by as far as the grid is
     *  met within, overlaps it, and each with no box. No curve that reaches the box is left out,
     *  and few others are in; what Window and Crossing pick from.
     *
     *  @param  box The box.
     *  @return     Their indices, in order; those with no box alone for an empty box.
     */
    [[nodiscard]] GEOXY_API auto Near(const Box &box) const -> std::vector<std::size_t>;

    /**
     *  The curves a window picks: each lying wholly within a box, as Inside tells it.
     *
     *  @param  box The box, its edges included.
     *  @return     Their indices, in order.
     */
    [[nodiscard]] GEOXY_API auto Window(const Box &box) const -> std::vector<std::size_t>;

    /**
     *  The curves a crossing window picks: each that meets a box, within it or across its edges, as
     *  Intersects tells it.
     *
     *  @param  box The box, its edges included.
     *  @return     Their indices, in order.
     */
    [[nodiscard]] GEOXY_API auto Crossing(const Box &box) const -> std::vector<std::size_t>;

    /**
     *  The curves a polygon window picks, WPOLYGON: each lying wholly within a region, as Inside
     *  tells it within the grid's tolerance.
     *
     *  @param  region The region, a polygon's or any.
     *  @return        Their indices, in order.
     */
    [[nodiscard]] GEOXY_API auto Window(const Region &region) const -> std::vector<std::size_t>;

    /**
     *  The curves a polygon crossing window picks, CPOLYGON: each that meets a region, within it or
     *  across its edges, as Intersects tells it within the grid's tolerance.
     *
     *  @param  region The region, a polygon's or any.
     *  @return        Their indices, in order.
     */
    [[nodiscard]] GEOXY_API auto Crossing(const Region &region) const -> std::vector<std::size_t>;

  private:
    friend struct Gridwork;

    /**
     *  A curve in a cell: its box grown, and the cell its grown box starts in.
     */
    struct Entry
    {
      Box           Grown  = {}; /**< Its box, grown by as far as it is met within. */
      std::size_t   Index  = 0;  /**< The curve's index.                            */
      std::uint32_t Column = 0;  /**< The column its grown box starts in.           */
      std::uint32_t Row    = 0;  /**< The row its grown box starts in.              */
    };

    /**
     *  A cell's run of entries, and the room it has before the next.
     */
    struct Cell
    {
      std::size_t First = 0; /**< Where its entries start.          */
      std::size_t Count = 0; /**< How many it holds.                */
      std::size_t Room  = 0; /**< How many it may hold where it is. */
    };

    /**
     *  A uniform tiling over grown boxes.
     */
    struct Tiling
    {
      Point                     Origin  = {};  /**< The least corner of the first cell.  */
      double                    Scale   = 0.0; /**< Cells to a unit of length.           */
      std::uint32_t             Columns = 1;   /**< Cells along x.                       */
      std::uint32_t             Rows    = 1;   /**< Cells along y.                       */
      std::vector<Cell>         Cells   = {};  /**< The cells, row by row.               */
      std::vector<Entry>        Entries = {};  /**< Every cell's entries, and room.      */
      std::size_t               Live    = 0;   /**< How many entries the cells hold.     */
      std::vector<std::uint8_t> Loose   = {};  /**< Whether each curve has no grown box. */
      std::vector<std::size_t>  Unboxed = {};  /**< The curves with none, in order.      */
    };

    std::vector<Curve>        m_Curves    = {}; /**< The curves.                                */
    std::vector<Prepared>     m_Prepared  = {}; /**< Each prepared.                             */
    Tolerance                 m_Tolerance = {}; /**< The tolerance the tiling is made for.      */
    std::vector<double>       m_MinX      = {}; /**< Each box's least x, as Intersect reads it. */
    std::vector<double>       m_MinY      = {}; /**< Each box's least y.                        */
    std::vector<double>       m_MaxX      = {}; /**< Each box's greatest x.                     */
    std::vector<double>       m_MaxY      = {}; /**< Each box's greatest y.                     */
    std::vector<double>       m_Grain     = {}; /**< Each curve's grain, 16 steps at its size.  */
    std::vector<std::uint8_t> m_Shapes    = {}; /**< Each curve's shape, telling pairs' kinds.  */
    Tiling                    m_Tiling    = {}; /**< The tiling, made for the tolerance.        */
  };

  /**
   *  Where the curves of a grid meet, each with each: every pair that meets, as intersecting the
   *  two prepared answers it, to the bit; no other.
   *
   *  A pair meets where its marking holds a mark or a stretch, or says a mark may be missing. Each
   *  is met once, the lower index first. The pairs come by their kinds, a horizontal line and an
   *  upright one first, then two other lines, a line and an arc, two arcs, and the rest; within a
   *  kind, by the first index, then the second. So the order is the curves' own, however the grid
   *  was made and moved, and on however many threads it is met.
   *
   *  It is met on a pool made for the call, a thread for every 1,024 curves up to as many as the
   *  hardware runs at once, and on the calling thread alone under 2,048 curves. The sink is told on
   *  the calling thread alone, one pair after another, as intersecting the grid alone tells it.
   *
   *  @param grid      The grid.
   *  @param tolerance The tolerance, as Intersect takes it.
   *  @param sink      Told of each pair that meets.
   */
  GEOXY_API void IntersectAll(const Grid &grid, const Tolerance &tolerance, const Grid::Sink &sink);

  /**
   *  Where the curves of a grid meet, each with each, met on a pool: the same pairs, marks and
   *  order, to the bit, as the grid met alone gives.
   *
   *  The pairs are found and met on every thread of the pool; the sink is told on the calling
   *  thread alone, in order, while the others meet the pairs it tells next. A pool of one thread,
   *  or one working another call, meets the grid on the calling thread alone. What the sink throws
   *  is thrown on once the pool is done with the call.
   *
   *  @param by        The grid, and the pool to meet it on.
   *  @param tolerance The tolerance, as Intersect takes it.
   *  @param sink      Told of each pair that meets, on the calling thread.
   */
  GEOXY_API void IntersectAll(const Grid::By &by, const Tolerance &tolerance,
                              const Grid::Sink &sink);

  /**
   *  Where a curve meets the curves of a grid: every curve of it the curve meets, as intersecting
   *  the two prepared answers it, the curve first, to the bit; no other. The edges a trim or an
   *  extension stops at, the points a snap picks from.
   *
   *  Each pair is told as Grid::Outside and the grid's curve's index, in the order of that index,
   *  its marks along the curve.
   *
   *  @param with      The curve prepared, and the grid.
   *  @param tolerance The tolerance, as Intersect takes it: any, finer or coarser than the grid's.
   *  @param sink      Told of each curve it meets.
   */
  GEOXY_API void IntersectWith(const Grid::With &with, const Tolerance &tolerance,
                               const Grid::Sink &sink);
} // namespace Geoxy

#endif /* GEOXY_GRID_HPP */
