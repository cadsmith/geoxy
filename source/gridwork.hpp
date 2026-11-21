/***************************************************************************************************
 *  @file      source/gridwork.hpp
 *  @brief     What is worked out of a grid: its curves' boxes grown, its tiling, and the pairs the
 *             tiling gives, sorted by kind; for the tests and benchmarks that look inside.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_GRIDWORK_HPP
#define GEOXY_GRIDWORK_HPP

#include "geoxy/box.hpp"
#include "geoxy/grid.hpp"
#include "geoxy/prepared.hpp"
#include "geoxy/tolerance.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

namespace Geoxy
{
  /**
   *  What is worked out of a grid's curves, its tiling, and the pairs it gives.
   */
  struct Gridwork
  {
    /**
     *  The kinds pairs are met by, in the order they are met: one kernel each, where one is.
     */
    enum class Kind : std::uint8_t
    {
      Square, /**< A horizontal line and an upright one: they cross at (x of one, y of other). */
      Lines,  /**< Two other lines.                                                            */
      Lined,  /**< A line and an arc.                                                          */
      Rounds, /**< Two arcs.                                                                   */
      Other,  /**< Any other pair.                                                             */
    };

    /**
     *  How many kinds there are.
     */
    static constexpr std::size_t Kinds = 5;

    /**
     *  The pairs to meet, a list for each kind, each in order of its first index, then its second.
     */
    using Pairs = std::array<std::vector<Grid::Pair>, Kinds>;

    /**
     *  A run of curves' indices: the first, and the one past the last.
     */
    using Span = std::array<std::size_t, 2>;

    /**
     *  How many spans of curves a pool's thread has to work, each, so that one slower than the
     *  others, or woken late, holds the rest up little.
     */
    static constexpr std::size_t Spans = 8;

    /**
     *  How many curves a span holds at the least, a power of two: some 470 pairs of a floor plan,
     *  16 µs of meeting, so that what each kind of a span costs to take and tell is a few percent.
     */
    static constexpr std::size_t Least = 256;

    /**
     *  How many curves of a grid it takes, at the least, for a thread more to pay for its start
     *  where IntersectAll is given no pool: about a tenth of a millisecond of meeting each.
     */
    static constexpr std::size_t Share = 1024;

    /**
     *  How many threads a pool has, at the least, for the calling thread to leave meeting to the
     *  others and tell what they meet as they go: with one or two others, it meets too, and tells
     *  once all is met, as they meet slower than it tells.
     */
    static constexpr std::size_t Helping = 4;

    /**
     *  A tiling over a grid's curves, and the tolerance their boxes in it are grown at.
     */
    struct Search
    {
      const Grid         *Of    = nullptr; /**< The grid.                                        */
      const Grid::Tiling *Tiles = nullptr; /**< Its tiling, or one made for a coarser tolerance. */
      Tolerance           Grown = {};      /**< The tolerance the tiling's boxes are grown at.   */
    };

    /**
     *  A run of pairs of one kind, and the tolerance they are met within.
     */
    struct Run
    {
      Kind                        Of     = Kind::Other; /**< Their kind.                      */
      std::span<const Grid::Pair> Pairs  = {};          /**< The pairs, in the order told.    */
      Tolerance                   Within = {}; /**< The tolerance, as Intersect takes it. */
    };

    /**
     *  The cells a box covers: the columns and rows from where it starts to where it ends.
     */
    struct Area
    {
      std::uint32_t Left  = 0; /**< The first column. */
      std::uint32_t Low   = 0; /**< The first row.    */
      std::uint32_t Right = 0; /**< The last column.  */
      std::uint32_t High  = 0; /**< The last row.     */
    };

    /**
     *  A curve's box and size read, as meeting it reads them, and its shape.
     *
     *  @param grid  The grid.
     *  @param index The curve's index.
     */
    static void Measured(Grid &grid, std::size_t index);

    /**
     *  A curve's box grown by as far as Intersect looks past it, at a tolerance, and a step of a
     *  double more on every side.
     *
     *  @param  grid      The grid.
     *  @param  index     The curve's index.
     *  @param  tolerance The tolerance.
     *  @return           The box grown; nothing where the curve has no box, or its box grown runs
     *                    past every double.
     */
    static auto Grown(const Grid &grid, std::size_t index, const Tolerance &tolerance)
        -> std::optional<Box>;

    /**
     *  A curve's box grown, for a curve that is not in the grid, at a tolerance that may be coarser
     *  than the grid's: grown by as much more, so that it reaches every box of the grid's tiling it
     *  would at its own.
     *
     *  @param  grid      The grid.
     *  @param  curve     The curve, made ready.
     *  @param  tolerance The tolerance.
     *  @return           The box grown; nothing as Grown gives none.
     */
    static auto Probed(const Grid &grid, const Prepared &curve, const Tolerance &tolerance)
        -> std::optional<Box>;

    /**
     *  A tiling over a grid's curves, their boxes grown at a tolerance.
     *
     *  @param  grid      The grid.
     *  @param  tolerance The tolerance.
     *  @return           The tiling.
     */
    static auto Tiled(const Grid &grid, const Tolerance &tolerance) -> Grid::Tiling;

    /**
     *  The cells of a tiling a box covers.
     *
     *  @param  tiling The tiling.
     *  @param  box    The box.
     *  @return        Its cells.
     */
    static auto Covered(const Grid::Tiling &tiling, const Box &box) -> Area;

    /**
     *  Every pair of a grid's curves that may meet within a tolerance: their boxes grown overlap,
     *  or either has none.
     *
     *  @param  grid      The grid.
     *  @param  tolerance The tolerance.
     *  @return           The pairs, by kind.
     */
    static auto Candidates(const Grid &grid, const Tolerance &tolerance) -> Pairs;

    /**
     *  Every curve of a grid a curve that is not in it may meet within a tolerance.
     *
     *  @param  grid      The grid.
     *  @param  curve     The curve, made ready.
     *  @param  tolerance The tolerance.
     *  @return           Their indices, in order.
     */
    static auto Near(const Grid &grid, const Prepared &curve, const Tolerance &tolerance)
        -> std::vector<std::size_t>;

    /**
     *  Every curve of a grid whose box grown, in its tiling, overlaps a box, and every curve with
     *  none.
     *
     *  @param  grid  The grid.
     *  @param  probe The box; every curve is near none.
     *  @return       Their indices, in order.
     */
    static auto Nearby(const Grid &grid, const std::optional<Box> &probe)
        -> std::vector<std::size_t>;

    /**
     *  A run of pairs of one kind met, each that meets told to a sink, in the run's order: as
     *  meeting the two made ready answers it, to the bit.
     *
     *  Two lines, horizontal and upright or any other, are met in a loop of their own by their
     *  exact pass alone where it tells them, as Plainly does: no meeting worked out, and no marking
     *  made for a pair that meets in nothing. Any other pair, and two lines the pass leaves, are
     *  met as Intersect meets them made ready.
     *
     *  @param grid The grid.
     *  @param run  The pairs, their kind, and the tolerance.
     *  @param sink Told of each pair that meets.
     */
    static void Met(const Grid &grid, const Run &run, const Grid::Sink &sink);

    /**
     *  Every pair of a grid that meets, met on the calling thread alone and told to a sink in
     *  IntersectAll's order.
     *
     *  @param grid      The grid.
     *  @param tolerance The tolerance, as Intersect takes it.
     *  @param sink      Told of each pair that meets.
     */
    static void Alone(const Grid &grid, const Tolerance &tolerance, const Grid::Sink &sink);

    /**
     *  Every pair of a grid that meets, met on a pool and told to a sink in IntersectAll's order,
     *  on the calling thread.
     *
     *  The curves are cut into spans of their indices, about eight for each thread. Every thread
     *  seeks the pairs of a span after another until none is left; then the threads but the caller
     *  meet the pairs of each kind of each span, kind after kind, each into a list of its own,
     *  while the caller tells those lists to the sink in that order as each is done. A kind of a
     *  span waits only for its span's pairs, and the caller for the list it tells next. In a pool
     *  of fewer than Helping threads the caller meets too, then tells.
     *
     *  @pre              The pool has a thread besides the caller.
     *  @param  by        The grid, and the pool.
     *  @param  tolerance The tolerance, as Intersect takes it.
     *  @param  sink      Told of each pair that meets.
     *  @return           Whether the pool was free: nothing is told where it was not.
     */
    static auto Together(const Grid::By &by, const Tolerance &tolerance, const Grid::Sink &sink)
        -> bool;

    /**
     *  A curve's entries taken out of the grid's tiling, as its box was.
     *
     *  @param grid  The grid.
     *  @param index The curve's index.
     */
    static void Ungridded(Grid &grid, std::size_t index);

    /**
     *  A curve's entries put in the grid's tiling, as its box now is.
     *
     *  @param grid  The grid.
     *  @param index The curve's index.
     */
    static void Regridded(Grid &grid, std::size_t index);

    /**
     *  The pairs whose boxes grown share a cell, among a run of the tiling's cells, each told in
     *  the one cell where both begin.
     *
     *  @param tiling The tiling.
     *  @param cells  The cells' places, row by row, from the first to the one past the last.
     *  @param into   Where the pairs go, the lower index first, in no order.
     */
    static void Cellwise(const Grid::Tiling &tiling, const Span &cells,
                         std::vector<Grid::Pair> &into);

    /**
     *  The pairs of which a curve with no box grown is one, it with every other curve, whose lower
     *  index lies in a span.
     *
     *  @param tiling The tiling, which tells which curves have none.
     *  @param span   The lower indices.
     *  @param into   Where the pairs go, the lower index first, in no order.
     */
    static void Loosely(const Grid::Tiling &tiling, const Span &span,
                        std::vector<Grid::Pair> &into);

    /**
     *  Pairs put in order and with their kinds.
     *
     *  @param  grid  The grid.
     *  @param  found The pairs, the lower index first, in no order.
     *  @param  span  The lower indices they have.
     *  @return       The pairs, by kind, each in order of its first index, then its second.
     */
    static auto Sorted(const Grid &grid, std::span<const Grid::Pair> found, const Span &span)
        -> Pairs;

  private:
    /**
     *  An entry taken out of a cell, the last put where it was.
     *
     *  @param tiling The tiling.
     *  @param place  The cell's place, row by row.
     *  @param index  The curve's index.
     */
    static void Removed(Grid::Tiling &tiling, std::size_t place, std::size_t index);

    /**
     *  An entry put in a cell, which moves to the end of the entries with twice its room where it
     *  has none left.
     *
     *  @param tiling The tiling.
     *  @param place  The cell's place, row by row.
     *  @param entry  The entry.
     */
    static void Inserted(Grid::Tiling &tiling, std::size_t place, const Grid::Entry &entry);

    /**
     *  Every cell's entries put one run after another again, the room left by cells that moved
     *  given back.
     *
     *  @param tiling The tiling.
     */
    static void Compacted(Grid::Tiling &tiling);
  };
} // namespace Geoxy

#endif /* GEOXY_GRIDWORK_HPP */
