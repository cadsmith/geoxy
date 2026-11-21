/***************************************************************************************************
 *  @file      source/grid.cpp
 *  @brief     A grid's curves made ready once and put in a uniform tiling by their boxes grown by
 *             as far as Intersect looks, so that only pairs Intersect may find anything for are
 *             met.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/grid.hpp"

#include "gridwork.hpp"
#include "plain.hpp"
#include "pooled.hpp"
#include "real.hpp"
#include "sane.hpp"
#include "scratch.hpp"
#include "thorough.hpp"

#include "geoxy/arc.hpp"
#include "geoxy/box.hpp"
#include "geoxy/contain.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/extents.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/line.hpp"
#include "geoxy/loop.hpp"
#include "geoxy/prepared.hpp"
#include "geoxy/region.hpp"
#include "geoxy/tolerance.hpp"
#include "geoxy/workers.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <numeric>
#include <optional>
#include <span>
#include <thread>
#include <utility>
#include <variant>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  The box of a region's edges, where a curve it picks lies.
     *
     *  @param  region The region.
     *  @return        The box of every edge's extents.
     */
    auto Encased(const Region &region) -> Box
    {
      Box box = {};
      for (const Loop &loop : region.Loops)
      {
        for (const Curve &edge : loop.Edges)
        {
          box = Joined(box, Extents(edge));
        }
      }
      return box;
    }

    /**
     *  The shapes a pair's kind is told from.
     */
    enum class Shape : std::uint8_t
    {
      Level,   /**< A horizontal line.                           */
      Upright, /**< An upright line.                             */
      Slanted, /**< Any other line, one of no length among them. */
      Round,   /**< An arc.                                      */
      Other,   /**< Any other curve.                             */
    };

    /**
     *  The kind of a pair of each two shapes.
     */
    constexpr std::array<std::array<Gridwork::Kind, Gridwork::Kinds>, Gridwork::Kinds> PairKinds = {
        {
            {Gridwork::Kind::Lines, Gridwork::Kind::Square, Gridwork::Kind::Lines,
             Gridwork::Kind::Lined, Gridwork::Kind::Other},
            {Gridwork::Kind::Square, Gridwork::Kind::Lines, Gridwork::Kind::Lines,
             Gridwork::Kind::Lined, Gridwork::Kind::Other},
            {Gridwork::Kind::Lines, Gridwork::Kind::Lines, Gridwork::Kind::Lines,
             Gridwork::Kind::Lined, Gridwork::Kind::Other},
            {Gridwork::Kind::Lined, Gridwork::Kind::Lined, Gridwork::Kind::Lined,
             Gridwork::Kind::Rounds, Gridwork::Kind::Other},
            {Gridwork::Kind::Other, Gridwork::Kind::Other, Gridwork::Kind::Other,
             Gridwork::Kind::Other, Gridwork::Kind::Other},
        }};

    /**
     *  A curve's shape.
     *
     *  @param  curve The curve.
     *  @return       Its shape: a line horizontal or upright where its ends share a coordinate
     *                exactly, and not both.
     */
    auto Shaped(const Curve &curve) -> Shape
    {
      if (const auto *line = std::get_if<Line>(&curve))
      {
        const bool level   = Exactly(line->Start.Y, line->End.Y);
        const bool upright = Exactly(line->Start.X, line->End.X);
        if (level == upright)
        {
          return Shape::Slanted;
        }
        return level ? Shape::Level : Shape::Upright;
      }
      return std::holds_alternative<Arc>(curve) ? Shape::Round : Shape::Other;
    }

    /**
     *  A double a step or more below another, so that any number the other was rounded from lies
     *  above it: the step up from its magnitude is no smaller than the step down from it.
     *
     *  @param  value The double.
     *  @return       The double below; not a number for one that is not finite.
     */
    auto Undershot(double value) -> double
    {
      return value - Above(std::abs(value));
    }

    /**
     *  A double a step or more above another, so that any number the other was rounded from lies
     *  below it.
     *
     *  @param  value The double.
     *  @return       The double above; not a number for one that is not finite.
     */
    auto Overshot(double value) -> double
    {
      return value + Above(std::abs(value));
    }

    /**
     *  A box grown by a margin on every side, and a step of a double more, so that it holds the box
     *  grown exactly.
     *
     *  @param  box    The box.
     *  @param  margin The margin.
     *  @return        The box grown; nothing where an edge of it is not finite.
     */
    auto Swollen(const Box &box, double margin) -> std::optional<Box>
    {
      const Box grown = {
          .Min = {.X = Undershot(box.Min.X - margin), .Y = Undershot(box.Min.Y - margin)},
          .Max = {.X = Overshot(box.Max.X + margin), .Y = Overshot(box.Max.Y + margin)},
      };
      /**
       *  A quarter of each edge, so that the sum of four finite ones is finite.
       */
      const double sum =
          (grown.Min.X / 4) + (grown.Min.Y / 4) + (grown.Max.X / 4) + (grown.Max.Y / 4);
      return std::isfinite(sum) ? std::optional<Box>{grown} : std::nullopt;
    }

    /**
     *  The largest magnitude of a box's edges.
     *
     *  @param  box The box.
     *  @return     The magnitude.
     */
    auto Outmost(const Box &box) -> double
    {
      return std::max(
          {std::abs(box.Min.X), std::abs(box.Max.X), std::abs(box.Min.Y), std::abs(box.Max.Y)});
    }

    /**
     *  The least length a curve is met within: 16 steps of a double at the larger of its box's
     *  reach and its numbers' size, as Intersect takes one or the other.
     *
     *  @param  curve The curve.
     *  @param  box   Its box, as Intersect reads it.
     *  @return       The length.
     */
    auto Fineness(const Curve &curve, const Box &box) -> double
    {
      return 16.0 * Above(std::max(Outmost(box), Sized(curve)));
    }

    /**
     *  Where along a tiling's columns or rows a place falls: its whole part, from the first to the
     *  last. It never decreases as the place grows, so that a box's cells hold every place of it.
     *
     *  @param  place The place, in cells from the tiling's origin.
     *  @param  count How many columns or rows there are.
     *  @return       The column or row; the first for a place that is not a number.
     */
    auto Binned(double place, std::uint32_t count) -> std::uint32_t
    {
      const auto last = static_cast<double>(count - 1);
      return static_cast<std::uint32_t>(std::min(last, std::max(0.0, place)));
    }

    /**
     *  How many columns or rows of cells span a length, no more than some.
     *
     *  @param  cells How many cells long it is.
     *  @param  most  The most there may be.
     *  @return       The count.
     */
    auto Rowed(double cells, double most) -> std::uint32_t
    {
      return static_cast<std::uint32_t>(std::min(most, std::floor(cells) + 1.0));
    }

    /**
     *  Whether meeting two curves found anything: a mark, a stretch, or that a mark may be missing.
     *
     *  @param  marking The marking.
     *  @return         Whether it did.
     */
    auto Met(const Marking &marking) -> bool
    {
      return marking.Marks.size()
               + marking.Overlaps.size()
               + static_cast<std::size_t>(marking.Spent)
           > 0;
    }

    /**
     *  Pairs of curves met one after another, each that meets told to a sink, one marking kept for
     *  them all.
     */
    class Teller
    {
    public:
      /**
       *  A teller of what pairs meet in, within a tolerance.
       *
       *  @param tolerance The tolerance, as Intersect takes it.
       *  @param sink      Told of each pair that meets.
       */
      Teller(const Tolerance &tolerance, const Grid::Sink &sink)
          : m_Tolerance(tolerance), m_Sane(Sane(tolerance)), m_Sink(&sink)
      {
      }

      /**
       *  Two curves met as Intersect meets them made ready, and told where they meet.
       *
       *  @param curves The two ready curves.
       *  @param pair   Their indices, as told.
       */
      void Paired(const Prepared::Pair &curves, const Grid::Pair &pair)
      {
        Intersect(curves, m_Tolerance, m_Marking);
        Told(pair);
      }

      /**
       *  Two lines met as Plainly meets them: a hit their exact passes tell as a mark of its point
       *  and kind, nothing told with no marking made, and anything else as Intersect meets them.
       *
       *  @param lines The two ready lines.
       *  @param pair  Their indices, as told.
       */
      void Lines(const Prepared::Pair &lines, const Grid::Pair &pair)
      {
        const Plain plain = Plainly(lines, m_Sane, m_Marking);
        if (!plain.Told)
        {
          Told(pair);
          return;
        }
        if (plain.Count == 0)
        {
          return;
        }
        m_Marking.Marks.assign(1, {.At = plain.At, .Kind = plain.Kind});
        m_Marking.Overlaps.clear();
        m_Marking.Spent = false;
        (*m_Sink)(pair, m_Marking);
      }

    private:
      /**
       *  A pair told to the sink, where what it met in, in the marking, is anything.
       *
       *  @param pair The pair's indices, as told.
       */
      void Told(const Grid::Pair &pair)
      {
        if (Met(m_Marking))
        {
          (*m_Sink)(pair, m_Marking);
        }
      }

      Tolerance         m_Tolerance = {};      /**< The tolerance, as Intersect takes it. */
      Tolerance         m_Sane      = {};      /**< The tolerance, as Sane makes it.      */
      const Grid::Sink *m_Sink      = nullptr; /**< Told of each pair that meets.         */
      Marking           m_Marking   = {};      /**< What the last pair met in.            */
    };

    /**
     *  A pair that met, kept to be told later: its marks and stretches lie in its list's, after
     *  those of the pairs kept before it.
     */
    struct Report
    {
      Grid::Pair  Curves   = {};    /**< The pair.                      */
      std::size_t Marks    = 0;     /**< How many marks it met in.      */
      std::size_t Overlaps = 0;     /**< How many stretches.            */
      bool        Spent    = false; /**< Whether a mark may be missing. */
    };

    /**
     *  The pairs of a run that met, kept in its order with what each met in, to be told to a sink
     *  on another thread.
     */
    struct Catch
    {
      std::vector<Report>  Pairs    = {}; /**< The pairs, in order.            */
      std::vector<Mark>    Marks    = {}; /**< Their marks, one after another. */
      std::vector<Overlap> Overlaps = {}; /**< Their stretches, likewise.      */
    };

    /**
     *  A sink that keeps what it is told.
     *
     *  @param  kept Where.
     *  @return      The sink.
     */
    auto Keeping(Catch &kept) -> Grid::Sink
    {
      return [&kept](const Grid::Pair &curves, const Marking &marking) -> void
      {
        kept.Pairs.push_back({.Curves   = curves,
                              .Marks    = marking.Marks.size(),
                              .Overlaps = marking.Overlaps.size(),
                              .Spent    = marking.Spent});
        kept.Marks.insert(kept.Marks.end(), marking.Marks.begin(), marking.Marks.end());
        kept.Overlaps.insert(kept.Overlaps.end(), marking.Overlaps.begin(), marking.Overlaps.end());
      };
    }

    /**
     *  What a run of pairs met in, as it was kept, told to a sink in its order, each pair's marks
     *  and stretches put in one marking.
     *
     *  @param kept    What was kept.
     *  @param sink    The sink.
     *  @param marking The marking.
     */
    void Retold(const Catch &kept, const Grid::Sink &sink, Marking &marking)
    {
      const std::span<const Mark>    marks    = kept.Marks;
      const std::span<const Overlap> overlaps = kept.Overlaps;
      std::size_t                    mark     = 0;
      std::size_t                    overlap  = 0;
      for (const Report &told : kept.Pairs)
      {
        const std::span<const Mark>    its     = marks.subspan(mark, told.Marks);
        const std::span<const Overlap> shared  = overlaps.subspan(overlap, told.Overlaps);
        mark                                  += told.Marks;
        overlap                               += told.Overlaps;
        marking.Marks.assign(its.begin(), its.end());
        marking.Overlaps.assign(shared.begin(), shared.end());
        marking.Spent = told.Spent;
        sink(told.Curves, marking);
      }
    }

    /**
     *  The next of some work to take, from a count that every thread taking it moves on.
     *
     *  @param  taken The count.
     *  @return       What is taken; the count of all there is, or more, where none is left.
     */
    auto Taken(std::atomic<std::size_t> &taken) -> std::size_t
    {
      return taken.fetch_add(1, std::memory_order_relaxed);
    }

    /**
     *  The pairs a run of a tiling's cells found, put in order of the span of curves their lower
     *  index lies in, and where each span's start.
     */
    struct Bucket
    {
      std::vector<Grid::Pair>  Pairs  = {}; /**< The pairs, span after span. */
      std::vector<std::size_t> Starts = {}; /**< Where each span's start, and one past the
                                                 last span's end.                            */
    };

    /**
     *  One grid met on a pool, as Together tells it: its tiling's cells cut into runs, each run's
     *  pairs found and put in order of the span of curves their lower index lies in; each span's
     *  pairs gathered from every run and sorted by kind; each kind of each span met and kept; and
     *  what was kept told in turn, kind after kind, span after span.
     */
    class Shift
    {
    public:
      /**
       *  The work of meeting a grid on some threads, cut up.
       *
       *  A span is a power of two of curves long, so that the span a curve lies in is its index
       *  shifted: as near as that allows to eight spans for each thread, and no fewer curves than
       *  pay for a span's lists.
       *
       *  @param search    The grid, its tiling and the tolerance the tiling's boxes are grown at.
       *  @param tolerance The tolerance, as Intersect takes it.
       *  @param threads   How many threads work it.
       */
      Shift(const Gridwork::Search &search, const Tolerance &tolerance, std::size_t threads)
          : m_Search(search), m_Tolerance(tolerance),
            m_Shift(static_cast<std::size_t>(std::countr_zero(std::max(
                Gridwork::Least, std::bit_ceil(search.Of->Size() / (Gridwork::Spans * threads)))))),
            m_Count(search.Of->Size()),
            m_Spans((m_Count + (std::size_t{1} << m_Shift) - 1) >> m_Shift),
            m_Cuts(Gridwork::Spans * threads), m_Cut(Gridwork::Spans * threads), m_Found(m_Spans),
            m_Sought(m_Spans), m_Kept(m_Spans * Gridwork::Kinds), m_Met(m_Spans * Gridwork::Kinds)
      {
      }

      /**
       *  The runs of cells' pairs found and put in order of their spans, one run after another,
       *  until every run is taken.
       */
      void Find()
      {
        const auto       &tiling = *m_Search.Tiles;
        const std::size_t cells  = tiling.Cells.size();
        const std::size_t runs   = m_Cuts.size();
        for (std::size_t run = Taken(m_Finding); run < runs; run = Taken(m_Finding))
        {
          const Pooled::Done      done{m_Cut[run]};
          std::vector<Grid::Pair> found;
          Gridwork::Cellwise(tiling, {run * cells / runs, (run + 1) * cells / runs}, found);
          Bucket &cut = m_Cuts[run];
          cut.Starts.assign(m_Spans + 1, 0);
          for (const Grid::Pair &pair : found)
          {
            ++cut.Starts[(pair[0] >> m_Shift) + 1];
          }
          std::partial_sum(cut.Starts.begin(), cut.Starts.end(), cut.Starts.begin());
          List<std::size_t> next(cut.Starts.begin(), cut.Starts.end());
          cut.Pairs.resize(found.size());
          for (const Grid::Pair &pair : found)
          {
            cut.Pairs[next[pair[0] >> m_Shift]] = pair;
            ++next[pair[0] >> m_Shift];
          }
        }
      }

      /**
       *  Each span's pairs gathered from every run once all are found, those of curves with no box
       *  added, and sorted by kind, one span after another, until every span is taken.
       */
      void Seek()
      {
        for (std::size_t span = Taken(m_Seeking); span < m_Spans; span = Taken(m_Seeking))
        {
          const Pooled::Done   done{m_Sought[span]};
          const Gridwork::Span curves = {span << m_Shift, std::min(m_Count, (span + 1) << m_Shift)};
          std::vector<Grid::Pair> found;
          for (std::size_t run = 0; run < m_Cuts.size(); ++run)
          {
            m_Cut[run].wait(false, std::memory_order_acquire);
            const Bucket                     &cut   = m_Cuts[run];
            const std::span<const Grid::Pair> pairs = std::span{cut.Pairs}.subspan(
                cut.Starts[span], cut.Starts[span + 1] - cut.Starts[span]);
            found.insert(found.end(), pairs.begin(), pairs.end());
          }
          Gridwork::Loosely(*m_Search.Tiles, curves, found);
          m_Found[span] = Gridwork::Sorted(*m_Search.Of, found, curves);
        }
      }

      /**
       *  The kinds of the spans met, kind after kind and span after span, each once its span's
       *  pairs are sought, until every one is taken.
       */
      void Intersect()
      {
        for (std::size_t run = Taken(m_Meeting); run < m_Met.size(); run = Taken(m_Meeting))
        {
          const Pooled::Done done{m_Met[run]};
          const std::size_t  kind = run / m_Spans;
          const std::size_t  span = run % m_Spans;
          m_Sought[span].wait(false, std::memory_order_acquire);
          const std::vector<Grid::Pair> &pairs = m_Found[span].at(kind);
          Catch                         &kept  = m_Kept[run];
          kept.Pairs.reserve(pairs.size());
          kept.Marks.reserve(pairs.size());
          Gridwork::Met(
              *m_Search.Of,
              {.Of = static_cast<Gridwork::Kind>(kind), .Pairs = pairs, .Within = m_Tolerance},
              Keeping(kept));
        }
      }

      /**
       *  What each kind of each span met in, told to a sink in turn, each once it is met.
       *
       *  @param sink The sink.
       */
      void Tell(const Grid::Sink &sink)
      {
        Marking marking;
        for (std::size_t run = 0; run < m_Met.size(); ++run)
        {
          m_Met[run].wait(false, std::memory_order_acquire);
          Retold(m_Kept[run], sink, marking);
          m_Kept[run] = {};
        }
      }

    private:
      Gridwork::Search               m_Search    = {}; /**< The grid and its tiling.              */
      Tolerance                      m_Tolerance = {}; /**< The tolerance, as Intersect takes it. */
      std::size_t                    m_Shift     = 0;  /**< A span's length, as a power of 2.     */
      std::size_t                    m_Count     = 0;  /**< How many curves.                      */
      std::size_t                    m_Spans     = 0;  /**< How many spans.                       */
      std::vector<Bucket>            m_Cuts;           /**< Each run of cells' pairs.             */
      std::vector<std::atomic<bool>> m_Cut;            /**< Whether each is found.                */
      std::vector<Gridwork::Pairs>   m_Found;          /**< Each span's pairs, by kind.           */
      std::vector<std::atomic<bool>> m_Sought;         /**< Whether each span's are.              */
      std::vector<Catch>             m_Kept;           /**< What each kind of each span met
                                                           in, kind after kind.               */
      std::vector<std::atomic<bool>> m_Met;            /**< Whether each is met.              */
      std::atomic<std::size_t>       m_Finding = 0;    /**< The next run of cells to find in. */
      std::atomic<std::size_t>       m_Seeking = 0;    /**< The next span to seek.            */
      std::atomic<std::size_t>       m_Meeting = 0;    /**< The next kind of a span to meet.  */
    };
  } // namespace

  void Gridwork::Measured(Grid &grid, std::size_t index)
  {
    const Curve &curve = grid.m_Curves[index];
    const Box    box   = Confines(curve).value_or(Box{});
    /**
     *  A curve with no box keeps one of edges not finite, which grown is none.
     */
    grid.m_MinX[index]   = box.Min.X;
    grid.m_MinY[index]   = box.Min.Y;
    grid.m_MaxX[index]   = box.Max.X;
    grid.m_MaxY[index]   = box.Max.Y;
    grid.m_Grain[index]  = Fineness(curve, box);
    grid.m_Shapes[index] = static_cast<std::uint8_t>(Shaped(curve));
  }

  /**
   *  Intersect tells two curves far apart from their boxes, as Intersect reads them, before
   *  anything else: far where along x or along y the gap between them is more than a clearance,
   *  four times the length it meets them within for two lines or arcs, sixteen for any other pair.
   *  That length is the tolerance's, or 16 steps of a double at the larger of the two curves' sizes
   *  where that is more: their boxes' reach for two lines or arcs, their numbers' size else. So the
   *  clearance of a pair is no more than 16 max(L1, L2), where each curve's L is the larger of the
   *  tolerance and its grain, 16 steps at the larger of its reach and size, as the step grows with
   *  the magnitude.
   *
   *  Each box is grown here by 32 L, and a step more on each side, so that it holds the box grown
   *  exactly. Two boxes grown that do not overlap have, along x or y, an exact gap g of more than
   *  32 L1 + 32 L2, at least twice the clearance c. The gap Intersect works out is g rounded, which
   *  is no less than 2c, a double, and so more than c where c is more than nothing; where c is
   *  nothing, g is more than nothing, and so is the difference of two doubles that are not equal.
   *  So the pair is far, as Intersect tells it, and meets in nothing: no pair Intersect finds
   *  anything for is left out. A box grown whose edges are not all finite, a tolerance that is not
   *  a number, a curve with no box, is none, and its curve is met with every other.
   */
  auto Gridwork::Grown(const Grid &grid, std::size_t index, const Tolerance &tolerance)
      -> std::optional<Box>
  {
    const Box box = {
        .Min = {.X = grid.m_MinX[index], .Y = grid.m_MinY[index]},
        .Max = {.X = grid.m_MaxX[index], .Y = grid.m_MaxY[index]},
    };
    return Swollen(box, 32.0 * std::max(tolerance.Length, grid.m_Grain[index]));
  }

  /**
   *  Each box of the grid's tiling is grown for the grid's tolerance. Met within one that is finer,
   *  a pair's clearance is no more than it was; within one that is coarser, by Δ, each L is no more
   *  than the grid's L plus Δ, so the curve's box grown by 32 Δ more reaches every box its own
   *  would.
   */
  auto Gridwork::Probed(const Grid &grid, const Prepared &curve, const Tolerance &tolerance)
      -> std::optional<Box>
  {
    const Curve &shape  = curve.Shape();
    const Box    box    = Confines(shape).value_or(Box{});
    const double coarse = std::max(0.0, tolerance.Length - grid.m_Tolerance.Length);
    return Swollen(box, 32.0 * (std::max(tolerance.Length, Fineness(shape, box)) + coarse));
  }

  /**
   *  A cell is twice the mean of the boxes' larger sides, as a floor plan's curves are met fastest
   *  in, each cell holding a few tens of entries; larger where the cells would be more than four
   *  for every curve, so that boxes scattered far apart leave few cells empty. The sides are worked
   *  in halves, which no finite box runs past every double in.
   */
  auto Gridwork::Tiled(const Grid &grid, const Tolerance &tolerance) -> Grid::Tiling
  {
    const std::size_t count = grid.m_Curves.size();
    Grid::Tiling      tiling;
    tiling.Loose.assign(count, 1);
    std::vector<Box> grown(count);
    Box              all;
    double           halves  = 0.0;
    double           gridded = 0.0;
    for (std::size_t index = 0; index < count; ++index)
    {
      const std::optional<Box> box = Grown(grid, index, tolerance);
      if (!box)
      {
        tiling.Unboxed.push_back(index);
        continue;
      }
      tiling.Loose[index] = 0;
      grown[index]        = *box;
      all                 = Joined(all, *box);
      halves  += std::max((box->Max.X / 2) - (box->Min.X / 2), (box->Max.Y / 2) - (box->Min.Y / 2));
      gridded += 1.0;
    }
    if (Zero(gridded))
    {
      return tiling;
    }
    const double most = std::min((4.0 * gridded) + 16.0, 0x1p30);
    const double wide = (all.Max.X / 2) - (all.Min.X / 2);
    const double high = (all.Max.Y / 2) - (all.Min.Y / 2);
    const double half = std::max({2.0 * halves / gridded, std::sqrt(wide) * std::sqrt(high / most),
                                  wide / most, high / most});
    tiling.Origin     = all.Min;
    tiling.Scale      = 0.5 / half;
    tiling.Columns    = Rowed(wide / half, most);
    tiling.Rows       = Rowed(high / half, most);
    tiling.Cells.assign(static_cast<std::size_t>(tiling.Columns) * tiling.Rows, {});
    /**
     *  Each cell counted, its run placed after the one before, and filled.
     */
    const auto each = [&tiling, &grown](const auto &visit) -> void
    {
      for (std::size_t index = 0; index < grown.size(); ++index)
      {
        if (tiling.Loose[index] != 0)
        {
          continue;
        }
        const Area area = Covered(tiling, grown[index]);
        for (std::size_t row = area.Low; row <= area.High; ++row)
        {
          for (std::size_t column = area.Left; column <= area.Right; ++column)
          {
            visit(
                Grid::Entry{
                    .Grown = grown[index], .Index = index, .Column = area.Left, .Row = area.Low},
                tiling.Cells[(row * tiling.Columns) + column]);
          }
        }
      }
    };
    each([](const Grid::Entry &, Grid::Cell &cell) -> void { ++cell.Count; });
    std::size_t first = 0;
    for (Grid::Cell &cell : tiling.Cells)
    {
      cell.First  = first;
      cell.Room   = cell.Count;
      first      += cell.Count;
      cell.Count  = 0;
    }
    tiling.Entries.resize(first);
    tiling.Live = first;
    each(
        [&tiling](const Grid::Entry &entry, Grid::Cell &cell) -> void
        {
          tiling.Entries[cell.First + cell.Count] = entry;
          ++cell.Count;
        });
    return tiling;
  }

  auto Gridwork::Covered(const Grid::Tiling &tiling, const Box &box) -> Area
  {
    return {.Left  = Binned((box.Min.X - tiling.Origin.X) * tiling.Scale, tiling.Columns),
            .Low   = Binned((box.Min.Y - tiling.Origin.Y) * tiling.Scale, tiling.Rows),
            .Right = Binned((box.Max.X - tiling.Origin.X) * tiling.Scale, tiling.Columns),
            .High  = Binned((box.Max.Y - tiling.Origin.Y) * tiling.Scale, tiling.Rows)};
  }

  /**
   *  The pairs in a cell are tested two by two; a pair is told only in the cell where both boxes
   *  begin, the column and row of the larger of their least corners, which both cover, as a cell's
   *  column and row never decrease along x and y. So each pair is told once, with no set of those
   *  told. Then the pairs are put in order of their first index by counting, then of their second
   *  among those of one first, and handed to their kinds in that order, so that the order is the
   *  curves' own.
   */
  auto Gridwork::Candidates(const Grid &grid, const Tolerance &tolerance) -> Pairs
  {
    const bool              kept   = tolerance.Length <= grid.m_Tolerance.Length;
    const Grid::Tiling      made   = kept ? Grid::Tiling{} : Tiled(grid, tolerance);
    const Grid::Tiling     &tiling = kept ? grid.m_Tiling : made;
    const Span              all    = {0, grid.m_Curves.size()};
    std::vector<Grid::Pair> found;
    Cellwise(tiling, {0, tiling.Cells.size()}, found);
    Loosely(tiling, all, found);
    return Sorted(grid, found, all);
  }

  void Gridwork::Cellwise(const Grid::Tiling &tiling, const Span &cells,
                          std::vector<Grid::Pair> &into)
  {
    const std::span<const Grid::Entry> entries = tiling.Entries;
    for (std::size_t place = cells[0]; place < cells[1]; ++place)
    {
      const Grid::Cell &cell                  = tiling.Cells[place];
      const auto        column                = static_cast<std::uint32_t>(place % tiling.Columns);
      const auto        row                   = static_cast<std::uint32_t>(place / tiling.Columns);
      const std::span<const Grid::Entry> held = entries.subspan(cell.First, cell.Count);
      for (auto one = held.begin(); one != held.end(); ++one)
      {
        const Box &box = one->Grown;
        for (auto other = std::next(one); other != held.end(); ++other)
        {
          /**
           *  Every test worked out, none skipped, as most pairs in a cell fail one at random.
           */
          const Box &next   = other->Grown;
          const bool across = std::max(box.Min.X, next.Min.X) <= std::min(box.Max.X, next.Max.X);
          const bool along  = std::max(box.Min.Y, next.Min.Y) <= std::min(box.Max.Y, next.Max.Y);
          const bool here   = std::max(one->Column, other->Column) == column;
          const bool there  = std::max(one->Row, other->Row) == row;
          if (static_cast<int>(across)
                  + static_cast<int>(along)
                  + static_cast<int>(here)
                  + static_cast<int>(there)
              == 4)
          {
            const auto [first, second] = std::minmax(one->Index, other->Index);
            into.push_back({first, second});
          }
        }
      }
    }
  }

  /**
   *  A pair of which a curve with no box is one has its lower index in the span where that is
   *  another curve before it, or where it is the curve itself, then paired with every curve after
   *  it. Two curves with no box are told once, from the lower.
   */
  void Gridwork::Loosely(const Grid::Tiling &tiling, const Span &span,
                         std::vector<Grid::Pair> &into)
  {
    const std::size_t count = tiling.Loose.size();
    for (const std::size_t loose : tiling.Unboxed)
    {
      for (std::size_t other = span[0]; other < std::min(span[1], loose); ++other)
      {
        if (tiling.Loose[other] == 0)
        {
          into.push_back({other, loose});
        }
      }
      if (span[0] <= loose && loose < span[1])
      {
        for (std::size_t other = loose + 1; other < count; ++other)
        {
          into.push_back({loose, other});
        }
      }
    }
  }

  auto Gridwork::Sorted(const Grid &grid, std::span<const Grid::Pair> found, const Span &span)
      -> Pairs
  {
    std::vector<std::size_t> starts(span[1] - span[0] + 1);
    for (const Grid::Pair &pair : found)
    {
      ++starts[pair[0] - span[0] + 1];
    }
    std::partial_sum(starts.begin(), starts.end(), starts.begin());
    std::vector<Grid::Pair> sorted(found.size());
    for (const Grid::Pair &pair : found)
    {
      sorted[starts[pair[0] - span[0]]] = pair;
      ++starts[pair[0] - span[0]];
    }
    /**
     *  Each start now ends its run; the runs of one first are short.
     */
    std::size_t begin = 0;
    for (std::size_t first = 0; first + 1 < starts.size(); ++first)
    {
      const std::size_t end = starts[first];
      if (end - begin > 1)
      {
        std::ranges::sort(std::span{sorted}.subspan(begin, end - begin));
      }
      begin = end;
    }
    Pairs pairs;
    for (const Grid::Pair &pair : sorted)
    {
      const Kind kind = PairKinds.at(grid.m_Shapes[pair[0]]).at(grid.m_Shapes[pair[1]]);
      pairs.at(static_cast<std::size_t>(kind)).push_back(pair);
    }
    return pairs;
  }

  /**
   *  The curve's box grown covers some cells; each is looked in once, and a curve of the grid is
   *  told in the one cell where both boxes grown begin, as for a pair of the grid's own.
   */
  auto Gridwork::Near(const Grid &grid, const Prepared &curve, const Tolerance &tolerance)
      -> std::vector<std::size_t>
  {
    return Nearby(grid, Probed(grid, curve, tolerance));
  }

  auto Gridwork::Nearby(const Grid &grid, const std::optional<Box> &probe)
      -> std::vector<std::size_t>
  {
    const Grid::Tiling      &tiling = grid.m_Tiling;
    std::vector<std::size_t> near;
    if (!probe)
    {
      near.resize(tiling.Loose.size());
      std::ranges::iota(near, std::size_t{0});
      return near;
    }
    if (!tiling.Cells.empty())
    {
      const Area area = Covered(tiling, *probe);
      for (std::size_t row = area.Low; row <= area.High; ++row)
      {
        for (std::size_t column = area.Left; column <= area.Right; ++column)
        {
          const Grid::Cell &cell = tiling.Cells[(row * tiling.Columns) + column];
          for (const Grid::Entry &entry : std::span{tiling.Entries}.subspan(cell.First, cell.Count))
          {
            const bool here =
                std::max(area.Left, entry.Column) == column && std::max(area.Low, entry.Row) == row;
            if (here && Overlaps(*probe, entry.Grown))
            {
              near.push_back(entry.Index);
            }
          }
        }
      }
    }
    near.insert(near.end(), tiling.Unboxed.begin(), tiling.Unboxed.end());
    std::ranges::sort(near);
    return near;
  }

  void Gridwork::Ungridded(Grid &grid, std::size_t index)
  {
    Grid::Tiling            &tiling = grid.m_Tiling;
    const std::optional<Box> grown  = Grown(grid, index, grid.m_Tolerance);
    if (!grown)
    {
      std::erase(tiling.Unboxed, index);
      return;
    }
    const Area area = Covered(tiling, *grown);
    for (std::size_t row = area.Low; row <= area.High; ++row)
    {
      for (std::size_t column = area.Left; column <= area.Right; ++column)
      {
        Removed(tiling, (row * tiling.Columns) + column, index);
      }
    }
  }

  /**
   *  A tiling with no cell, made where no curve had a box, is made again: the cells are sized from
   *  the boxes there are.
   */
  void Gridwork::Regridded(Grid &grid, std::size_t index)
  {
    Grid::Tiling &tiling = grid.m_Tiling;
    if (tiling.Cells.empty())
    {
      tiling = Tiled(grid, grid.m_Tolerance);
      return;
    }
    const std::optional<Box> grown = Grown(grid, index, grid.m_Tolerance);
    tiling.Loose[index]            = grown ? 0 : 1;
    if (!grown)
    {
      tiling.Unboxed.insert(std::ranges::upper_bound(tiling.Unboxed, index), index);
      return;
    }
    const Area area = Covered(tiling, *grown);
    for (std::size_t row = area.Low; row <= area.High; ++row)
    {
      for (std::size_t column = area.Left; column <= area.Right; ++column)
      {
        Inserted(tiling, (row * tiling.Columns) + column,
                 {.Grown = *grown, .Index = index, .Column = area.Left, .Row = area.Low});
      }
    }
    /**
     *  The room left by cells that moved is given back once it is more than the entries held.
     */
    if (tiling.Entries.size() > (2 * tiling.Live) + 64)
    {
      Compacted(tiling);
    }
  }

  void Gridwork::Removed(Grid::Tiling &tiling, std::size_t place, std::size_t index)
  {
    Grid::Cell                  &cell  = tiling.Cells[place];
    const std::span<Grid::Entry> held  = std::span{tiling.Entries}.subspan(cell.First, cell.Count);
    const auto                   found = std::ranges::find(held, index, &Grid::Entry::Index);
    *found                             = held.back();
    --cell.Count;
    --tiling.Live;
  }

  void Gridwork::Inserted(Grid::Tiling &tiling, std::size_t place, const Grid::Entry &entry)
  {
    Grid::Cell &cell = tiling.Cells[place];
    if (cell.Count == cell.Room)
    {
      const std::size_t first = tiling.Entries.size();
      const std::size_t room  = (2 * cell.Room) + 4;
      tiling.Entries.resize(first + room);
      const std::span<Grid::Entry> entries = tiling.Entries;
      std::ranges::copy(entries.subspan(cell.First, cell.Count), entries.subspan(first).begin());
      cell.First = first;
      cell.Room  = room;
    }
    tiling.Entries[cell.First + cell.Count] = entry;
    ++cell.Count;
    ++tiling.Live;
  }

  void Gridwork::Compacted(Grid::Tiling &tiling)
  {
    std::vector<Grid::Entry> entries;
    entries.reserve(tiling.Live);
    const std::span<const Grid::Entry> held = tiling.Entries;
    for (Grid::Cell &cell : tiling.Cells)
    {
      const std::span<const Grid::Entry> run = held.subspan(cell.First, cell.Count);
      cell.First                             = entries.size();
      cell.Room                              = cell.Count;
      entries.insert(entries.end(), run.begin(), run.end());
    }
    tiling.Entries = std::move(entries);
  }

  Grid::Grid(std::span<const Curve> curves, const Tolerance &tolerance)
      : m_Curves(curves.begin(), curves.end()), m_Tolerance(tolerance)
  {
    const std::size_t count = m_Curves.size();
    m_Prepared.reserve(count);
    for (const Curve &curve : m_Curves)
    {
      m_Prepared.emplace_back(curve);
    }
    for (std::vector<double> *values : {&m_MinX, &m_MinY, &m_MaxX, &m_MaxY, &m_Grain})
    {
      values->resize(count);
    }
    m_Shapes.resize(count);
    for (std::size_t index = 0; index < count; ++index)
    {
      Gridwork::Measured(*this, index);
    }
    m_Tiling = Gridwork::Tiled(*this, m_Tolerance);
  }

  Grid::Grid(const Grid &other) : Grid(other.m_Curves, other.m_Tolerance) {}

  auto Grid::operator=(const Grid &other) -> Grid &
  {
    Grid copy{other};
    *this = std::move(copy);
    return *this;
  }

  void Grid::Moved(std::size_t index, Curve curve)
  {
    Curve &kept = m_Curves.at(index);
    Gridwork::Ungridded(*this, index);
    kept              = std::move(curve);
    m_Prepared[index] = Prepared{kept};
    Gridwork::Measured(*this, index);
    Gridwork::Regridded(*this, index);
  }

  auto Grid::Near(const Box &box) const -> std::vector<std::size_t>
  {
    return Gridwork::Nearby(*this, box);
  }

  auto Grid::Window(const Box &box) const -> std::vector<std::size_t>
  {
    std::vector<std::size_t> near = Near(box);
    std::erase_if(near, [this, &box](std::size_t index) -> bool
                  { return !Inside(m_Curves[index], box); });
    return near;
  }

  auto Grid::Crossing(const Box &box) const -> std::vector<std::size_t>
  {
    std::vector<std::size_t> near = Near(box);
    std::erase_if(near, [this, &box](std::size_t index) -> bool
                  { return !Intersects(m_Curves[index], box); });
    return near;
  }

  auto Grid::Window(const Region &region) const -> std::vector<std::size_t>
  {
    std::vector<std::size_t> near = Near(Encased(region));
    std::erase_if(near, [this, &region](std::size_t index) -> bool
                  { return !Inside(m_Curves[index], region, m_Tolerance); });
    return near;
  }

  auto Grid::Crossing(const Region &region) const -> std::vector<std::size_t>
  {
    std::vector<std::size_t> near = Near(Encased(region));
    std::erase_if(near, [this, &region](std::size_t index) -> bool
                  { return !Intersects(m_Curves[index], region, m_Tolerance); });
    return near;
  }

  /**
   *  Lines are told from their kind alone, as the grid's shapes say: no curve is asked what it is
   *  again, and the loop over them calls nothing else.
   */
  void Gridwork::Met(const Grid &grid, const Run &run, const Grid::Sink &sink)
  {
    Teller     teller{run.Within, sink};
    const bool lines = run.Of == Kind::Square || run.Of == Kind::Lines;
    if (lines)
    {
      for (const Grid::Pair &pair : run.Pairs)
      {
        teller.Lines({grid[pair[0]], grid[pair[1]]}, pair);
      }
      return;
    }
    for (const Grid::Pair &pair : run.Pairs)
    {
      teller.Paired({grid[pair[0]], grid[pair[1]]}, pair);
    }
  }

  void Gridwork::Alone(const Grid &grid, const Tolerance &tolerance, const Grid::Sink &sink)
  {
    const Pairs pairs = Candidates(grid, tolerance);
    for (std::size_t kind = 0; kind < Kinds; ++kind)
    {
      Met(grid, {.Of = static_cast<Kind>(kind), .Pairs = pairs.at(kind), .Within = tolerance},
          sink);
    }
  }

  auto Gridwork::Together(const Grid::By &by, const Tolerance &tolerance, const Grid::Sink &sink)
      -> bool
  {
    const Grid        &grid   = by.first;
    Workers           &pool   = by.second;
    const bool         kept   = tolerance.Length <= grid.m_Tolerance.Length;
    const Grid::Tiling made   = kept ? Grid::Tiling{} : Tiled(grid, tolerance);
    const Search       search = {.Of    = &grid,
                                 .Tiles = kept ? &grid.m_Tiling : &made,
                                 .Grown = kept ? grid.m_Tolerance : tolerance};
    Shift              shift{search, tolerance, pool.Size()};
    const bool         helps = pool.Size() < Helping;
    return Pooled::Run(pool,
                       [&shift, &sink, helps](std::size_t slot) -> void
                       {
                         shift.Find();
                         shift.Seek();
                         if (slot != 0)
                         {
                           shift.Intersect();
                           return;
                         }
                         if (helps)
                         {
                           shift.Intersect();
                         }
                         shift.Tell(sink);
                       });
  }

  void IntersectAll(const Grid &grid, const Tolerance &tolerance, const Grid::Sink &sink)
  {
    const std::size_t threads =
        std::min<std::size_t>(std::thread::hardware_concurrency(), grid.Size() / Gridwork::Share);
    if (threads < 2)
    {
      Gridwork::Alone(grid, tolerance, sink);
      return;
    }
    Workers pool{threads};
    IntersectAll({grid, pool}, tolerance, sink);
  }

  void IntersectAll(const Grid::By &by, const Tolerance &tolerance, const Grid::Sink &sink)
  {
    const Workers &pool = by.second;
    if (pool.Size() < 2 || !Gridwork::Together(by, tolerance, sink))
    {
      Gridwork::Alone(by.first, tolerance, sink);
    }
  }

  void IntersectWith(const Grid::With &with, const Tolerance &tolerance, const Grid::Sink &sink)
  {
    const Prepared &curve = with.first;
    const Grid     &grid  = with.second;
    Teller          teller{tolerance, sink};
    const bool      line = std::holds_alternative<Line>(curve.Shape());
    for (const std::size_t index : Gridwork::Near(grid, curve, tolerance))
    {
      const Prepared::Pair curves = {curve, grid[index]};
      const Grid::Pair     pair   = {Grid::Outside, index};
      if (line && std::holds_alternative<Line>(grid[index].Shape()))
      {
        teller.Lines(curves, pair);
        continue;
      }
      teller.Paired(curves, pair);
    }
  }
} // namespace Geoxy
