/***************************************************************************************************
 *  @file      source/weave/merge.cpp
 *  @brief     Cuts carried across shared stretches, and vertices within the tolerance merged in a
 *             stated order.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "weave/merge.hpp"

#include "scratch.hpp"
#include "unions.hpp"
#include "weave/found.hpp"

#include "geoxy/box.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/intersect.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <span>
#include <vector>

namespace Geoxy::Weave
{
  namespace
  {
    /**
     *  The cuts of one curve of a shared stretch inside it carried to the other, each at the
     *  other's place the stretch maps it to, with the same vertex.
     *
     *  @param  grouped Each curve's cuts.
     *  @param  common  The stretch.
     *  @param  back    Whether from the second curve to the first.
     *  @return         Whether a cut was added.
     */
    auto Carried(std::vector<std::vector<Notch>> &grouped, const Common &common, bool back) -> bool
    {
      const std::size_t from  = back ? common.Second : common.First;
      const std::size_t onto  = back ? common.First : common.Second;
      const Range       span  = back ? common.Shared.Second : common.Shared.First;
      bool              added = false;
      for (const Notch notch : std::vector<Notch>(grouped[from]))
      {
        const bool inside = span.Start < notch.Along && notch.Along < span.End;
        const bool known  = std::ranges::any_of(grouped[onto], [&notch](const Notch &other) -> bool
                                                { return other.Vertex == notch.Vertex; });
        if (!inside || known)
        {
          continue;
        }
        grouped[onto].push_back({.Source = onto,
                                 .Along  = Reflected(common, notch.Along, back),
                                 .Vertex = notch.Vertex});
        added = true;
      }
      return added;
    }

    /**
     *  The column or row of a grid a coordinate falls in.
     *
     *  @param  coordinate The coordinate.
     *  @param  cell       The cells' side, no less than 16 steps of a double at the coordinates'
     *                     size, so that the quotient is far inside what an integer holds.
     *  @return            The column or row.
     */
    auto Celled(double coordinate, double cell) -> std::int64_t
    {
      return static_cast<std::int64_t>(std::floor(coordinate / cell));
    }

    /**
     *  A cell of a grid, as a key.
     *
     *  @param  column The cell's column.
     *  @param  row    Its row.
     *  @return        The key.
     */
    auto Keyed(std::int64_t column, std::int64_t row) -> std::uint64_t
    {
      return (static_cast<std::uint64_t>(column) * 0x9E37'79B9'7F4A'7C15U)
           ^ static_cast<std::uint64_t>(row);
    }

    /**
     *  The vertices in each cell of a grid: a table of the cells' keys, open and probed in turn,
     *  each slot holding its key beside the cell's last vertex, and each cell's vertices linked
     *  from the last added; all of it in the thread's scratch memory, asked for once.
     */
    class Cells
    {
    public:
      /**
       *  A table for up to some vertices of some more, at most half full.
       *
       *  @param count How many vertices there are.
       *  @param room  How many of them it may hold.
       *  @param cell  The cells' side.
       */
      Cells(std::size_t count, std::size_t room, double cell)
          : m_Slots(std::bit_ceil((2 * room) + 2), Slot{.Last = count}),
            m_Marks(std::max<std::size_t>(m_Slots.size() / 4, 1), 0), m_Links(count, count),
            m_None(count), m_Cell(cell)
      {
      }

      /**
       *  The first earlier vertex that stands for itself within a length of a vertex, looked for in
       *  the cells the square of its length about it reaches; or the vertex itself, added to its
       *  cell, where there is none.
       *
       *  @param  points The vertices' points.
       *  @param  index  The vertex.
       *  @param  length The length, a quarter of a cell's side.
       *  @return        The first such vertex, or the vertex.
       */
      auto Merged(std::span<const Point> points, std::size_t index, double length) -> std::size_t
      {
        const Point        point  = points[index];
        const std::int64_t left   = Celled(point.X - length, m_Cell);
        const std::int64_t right  = Celled(point.X + length, m_Cell);
        const std::int64_t bottom = Celled(point.Y - length, m_Cell);
        const std::int64_t top    = Celled(point.Y + length, m_Cell);
        std::size_t        first  = index;
        for (std::int64_t column = left; column <= right; ++column)
        {
          for (std::int64_t row = bottom; row <= top; ++row)
          {
            for (std::size_t other = First(Keyed(column, row)); other < index;
                 other             = m_Links[other])
            {
              if (Distance(points[other], point) <= length)
              {
                first = std::min(first, other);
              }
            }
          }
        }
        if (first == index)
        {
          /**
           *  A point's own cell lies between those of the point less and more the length, rounded
           *  each way, as rounding and the floor of a quotient keep the order: so it is theirs
           *  where they are one.
           */
          const std::int64_t column = left == right ? left : Celled(point.X, m_Cell);
          const std::int64_t row    = bottom == top ? bottom : Celled(point.Y, m_Cell);
          Add(Keyed(column, row), index);
        }
        return first;
      }

    private:
      /**
       *  The vertex added last to a cell: none where the cell's bit in the sieve is clear, as no
       *  vertex was added to it, else the table's.
       *
       *  @param  key The cell's key.
       *  @return     The vertex; past every vertex for a cell with none.
       */
      [[nodiscard]] auto First(std::uint64_t key) const -> std::size_t
      {
        const std::size_t bit = Sieved(key);
        if ((m_Marks[bit / 64] & (std::uint64_t{1} << (bit % 64))) == 0)
        {
          return m_None;
        }
        return m_Slots[Probed(key)].Last;
      }

      /**
       *  A vertex added to a cell.
       *
       *  @param key    The cell's key.
       *  @param vertex The vertex.
       */
      void Add(std::uint64_t key, std::size_t vertex)
      {
        const std::size_t bit  = Sieved(key);
        m_Marks[bit / 64]     |= std::uint64_t{1} << (bit % 64);
        Slot &slot             = m_Slots[Probed(key)];
        slot.Key               = key;
        m_Links[vertex]        = slot.Last;
        slot.Last              = vertex;
      }

      /**
       *  A slot of the table: a cell's key and its last vertex, or no vertex.
       */
      struct Slot
      {
        std::uint64_t Key  = 0; /**< The cell's key, where it holds one. */
        std::size_t   Last = 0; /**< The cell's last vertex, or none.    */
      };

      /**
       *  The slot of a cell: where its key is, or the empty slot it would go in.
       *
       *  @param  key The cell's key.
       *  @return     The slot.
       */
      [[nodiscard]] auto Probed(std::uint64_t key) const -> std::size_t
      {
        const std::size_t mask = m_Slots.size() - 1;
        std::size_t       slot = (key ^ (key >> 29U)) & mask;
        while (m_Slots[slot].Last != m_None && m_Slots[slot].Key != key)
        {
          slot = (slot + 1) & mask;
        }
        return slot;
      }

      /**
       *  A cell's bit in the sieve: its key mixed, its top bits, sixteen bits for each slot.
       *
       *  @param  key The cell's key.
       *  @return     The bit.
       */
      [[nodiscard]] auto Sieved(std::uint64_t key) const -> std::size_t
      {
        const auto bits = static_cast<unsigned>(std::countr_zero(m_Marks.size() * 64));
        return static_cast<std::size_t>((key * 0xD6E8'FEB8'6659'FD93U) >> (64U - bits));
      }

      List<Slot>          m_Slots = {}; /**< The slots. */
      List<std::uint64_t> m_Marks = {}; /**< A bit for each of many cells, set where
                                             a vertex was added to one of them.   */
      List<std::size_t> m_Links = {};   /**< Each vertex's one before, in its cell. */
      std::size_t       m_None  = 0;    /**< The index that stands for no vertex.   */
      double            m_Cell  = 0.0;  /**< The cells' side.                       */
    };

    /**
     *  A grid of as many cells along each side over a square from a corner: a point's cell by its
     *  offsets from the corner times a scale, floored and held within the grid, which never falls
     *  as a coordinate rises.
     */
    struct Lattice
    {
      Point       Origin = {};  /**< The least corner.            */
      double      Scale  = 0.0; /**< Cells for each unit along.   */
      std::size_t Side   = 1;   /**< How many cells along a side. */

      /**
       *  The column or row of an offset from the corner.
       *
       *  @param  offset The offset.
       *  @return        The column or row.
       */
      [[nodiscard]] auto Along(double offset) const -> std::size_t
      {
        const double place = std::min(offset * Scale, static_cast<double>(Side - 1));
        return static_cast<std::size_t>(std::max(place, 0.0));
      }

      /**
       *  A point's cell, row by row.
       *
       *  @param  point The point.
       *  @return       The cell.
       */
      [[nodiscard]] auto Cell(Point point) const -> std::size_t
      {
        return (Along(point.Y - Origin.Y) * Side) + Along(point.X - Origin.X);
      }

      /**
       *  Whether the square of a reach about a point, each side its coordinates less and more the
       *  reach, lies in more than one cell.
       *
       *  @param  point The point.
       *  @param  reach The reach.
       *  @return       Whether it does.
       */
      [[nodiscard]] auto Straddles(Point point, double reach) const -> bool
      {
        return Along(point.X - reach - Origin.X) != Along(point.X + reach - Origin.X)
            || Along(point.Y - reach - Origin.Y) != Along(point.Y + reach - Origin.Y);
      }
    };

    /**
     *  A grid over a box of about some cells, as many along each side: the box's wider side split
     *  so, or no less than the least normal double, its scale held finite, so that no offset times
     *  it is undefined.
     *
     *  @param  box   The box.
     *  @param  cells About how many cells.
     *  @return       The grid.
     */
    auto Latticed(const Box &box, std::size_t cells) -> Lattice
    {
      const double wide  = std::max(box.Max.X - box.Min.X, box.Max.Y - box.Min.Y);
      const auto   side  = static_cast<double>(std::max<std::size_t>(
          1, static_cast<std::size_t>(std::sqrt(static_cast<double>(cells)))));
      const double least = std::numeric_limits<double>::min();
      return {.Origin = box.Min,
              .Scale  = std::min(side / std::max(wide, least), std::numeric_limits<double>::max()),
              .Side   = static_cast<std::size_t>(side)};
    }

    /**
     *  How many vertices, at the least, Crowds tells apart: fewer keep Merge's table of cells
     *  within the caches.
     */
    constexpr std::size_t Throng = std::size_t{1} << 16U;

    /**
     *  About how many vertices Crowds puts in a bucket: a bucket's grid stays within the caches.
     */
    constexpr std::size_t Flock = 1024;

    /**
     *  What Huddled works with: the reach two vertices are told near within, and room for a
     *  bucket's grid.
     */
    struct Huddle
    {
      double            Reach  = 0.0; /**< The reach.                        */
      List<std::size_t> Starts = {};  /**< Where each cell's vertices start. */
      List<std::size_t> Cells  = {};  /**< Each vertex's cell.               */
      List<std::size_t> Dealt  = {};  /**< The vertices, cell by cell.       */
    };

    /**
     *  A bucket's vertices marked near where another of them may lie within the reach: laid in a
     *  grid of about a cell each, and marked where the square of the reach about one lies in more
     *  than one cell, where its cell holds more than eight, or where another in its cell lies
     *  within the reach along both axes.
     *
     *  @param spots  The bucket's points.
     *  @param near   Whether each is near, marked.
     *  @param huddle The reach, and room for the grid.
     */
    void Huddled(std::span<const Point> spots, std::span<std::uint8_t> near, Huddle &huddle)
    {
      Box box = {.Min = spots.front(), .Max = spots.front()};
      for (const Point spot : spots)
      {
        box = Grown(box, spot);
      }
      const Lattice      fine   = Latticed(box, spots.size());
      const std::size_t  cells  = fine.Side * fine.Side;
      List<std::size_t> &starts = huddle.Starts;
      starts.assign(cells + 1, 0);
      huddle.Cells.resize(spots.size());
      huddle.Dealt.resize(spots.size());
      for (std::size_t slot = 0; slot < spots.size(); ++slot)
      {
        huddle.Cells[slot] = fine.Cell(spots[slot]);
        ++starts[huddle.Cells[slot] + 1];
      }
      std::partial_sum(starts.begin(), starts.end(), starts.begin());
      for (std::size_t slot = 0; slot < spots.size(); ++slot)
      {
        huddle.Dealt[starts[huddle.Cells[slot]]++] = slot;
      }
      std::size_t from = 0;
      for (std::size_t cell = 0; cell < cells; ++cell)
      {
        const std::span<const std::size_t> held =
            std::span<const std::size_t>(huddle.Dealt).subspan(from, starts[cell] - from);
        from = starts[cell];
        for (std::size_t one = 0; one < held.size(); ++one)
        {
          const Point spot = spots[held[one]];
          near[held[one]] |=
              static_cast<std::uint8_t>(held.size() > 8 || fine.Straddles(spot, huddle.Reach));
          for (std::size_t two = one + 1; two < held.size() && held.size() <= 8; ++two)
          {
            const Vector gap = spots[held[two]] - spot;
            if (std::abs(gap.X) <= huddle.Reach && std::abs(gap.Y) <= huddle.Reach)
            {
              near[held[one]] = 1;
              near[held[two]] = 1;
            }
          }
        }
      }
    }

    /**
     *  Which vertices no twin may have another such within a length, of Throng or more: those
     *  Huddled marks in buckets of about Flock each, by a coarse grid over them all, and those
     *  whose square of twice the length lies in more than one bucket.
     *
     *  Two whose distance Merge finds within the length lie within twice it along both axes, as a
     *  distance is no less than either gap but for its rounding: where they are in two buckets, or
     *  two cells of one, each one's square lies in more than one, as the grids' cells never fall as
     *  a coordinate rises, and else Huddled weighs them, so both are marked. Twins are never looked
     *  for, nor looked at.
     *
     *  @param  found  The findings: the vertices and their twins.
     *  @param  length The length.
     *  @return        Whether each vertex is marked: none of the twins.
     */
    auto Crowds(const Findings &found, double length) -> List<std::uint8_t>
    {
      const std::span<const Point> points = found.Points;
      const std::size_t            count  = points.size();
      List<std::uint8_t>           crowded(count, 0);
      Box                          bounds = {};
      std::size_t                  loners = 0;
      for (std::size_t index = 0; index < count; ++index)
      {
        if (found.Twins[index] >= index)
        {
          crowded[index] = 1;
          bounds         = loners == 0 ? Box{.Min = points[index], .Max = points[index]}
                                       : Grown(bounds, points[index]);
          ++loners;
        }
      }
      const Lattice     coarse = Latticed(bounds, loners / Flock);
      List<std::size_t> starts((coarse.Side * coarse.Side) + 1, 0);
      List<std::size_t> cells(count);
      for (std::size_t index = 0; index < count; ++index)
      {
        if (crowded[index] != 0)
        {
          cells[index] = coarse.Cell(points[index]);
          ++starts[cells[index] + 1];
        }
      }
      std::partial_sum(starts.begin(), starts.end(), starts.begin());
      Huddle             huddle = {.Reach = 2.0 * length};
      List<std::size_t>  order(loners);
      List<Point>        spots(loners);
      List<std::uint8_t> near(loners, 0);
      for (std::size_t index = 0; index < count; ++index)
      {
        if (crowded[index] != 0)
        {
          const std::size_t slot = starts[cells[index]]++;
          order[slot]            = index;
          spots[slot]            = points[index];
          near[slot] = static_cast<std::uint8_t>(coarse.Straddles(points[index], huddle.Reach));
        }
      }
      std::size_t from = 0;
      for (std::size_t bucket = 0; bucket + 1 < starts.size(); ++bucket)
      {
        const std::size_t size = starts[bucket] - from;
        if (size > 1)
        {
          Huddled(std::span<const Point>(spots).subspan(from, size),
                  std::span<std::uint8_t>(near).subspan(from, size), huddle);
        }
        from = starts[bucket];
      }
      for (std::size_t slot = 0; slot < loners; ++slot)
      {
        crowded[order[slot]] = near[slot];
      }
      return crowded;
    }
  } // namespace

  void Carry(Findings &found)
  {
    if (found.Commons.empty())
    {
      return;
    }
    std::vector<std::vector<Notch>> grouped(found.Count);
    for (std::size_t curve = 0; curve < found.Count; ++curve)
    {
      const auto [ends, count]         = Ended(found, curve);
      const std::span<const Notch> own = std::span<const Notch>(ends).first(count);
      grouped[curve].assign(own.begin(), own.end());
    }
    for (const Notch &notch : found.Notches)
    {
      grouped[notch.Source].push_back(notch);
    }
    List<std::size_t> held(found.Count);
    for (std::size_t curve = 0; curve < found.Count; ++curve)
    {
      held[curve] = grouped[curve].size();
    }
    std::size_t added = 1;
    while (added > 0)
    {
      added = 0;
      for (const Common &common : found.Commons)
      {
        added += static_cast<std::size_t>(Carried(grouped, common, false));
        added += static_cast<std::size_t>(Carried(grouped, common, true));
      }
    }
    for (std::size_t curve = 0; curve < found.Count; ++curve)
    {
      const std::vector<Notch> &notches = grouped[curve];
      found.Notches.insert(found.Notches.end(),
                           notches.begin() + static_cast<std::ptrdiff_t>(held[curve]),
                           notches.end());
    }
  }

  auto Merge(const Findings &found, double length, Unions &sets) -> Merging
  {
    const std::span<const Point> points = found.Points;
    const double                 cell   = 4.0 * length;
    const auto                   room   = static_cast<std::size_t>(std::ranges::count_if(
        found.Twins, [&points](std::size_t twin) -> bool { return twin >= points.size(); }));
    const List<std::uint8_t> crowded = room < Throng ? List<std::uint8_t>{} : Crowds(found, length);
    Cells   cells{points.size(),
                  crowded.empty() ? room : static_cast<std::size_t>(std::ranges::count(crowded, 1)),
                  cell};
    Merging merging          = {.Ids = List<std::size_t>(points.size())};
    List<std::size_t> &ids   = merging.Ids;
    std::size_t        count = 0;
    for (std::size_t index = 0; index < points.size(); ++index)
    {
      const std::size_t twin = found.Twins[index];
      if (twin < index)
      {
        sets.Join(twin, index);
        ids[index] = ids[twin];
        continue;
      }
      const std::size_t first =
          crowded.empty() || crowded[index] != 0 ? cells.Merged(points, index, length) : index;
      if (first == index)
      {
        ids[index] = count++;
      }
      else
      {
        sets.Join(first, index);
        ids[index] = ids[first];
        if (!Twinned(points[first], points[index]))
        {
          merging.Moved.push_back(index);
        }
      }
    }
    return merging;
  }
} // namespace Geoxy::Weave
