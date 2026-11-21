/***************************************************************************************************
 *  @file      source/belt/work.hpp
 *  @brief     What a polyline's offset is worked out in: its pieces, the curves made from them,
 *             where they are cut and the points they share.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_BELT_WORK_HPP
#define GEOXY_BELT_WORK_HPP

#include "belt/carrier.hpp"
#include "belt/piece.hpp"
#include "scratch.hpp"

#include "geoxy/box.hpp"
#include "geoxy/offset.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

namespace Geoxy::Belt
{
  /**
   *  How many pieces of a shape at the most are each measured from a point, with no grid.
   */
  inline constexpr std::size_t Handful = 24;

  /**
   *  A grid of cells over a box, each holding the numbers of the boxes that touch it.
   */
  class Stalls
  {
  public:
    /**
     *  The grid laid over a box, its cells about as many as the boxes, no smaller than a 512th of
     *  the box's longer side, and every box put in the cells it touches.
     *
     *  @param all   The box it covers.
     *  @param boxes The boxes.
     */
    void Lay(const Box &all, std::span<const Box> boxes)
    {
      const double wide = std::max(all.Max.X - all.Min.X, std::numeric_limits<double>::min());
      const double high = std::max(all.Max.Y - all.Min.Y, std::numeric_limits<double>::min());
      m_Cell            = std::max(std::sqrt((wide * high) / static_cast<double>(boxes.size())),
                                   std::max(wide, high) / 512.0);
      m_Min             = all.Min;
      m_Wide            = static_cast<std::size_t>(std::max(1.0, std::ceil(wide / m_Cell)));
      m_High            = static_cast<std::size_t>(std::max(1.0, std::ceil(high / m_Cell)));
      m_Starts.assign((m_Wide * m_High) + 1, 0);
      m_Spans.clear();
      m_Spans.reserve(boxes.size());
      for (const Box &box : boxes)
      {
        m_Spans.push_back({static_cast<std::uint32_t>(Column(box.Min.X)),
                           static_cast<std::uint32_t>(Column(box.Max.X)),
                           static_cast<std::uint32_t>(Row(box.Min.Y)),
                           static_cast<std::uint32_t>(Row(box.Max.Y))});
      }
      const auto each = [this](auto &&put) -> void
      {
        for (std::size_t index = 0; index < m_Spans.size(); ++index)
        {
          const std::array<std::uint32_t, 4> &span = m_Spans[index];
          for (std::size_t row = span[2]; row <= span[3]; ++row)
          {
            for (std::size_t column = span[0]; column <= span[1]; ++column)
            {
              put((row * m_Wide) + column, index);
            }
          }
        }
      };
      each([this](std::size_t cell, std::size_t /*index*/) -> void { ++m_Starts[cell + 1]; });
      for (std::size_t cell = 0; cell + 1 < m_Starts.size(); ++cell)
      {
        m_Starts[cell + 1] += m_Starts[cell];
      }
      m_Items.resize(m_Starts.back());
      List<std::size_t> next(m_Starts.begin(), m_Starts.end() - 1);
      each([this, &next](std::size_t cell, std::size_t index) -> void
           { m_Items[next[cell]++] = index; });
    }

    /**
     *  The column a coordinate falls in, the first or the last for one off the grid.
     *
     *  @param  x The coordinate.
     *  @return   The column.
     */
    [[nodiscard]] auto Column(double x) const -> std::size_t
    {
      return Cell((x - m_Min.X) / m_Cell, m_Wide);
    }

    /**
     *  The row a coordinate falls in, the first or the last for one off the grid.
     *
     *  @param  y The coordinate.
     *  @return   The row.
     */
    [[nodiscard]] auto Row(double y) const -> std::size_t
    {
      return Cell((y - m_Min.Y) / m_Cell, m_High);
    }

    /**
     *  The boxes that touch a cell.
     *
     *  @param  row    The cell's row.
     *  @param  column Its column.
     *  @return        Their numbers.
     */
    [[nodiscard]] auto Held(std::size_t row, std::size_t column) const
        -> std::span<const std::size_t>
    {
      const std::size_t cell = (row * m_Wide) + column;
      return std::span{m_Items}.subspan(m_Starts[cell], m_Starts[cell + 1] - m_Starts[cell]);
    }

    /**
     *  The boxes that touch the cell a point lies in.
     *
     *  @param  point The point.
     *  @return       Their numbers; for a point off the grid, those of the cell nearest it.
     */
    [[nodiscard]] auto Around(Point point) const -> std::span<const std::size_t>
    {
      return Held(Row(point.Y), Column(point.X));
    }

    /**
     *  Whether a cell is the first, by row and then by column, that two boxes both touch: so that
     *  two that share several cells are tried in one alone.
     *
     *  @param  pair The two boxes' numbers.
     *  @param  cell The cell's row and column.
     *  @return      Whether it is.
     */
    [[nodiscard]] auto First(const std::array<std::size_t, 2> &pair,
                             const std::array<std::size_t, 2> &cell) const -> bool
    {
      const std::array<std::uint32_t, 4> &one   = m_Spans[pair[0]];
      const std::array<std::uint32_t, 4> &other = m_Spans[pair[1]];
      return std::max(one[2], other[2]) == cell[0] && std::max(one[0], other[0]) == cell[1];
    }

    /**
     *  How many columns the grid has.
     *
     *  @return The count.
     */
    [[nodiscard]] auto Wide() const -> std::size_t { return m_Wide; }

    /**
     *  How many rows the grid has.
     *
     *  @return The count.
     */
    [[nodiscard]] auto High() const -> std::size_t { return m_High; }

  private:
    /**
     *  A cell's number along one axis.
     *
     *  @param  along How many cells along the axis the coordinate lies.
     *  @param  count How many cells the axis has.
     *  @return       The cell, the first or the last for a coordinate off the grid.
     */
    static auto Cell(double along, std::size_t count) -> std::size_t
    {
      return static_cast<std::size_t>(
          std::clamp(std::floor(along), 0.0, static_cast<double>(count - 1)));
    }

    Point             m_Min    = {};  /**< The corner the cells are counted from. */
    double            m_Cell   = 1.0; /**< A cell's side.                         */
    std::size_t       m_Wide   = 1;   /**< How many columns.                      */
    std::size_t       m_High   = 1;   /**< How many rows.                         */
    List<std::size_t> m_Starts = {};  /**< Where each cell's boxes start among the items,
                                           and one past the last cell's.                 */
    List<std::size_t>                  m_Items = {}; /**< The boxes' numbers, cell after cell. */
    List<std::array<std::uint32_t, 4>> m_Spans = {}; /**< Each box's first and last column,
                                                          and first and last row.          */
  };

  /**
   *  Points taken for one when within a reach of each other along both axes, each given a number. A
   *  table finds them by the cell they fall in, the cells a thousand reaches wide: so a point is
   *  looked for in its own cell, and in the next only where it lies within a reach of it.
   */
  class Nodes
  {
  public:
    /**
     *  No points yet, and room for about as many as are likely.
     *
     *  @param reach  The reach.
     *  @param likely About how many points there will be.
     */
    Nodes(double reach, std::size_t likely) : m_Reach(reach), m_Per(1.0 / (reach * 1024.0))
    {
      m_At.reserve(likely);
      m_Slots.assign(std::max(Least, std::bit_ceil(2 * likely)), 0);
    }

    /**
     *  The number of a point: that of one already numbered within the reach of it, else a new one.
     *
     *  @param  point The point.
     *  @return       Its number.
     */
    auto Of(Point point) -> std::size_t
    {
      const double across = point.X * m_Per;
      const double upward = point.Y * m_Per;
      const double column = std::floor(across);
      const double row    = std::floor(upward);
      /**
       *  Which neighbouring cells the point is within a reach of: the one before or after along
       *  each axis, or neither.
       */
      const double edge  = 1.0 / 1024.0;
      const double aside = Next(across - column, edge);
      const double above = Next(upward - row, edge);
      std::size_t  found = Sought(point, {column, row});
      if (found == None && !(std::abs(aside) < 1.0))
      {
        found = Sought(point, {column + aside, row});
      }
      if (found == None && !(std::abs(above) < 1.0))
      {
        found = Sought(point, {column, row + above});
      }
      if (found == None && !(std::abs(aside) < 1.0) && !(std::abs(above) < 1.0))
      {
        found = Sought(point, {column + aside, row + above});
      }
      if (found != None)
      {
        return found;
      }
      m_At.push_back(point);
      if (m_At.size() * 2 > m_Slots.size())
      {
        Grown();
      }
      else
      {
        Kept(m_At.size() - 1);
      }
      return m_At.size() - 1;
    }

    /**
     *  A numbered point put where it is found to lie, within the reach of where it was.
     *
     *  @param node  Its number.
     *  @param point Where it lies.
     */
    void Move(std::size_t node, Point point) { m_At[node] = point; }

    /**
     *  The point of a number.
     *
     *  @param  node The number.
     *  @return      The point first given it.
     */
    [[nodiscard]] auto At(std::size_t node) const -> Point { return m_At[node]; }

    /**
     *  How many points are numbered.
     *
     *  @return The count.
     */
    [[nodiscard]] auto Count() const -> std::size_t { return m_At.size(); }

  private:
    /**
     *  How many slots the table starts with, a power of two.
     */
    static constexpr std::size_t Least = 256;

    /**
     *  Which neighbouring cell along an axis a place in a cell is within a reach of.
     *
     *  @param  place Where in the cell, 0 to 1.
     *  @param  edge  The reach, as a share of the cell.
     *  @return       -1 for the cell before, 1 for the one after, 0 for neither.
     */
    static auto Next(double place, double edge) -> double
    {
      if (place <= edge)
      {
        return -1.0;
      }
      return place >= 1.0 - edge ? 1.0 : 0.0;
    }

    /**
     *  The slot a cell's points are looked for from.
     *
     *  @param  cell The cell's column and row.
     *  @return      The slot.
     */
    [[nodiscard]] auto Slot(const std::array<double, 2> &cell) const -> std::size_t
    {
      const auto          first  = static_cast<std::uint64_t>(static_cast<std::int64_t>(cell[0]));
      const auto          second = static_cast<std::uint64_t>(static_cast<std::int64_t>(cell[1]));
      const std::uint64_t mixed =
          ((first * 0x9E3779B97F4A7C15ULL) ^ (second * 0xC2B2AE3D27D4EB4FULL)) >> 20U;
      return static_cast<std::size_t>(mixed) & (m_Slots.size() - 1);
    }

    /**
     *  The cell a point falls in.
     *
     *  @param  point The point.
     *  @return       Its column and row.
     */
    [[nodiscard]] auto Celled(Point point) const -> std::array<double, 2>
    {
      return {std::floor(point.X * m_Per), std::floor(point.Y * m_Per)};
    }

    /**
     *  The number of a point already numbered that lies within the reach of a point, among those of
     *  one cell.
     *
     *  @param  point The point.
     *  @param  cell  The cell's column and row.
     *  @return       The number, or none.
     */
    [[nodiscard]] auto Sought(Point point, const std::array<double, 2> &cell) const -> std::size_t
    {
      for (std::size_t slot = Slot(cell); m_Slots[slot] != 0;
           slot             = (slot + 1) & (m_Slots.size() - 1))
      {
        const Point other = m_At[m_Slots[slot] - 1];
        if (std::abs(other.X - point.X) <= m_Reach && std::abs(other.Y - point.Y) <= m_Reach)
        {
          return m_Slots[slot] - 1;
        }
      }
      return None;
    }

    /**
     *  A numbered point put in the table.
     *
     *  @param node Its number.
     */
    void Kept(std::size_t node)
    {
      std::size_t slot = Slot(Celled(m_At[node]));
      while (m_Slots[slot] != 0)
      {
        slot = (slot + 1) & (m_Slots.size() - 1);
      }
      m_Slots[slot] = static_cast<std::uint32_t>(node + 1);
    }

    /**
     *  The table made twice as large, every point put in it again.
     */
    void Grown()
    {
      m_Slots.assign(m_Slots.size() * 2, 0);
      for (std::size_t node = 0; node < m_At.size(); ++node)
      {
        Kept(node);
      }
    }

    double              m_Reach = 0.0; /**< The reach.              */
    double              m_Per   = 0.0; /**< One over a cell's side. */
    List<Point>         m_At    = {};  /**< The points, by number.  */
    List<std::uint32_t> m_Slots = {};  /**< The table: a point's number and one, 0 for
                                            an empty slot.                            */
  };

  /**
   *  A cut on a curve: where along it, and the number of the point there.
   */
  struct Cutting
  {
    std::size_t Curve = 0;   /**< The curve.          */
    double      Along = 0.0; /**< Where along it.     */
    std::size_t Node  = 0;   /**< The point's number. */
  };

  /**
   *  A stretch of a curve between two cuts, kept as a part of the region's edge.
   */
  struct Stretch
  {
    std::size_t Curve = 0;   /**< The curve.                */
    double      From  = 0.0; /**< Where along it it starts. */
    double      Onto  = 0.0; /**< Where it ends.            */
    std::size_t Start = 0;   /**< The number of its start.  */
    std::size_t End   = 0;   /**< The number of its end.    */
  };

  /**
   *  A corner's mitre or an end's square: a part of the band past what lies within the distance of
   *  the shape.
   */
  struct Wedge
  {
    List<Piece>       Rim   = {};     /**< Its edge, a loop of pieces. */
    Box               Hull  = {};     /**< A box that holds it.        */
    List<std::size_t> Edges = {};     /**< A mitre: the curve along each piece of its edge, none
                                           for the two from the corner.                       */
    bool Coned = false;               /**< Whether it is a corner's mitre, which lies between the
                                           ways from the corner to its two ends.                   */
    Point                 About = {}; /**< That corner. */
    std::array<Vector, 2> Span  = {}; /**< Those ways, the second counterclockwise from the
                                           first by no more than half a turn.          */
    bool Left = true;                 /**< Whether it lies on the left of its edge as that runs. */
  };

  /**
   *  How a numbered point is made, that it may be worked out again in two doubles a number.
   */
  enum class Making : std::uint8_t
  {
    Unknown, /**< Not told yet.                                                  */
    Hit,     /**< Where two curves' lines or circles meet.                       */
    Given,   /**< The doubles themselves are the point.                          */
    Foot,    /**< A vertex moved the distance square to a piece that ends there. */
    Held,    /**< Worked out when its curve was drawn, and kept.                 */
  };

  /**
   *  How a numbered point is made.
   */
  struct Recipe
  {
    Making      Kind   = Making::Unknown; /**< Which way.             */
    std::size_t First  = None;            /**< A curve; or the piece. */
    std::size_t Second = None;            /**< The curve it is met by; or 1 where the vertex
                                               is the piece's end, 0 its start.               */
    Point About = {};                     /**< The vertex, or the point kept. */
  };

  /**
   *  Where a numbered point comes from: an end of a curve, or where two curves meet. Of the points
   *  taken for one, the first numbered tells: the curves' ends, numbered before any meeting.
   */
  struct Birth
  {
    std::uint32_t Curve = 0; /**< The curve.                                             */
    std::uint32_t Other = 0; /**< The curve it meets; or Starts or Ends for the curve's. */
  };

  /**
   *  A curve's start, as where a point comes from.
   */
  inline constexpr std::uint32_t Starts = 0xFFFFFFFEU;

  /**
   *  A curve's end.
   */
  inline constexpr std::uint32_t Ends = 0xFFFFFFFFU;

  /**
   *  A polyline's offset as it is worked out: the shape, the curves, their cuts and what is kept.
   */
  struct Work
  {
    std::span<const Vertex> Vertices;       /**< The polyline's vertices, about its first. */
    bool                    Closed = false; /**< Whether it is closed.                     */
    Buffering               Asked  = {};    /**< The distance, the side and the filling.   */
    List<Piece>             Shape  = {};    /**< Its pieces, those of no length left out.  */
    List<Piece>             Curves = {};    /**< The curves the edge can run along, each
                                                 with the band on its left.                 */
    List<double> Ways = {};                 /**< How far to its piece's left each was moved,
                                                 0 for a circle about a vertex.             */
    List<std::size_t> Rights = {};          /**< Each piece's curve moved right, or none. */
    List<std::size_t> Lefts  = {};          /**< Each piece's curve moved left, or none.  */
    Box               All    = {};          /**< A box that holds every curve.            */
    double            Scale  = 0.0;         /**< The size of the shape and the distance.  */
    double            Reach  = 0.0;         /**< Within which two points are one.         */
    double            Pad    = 0.0;         /**< How far past a piece's box a point within
                                                 the distance of it may lie.                */
    List<Carrier>      WideShape = {};      /**< The pieces' carriers, each once asked.      */
    List<std::uint8_t> Carried   = {};      /**< Which of them are made.                     */
    List<Birth>        By        = {};      /**< Where each numbered point comes from.       */
    List<std::uint8_t> Fined     = {};      /**< Which numbered points are worked out again. */
    List<std::array<std::uint32_t, 3>> Feet = {}; /**< For each circle about a vertex, in the
                                             order they are put: the curve, and for its
                                             start and its end twice the piece that
                                             ends there, and 1 more at that piece's
                                             end.                                   */
    List<std::array<Recipe, 2>> Told = {}; /**< How each curve's start and end are made, where
                                                told as the curve was drawn.                 */
    List<Carrier>      WideCurves = {};    /**< The curves' carriers, each once asked.  */
    List<std::uint8_t> Made       = {};    /**< Which of the curves' carriers are made. */
    List<Cutting>      Cuts       = {};    /**< The cuts.                               */
    Stalls             Near       = {};    /**< The pieces' boxes grown by the pad, where
                                              there are more than a handful.             */
    List<std::size_t> Everyone = {};       /**< Every piece's number, where there are no
                                                more than a handful.                       */
    List<std::uint8_t> Even = {};          /**< Whether the shape runs on at each piece's
                                                start with no turn doubles can tell.       */
    List<Box> Hulls = {};                  /**< The pieces' boxes. */
    List<Box> Pads  = {};                  /**< The pieces' boxes grown by the pad: a point
                                                outside one is farther from that piece
                                                than the distance.                         */
    int Keeps =
        0; /**< Which side alone is offset, the right above 0 and the left below: 0 for both. */
    int Flank = 0; /**< One side of an open shape asked, not a region: 1 its right, -1 its
                        left; 0 for a region.                                              */
    List<Wedge>   Wedges    = {}; /**< The mitres and squares. */
    List<Stretch> Stretches = {}; /**< What is kept.           */
  };
} // namespace Geoxy::Belt

#endif /* GEOXY_BELT_WORK_HPP */
