/***************************************************************************************************
 *  @file      source/arranged.hpp
 *  @brief     Lines and arcs cut where they meet into the edges of a plane graph, its faces traced,
 *             and what loops, regions and booleans read from them.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_ARRANGED_HPP
#define GEOXY_ARRANGED_HPP

#include "geoxy/box.hpp"
#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"
#include "geoxy/error.hpp"
#include "geoxy/line.hpp"
#include "geoxy/loop.hpp"
#include "geoxy/point.hpp"
#include "geoxy/region.hpp"
#include "geoxy/tolerance.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <memory_resource>
#include <optional>
#include <span>
#include <vector>

namespace Geoxy
{
  /**
   *  What a curve weighs in each of two operands: how many times a region's loops run along it.
   */
  using Tally = std::array<int, 2>;

  /**
   *  The curves an arrangement is made of, run by run, each run with its tally; read where they
   *  lie, none copied.
   */
  struct Inputs
  {
    std::vector<std::span<const Curve>> Runs    = {}; /**< The runs of curves.         */
    std::vector<Tally>                  Tallies = {}; /**< Each run's tally, in order. */
    std::pmr::memory_resource          *Memory  = std::pmr::get_default_resource(); /**< What the
                                                     arrangement's lists draw on. */
    bool Chained = false;   /**< Whether a run of lines that meet nothing but each other may stand
                                 as one edge: for an arrangement only traced, wound and stitched. */
    bool Untwisted = false; /**< Whether a cluster of vertices merged by distance through which
                                 two lines pass that meet nowhere in it is split into its points,
                                 as Twisted tells of one. */
    bool Spare = false;     /**< Whether, where such a cluster is merged, the arrangement with
                                 every one split is made too, from the same cuts, as Split: for
                                 one whose clusters are not split. */
  };

  /**
   *  Curves cut where they meet: the vertices, the edges between them, and, once traced, the cycles
   *  of the faces, each half followed by its next.
   *
   *  Each edge has two halves, 2e running from its first vertex to its second and 2e + 1 back; the
   *  face of a half is the one on its left. An edge is its source's stretch alone, held in a few
   *  numbers; the pieces as curves and their boxes are made only where they are asked for, by
   *  Outfitted.
   */
  struct Arranged
  {
    /**
     *  An edge: a piece of a line or an arc, from one vertex to another.
     */
    struct Edge
    {
      std::size_t                Source  = 0;  /**< The line or arc it is a piece of.   */
      Range                      Stretch = {}; /**< Its stretch of that curve, forward. */
      std::array<std::size_t, 2> Ends    = {}; /**< The vertices at its start and end.  */
      Tally                      Weighs  = {}; /**< Its tally, along it.                */
      std::size_t                Lines   = 1;  /**< How many whole lines it runs along,
                                                    its source and those after it, each
                                                    ending where the next starts: more
                                                    than one only in a chained one.     */
    };

    /**
     *  A cycle of halves, each followed by the one just clockwise of its twin about its end: a
     *  face's boundary, or a component's outer one.
     */
    struct Cycle
    {
      std::size_t First = 0;      /**< Its first half.              */
      double      Area  = 0.0;    /**< The signed area it encloses. */
      double      Slack = 0.0;    /**< How far the area may lie from the exact
                                       one: 0 but where worked out plainly.    */
      Box         Bounds    = {}; /**< Its box, once Outfitted makes it.         */
      std::size_t Component = 0;  /**< The component of the graph it belongs to. */
    };

    /**
     *  An arrangement of nothing, its lists drawing on a memory resource: the thread's arena's for
     *  one made, used and dropped in one call, as a boolean's is, the default's for one kept.
     *
     *  @param memory What its lists draw on.
     */
    explicit Arranged(std::pmr::memory_resource *memory = std::pmr::get_default_resource())
        : Sources(memory), Points(memory), Edges(memory), Shapes(memory), Boxes(memory),
          Next(memory), Faces(memory), Cycles(memory), Outer(memory), Samples(memory), Runs(memory),
          Tallies(memory), Lone(memory)
    {
    }

    /**
     *  What Lone holds for a vertex where one curve ends and the next starts and nothing else is.
     */
    static constexpr std::uint32_t Joint = 0xFFFF'FFFEU;

    /**
     *  What Lone holds for a vertex that is neither a lone crossing nor a lone joint.
     */
    static constexpr std::uint32_t Crowded = 0xFFFF'FFFFU;

    Arranged(const Arranged &)                     = default;
    Arranged(Arranged &&) noexcept                 = default;
    auto operator=(const Arranged &) -> Arranged & = delete;
    auto operator=(Arranged &&) -> Arranged &      = delete;
    ~Arranged()                                    = default;

    std::pmr::vector<Curve>       Sources; /**< The lines and arcs.                */
    std::pmr::vector<Point>       Points;  /**< The vertices.                      */
    std::pmr::vector<Edge>        Edges;   /**< The edges.                         */
    std::pmr::vector<Curve>       Shapes;  /**< Each edge's piece, once Outfitted. */
    std::pmr::vector<Box>         Boxes;   /**< Each edge's box, once Outfitted.   */
    std::pmr::vector<std::size_t> Next;    /**< Each half's next in its cycle: the half
                                                just clockwise of its twin about the
                                                vertex it ends at.                     */
    std::pmr::vector<std::size_t> Faces;   /**< Each half's cycle.             */
    std::pmr::vector<Cycle>       Cycles;  /**< The cycles.                    */
    std::pmr::vector<std::size_t> Outer;   /**< Each component's outer cycle.  */
    std::pmr::vector<std::size_t> Samples; /**< Each component's least vertex. */
    std::pmr::vector<std::size_t> Runs;    /**< Each run's first curve, and past the last,
                                                how many curves there are.             */
    std::pmr::vector<Tally>                        Tallies; /**< Each curve's tally. */
    std::pmr::vector<std::array<std::uint32_t, 2>> Lone;    /**< Each vertex's two lines, the lesser
                                                       first, where they cross there and no other
                                                       point was merged into it; Joint twice where
                                                       one curve ends and the next starts there and
                                                       nothing else; Crowded twice elsewhere. Told
                                                       only where made for a boolean of lines, no
                                                       stretch shared; none elsewhere.          */
    Tolerance            Within = {};                       /**< The tolerance, made fit.  */
    std::optional<Error> Fault  = {};                       /**< Why it could not be made. */
    bool Chained = false; /**< Whether made for a boolean alone: a run of lines may be one edge,
                               and a cycle's area is told only as far as which is least. */
    bool Twisted = false; /**< Whether two lines pass through vertices merged by distance into
                               one and meet nowhere in it: lines through one point cross there,
                               so they would cross twice, or where they do not. Told, and such
                               clusters split, only where the inputs asked. */
    bool Banded = false;  /**< Whether the halves about a vertex whose curves stay within twice
                               the length of each other until they cross, as Intersect takes for one
                               hit, are ordered by where they part, not where they leave. */
    std::shared_ptr<Arranged> Split = {}; /**< The same curves with every such cluster split,
                                               where the inputs asked for it to spare and one
                                               is merged; none elsewhere. */
  };

  /**
   *  Curves cut where they meet into edges: lines, arcs of up to a turn and polylines' segments;
   *  every crossing and shared stretch from a drawing of them met all at once, cuts in a shared
   *  stretch made in both its curves, vertices merged in order, and pieces of one shared stretch
   *  made one edge.
   *
   *  @param  inputs    The curves and their tallies.
   *  @param  tolerance Within which two points are one.
   *  @return           The vertices and edges, not traced; with a fault where a curve is of another
   *                    kind, not finite, or a meeting's budget ran out.
   */
  auto Arrange(const Inputs &inputs, const Tolerance &tolerance) -> Arranged;

  /**
   *  Which faces a boolean keeps, by whether each operand fills them.
   */
  using Keeps = auto (*)(bool first, bool second) -> bool;

  /**
   *  What reading a boolean's loops from windings ferried along its curves came to.
   */
  enum class Ferry : std::uint8_t
  {
    Done,    /**< The loops were read.                                                 */
    Twisted, /**< The windings ferried disagree at a vertex: the faces are no plane's. */
    Unfit,   /**< Not read: no lone tells, more than one component, or the outer side
                  of the arrangement not told alike from its four farthest points.    */
  };

  /**
   *  The loops that bound the faces a boolean keeps, read with no face traced: each run's pieces
   *  walked in turn, the windings on their left ferried from one to the next across each vertex,
   *  changed by the tallies of the halves that leave it on that side; then only the pieces whose
   *  sides a boolean keeps one of are ordered about their vertices and stitched.
   *
   *  At a lone crossing of two lines the other line's halves lie one on each side, as the sign of
   *  their ways' cross product tells, where they are farther apart than twice the angle; at a lone
   *  joint nothing changes; about any other vertex the halves are ordered as Trace orders them.
   *  Each run's windings are known up to a constant, and one more for each such other vertex they
   *  pass; every vertex ties the constants of the two sides of each of its sectors, and a tie that
   *  disagrees with those before means the faces are no plane's: the halves about a vertex are
   *  ordered as no drawing leaves them. Where every tie agrees, the windings are those Windings
   *  gives the traced faces, the outer side of the one component its windings' zero, as told at its
   *  farthest point each way along the axes; four that disagree leave it unread.
   *
   *  @param  arranged The arrangement, made for a boolean, not traced.
   *  @param  keeps    The rule the faces are kept by.
   *  @param  region   The loops, where read.
   *  @return          What it came to.
   */
  auto Ferried(const Arranged &arranged, Keeps keeps, Region &region) -> Ferry;

  /**
   *  A closed run of halves as a loop, the halves that run on along one curve joined: started after
   *  a half that does not run on into its next, so that no run is cut by where it starts. A half
   *  along a run of lines stands for its lines' halves, which run on into none; it never comes
   *  first, as the piece before it, of a lower number, runs into it.
   *
   *  @param  arranged The arrangement.
   *  @param  halves   The halves, in order.
   *  @return          The loop.
   */
  auto Sewn(const Arranged &arranged, std::span<std::size_t> halves) -> Loop;

  /**
   *  A point an edge passes from its start: its start's vertex at 0, its end's after its last line,
   *  and between them where each of its lines ends, as At gives it, the vertex that stands there,
   *  which merges with nothing but the next line's start.
   *
   *  @param  arranged The arrangement.
   *  @param  edge     The edge.
   *  @param  joint    How many of its lines lie before the point, up to all of them.
   *  @return          The point.
   */
  inline auto Beaded(const Arranged &arranged, std::size_t edge, std::size_t joint) -> Point
  {
    const Arranged::Edge &piece = arranged.Edges[edge];
    if (joint == 0 || joint == piece.Lines)
    {
      return arranged.Points[piece.Ends.at(joint == 0 ? 0 : 1)];
    }
    const Line &line = std::get<Line>(arranged.Sources[piece.Source + joint - 1]);
    return Lerp(line.Start, line.End, 1.0);
  }

  /**
   *  An edge's piece of its line or arc as a curve of its own: a line between the vertices it
   *  joins, an arc's stretch of the same circle.
   *
   *  @param  arranged The arrangement.
   *  @param  edge     The edge's index.
   *  @return          The curve.
   */
  auto Carved(const Arranged &arranged, std::size_t edge) -> Curve;

  /**
   *  Every edge's piece made a curve of its own and boxed, and, once traced, each cycle's box
   *  joined from its edges' boxes: what picking a point and bridging gaps ask of an arrangement,
   *  and booleans do not.
   *
   *  @param arranged The arrangement, pruned where it is to be.
   */
  void Outfitted(Arranged &arranged);

  /**
   *  The loose ends of an arrangement not traced: the vertices only one edge reaches.
   *
   *  @param  arranged The arrangement.
   *  @return          The vertices, in order.
   */
  auto Frayed(const Arranged &arranged) -> std::vector<std::size_t>;

  /**
   *  Every edge with a loose end taken away, and every edge it leaves loose, until none is.
   *
   *  @param arranged The arrangement, not traced nor outfitted.
   */
  void Prune(Arranged &arranged);

  /**
   *  The halves of an arrangement's edges ordered about each vertex, and its cycles traced.
   *
   *  @param arranged The arrangement.
   */
  void Trace(Arranged &arranged);

  /**
   *  Whether a traced arrangement's faces are a plane's: in each component the vertices less the
   *  edges and with the cycles, the outer one with them, come to 2. Where they come to less, the
   *  halves about some vertex are not in the order a drawing in the plane would leave them in, and
   *  windings carried across its edges are not those of its faces.
   *
   *  @param  arranged The traced arrangement.
   *  @return          Whether they are.
   */
  auto Planar(const Arranged &arranged) -> bool;

  /**
   *  Each two pieces between the same two vertices that cross where they have no vertex made one
   *  edge, their tallies summed along the first, and the arrangement traced again: two pieces of
   *  lines, drawn between the same two points, one segment; a piece of an arc and another, halfway
   *  along within four lengths of each other, that bound no face of their own.
   *
   *  Where a cluster merged by distance holds a crossing, a touch or an end that lies farther than
   *  the length from a line or an arc through it, the halves about it, ordered by their tangents
   *  there, may have two curves cross both at it and at another vertex, though they meet once: the
   *  pieces between the two then cross, and the faces are no plane's. Two such pieces lie within
   *  the tolerance of each other at both ends, and as one they cross nowhere.
   *
   *  @param  arranged The traced arrangement, not outfitted.
   *  @return          Whether any two were made one.
   */
  auto Melded(Arranged &arranged) -> bool;

  /**
   *  A traced arrangement whose faces are a plane's, from one: itself where its are; else the one
   *  with every twisted cluster split, spared with it, traced, where that one's are; else itself
   *  with the pieces that cross between two vertices made one, as Melded makes them, where that
   *  one's are; else the same traced banded, before and after, where either's are; else itself.
   *
   *  @param  arranged The traced arrangement.
   *  @param  pruned   Whether it was pruned, and so is the split one before it is traced.
   *  @return          The one chosen, traced.
   */
  auto Mended(Arranged arranged, bool pruned) -> Arranged;

  /**
   *  Each half of a cycle visited, from one of them, each followed by its next until it comes
   *  round.
   *
   *  @tparam VisitT   The kind of the visit.
   *  @param  arranged The traced arrangement.
   *  @param  first    The half it starts from.
   *  @param  visit    What is done with each half.
   */
  template <typename VisitT>
  void Toured(const Arranged &arranged, std::size_t first, const VisitT &visit)
  {
    visit(first);
    for (std::size_t half = arranged.Next[first]; half != first; half = arranged.Next[half])
    {
      visit(half);
    }
  }

  /**
   *  The half just clockwise of a half about the vertex it leaves: the next of its twin.
   *
   *  @param  arranged The traced arrangement.
   *  @param  half     The half.
   *  @return          The half before it in its vertex's counterclockwise order, the last for the
   *                   first.
   */
  inline auto Clockwise(const Arranged &arranged, std::size_t half) -> std::size_t
  {
    return arranged.Next[half ^ 1U];
  }

  /**
   *  How many times each operand's loops wind about each face of a traced arrangement.
   *
   *  @param  arranged The arrangement.
   *  @return          Each cycle's windings, the face on the left of its halves.
   */
  auto Windings(const Arranged &arranged) -> std::vector<Tally>;

  /**
   *  The loops that bound the faces kept: every half with a kept face on its left and a dropped one
   *  on its right, each followed by the first such half clockwise from its twin.
   *
   *  @param  arranged The traced arrangement.
   *  @param  kept     Whether each cycle's face is kept.
   *  @return          The loops, pieces of one curve that follow each other joined.
   */
  auto Stitched(const Arranged &arranged, const std::vector<bool> &kept) -> Region;

  /**
   *  A cycle as a loop.
   *
   *  @param  arranged The traced arrangement.
   *  @param  cycle    The cycle.
   *  @return          Its loop, pieces of one curve that follow each other joined.
   */
  auto Ringed(const Arranged &arranged, std::size_t cycle) -> Loop;

  /**
   *  The cycle of least area that runs counterclockwise and winds about a point, among those of
   *  every component but one.
   *
   *  @param  arranged The traced and outfitted arrangement.
   *  @param  point    The point, farther than the tolerance from every edge.
   *  @param  skipped  The component left out, or one past the last for none.
   *  @return          The cycle; none where none winds about it.
   */
  auto Enclosing(const Arranged &arranged, Point point, std::size_t skipped)
      -> std::optional<std::size_t>;

  /**
   *  Whether a point is within an arrangement's tolerance of an edge.
   *
   *  @param  arranged The outfitted arrangement.
   *  @param  point    The point.
   *  @return          Whether it is.
   */
  auto Grazes(const Arranged &arranged, Point point) -> bool;
} // namespace Geoxy

#endif /* GEOXY_ARRANGED_HPP */
