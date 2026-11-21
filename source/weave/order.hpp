/***************************************************************************************************
 *  @file      source/weave/order.hpp
 *  @brief     The pieces between a curve's notches, in order along it, and the pieces of one
 *             stretch made one edge.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#ifndef GEOXY_WEAVE_ORDER_HPP
#define GEOXY_WEAVE_ORDER_HPP

#include "arranged.hpp"
#include "scratch.hpp"
#include "unions.hpp"
#include "weave/found.hpp"

#include "geoxy/curve.hpp"
#include "geoxy/domain.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

namespace Geoxy::Weave
{
  /**
   *  A piece of a curve between two cuts.
   */
  struct Piece
  {
    std::size_t                Source  = 0;  /**< The curve.                */
    Range                      Stretch = {}; /**< Its stretch of it.        */
    std::array<std::size_t, 2> Ends    = {}; /**< The vertices at its ends. */
    std::size_t                Lines   = 1;  /**< How many whole lines in a row
                                                  it runs along, from its curve. */
  };

  /**
   *  The cuts of the curves that have more than their ends', each such curve's in order along it,
   *  and whether any two at one place joined two sets of vertices.
   */
  struct Sequence
  {
    List<std::size_t> Starts = {}; /**< Where each curve's cuts start, and past the last,
                                        their end; a curve of its ends' alone has none. */
    List<Notch>  Cuts   = {};      /**< The cuts, curve by curve.                      */
    bool         Joined = false;   /**< Whether two cuts at one place joined two sets. */
    List<Couple> Joins  = {};      /**< The vertices of each two cuts at one place.    */
  };

  /**
   *  Every curve's cuts in order along it, the vertices of cuts at one place made one: a curve's
   *  ends' cuts, at 0, 0.5 and 1, in order already, and the rest put curve by curve by counting and
   *  each such curve's few sorted with them, its Bins or more put in order as Bucketed puts them.
   *  Sooner orders them wholly, no two of a curve's being at one vertex, so the order is the one
   *  sort of them all; it is asked through a lambda, which the sort inlines, as it would not a
   *  function's address.
   *
   *  @param  found The findings.
   *  @param  sets  The merged vertices.
   *  @param  kept  Whether the vertices of each two cuts at one place are kept.
   *  @return       The cuts of the curves that have more.
   */
  auto Sequenced(const Findings &found, Unions &sets, bool kept) -> Sequence;

  /**
   *  Every vertex numbered by the one its set stands for, every merge made: the least of the set,
   *  so numbered before it; or, for a vertex of a cluster split, the least of its point's.
   *
   *  @param  found   The findings.
   *  @param  sets    The merged vertices, every merge made.
   *  @param  leaders Each vertex's leader where its cluster was split, past the last vertex
   *                  elsewhere; or none.
   *  @return         Each vertex's number among those that stand for themselves.
   */
  auto Numbered(const Findings &found, Unions &sets, std::span<const std::size_t> leaders)
      -> List<std::size_t>;

  /**
   *  What tells a cluster of merged vertices whose lines cross where they do not: the curves, the
   *  vertices merged into one at another point, and the vertices of cuts at one place.
   */
  struct Tangle
  {
    std::span<const Curve>       Sources = {}; /**< The lines and arcs.            */
    std::span<const std::size_t> Moved   = {}; /**< The vertices merged elsewhere. */
    std::span<const Couple>      Joins   = {}; /**< Cuts' vertices at one place.   */
  };

  /**
   *  What untwisting found: whether a cluster is twisted, and each vertex's leader where one was
   *  split.
   */
  struct Untwist
  {
    bool              Twisted = false; /**< Whether a cluster is twisted. */
    List<std::size_t> Leaders = {};    /**< Each vertex's leader where its cluster was split,
                                            past the last vertex elsewhere; or none.   */
  };

  /**
   *  The clusters of vertices merged by distance through which two lines pass that meet nowhere in
   *  them, told and split: each vertex then stands with those at its very point, its twins and the
   *  cuts at its place along a curve, an atom, led by their least.
   *
   *  Lines through one point cross there, as their halves are put in order about it by their ways;
   *  two lines that pass through a cluster and meet at another vertex would cross twice, which no
   *  two lines do, and the halves about the two vertices are no plane's. Within an atom every two
   *  lines meet, at its point. A cluster whose every two lines meet in an atom of it has its lines
   *  cross there alone, and is kept.
   *
   *  @param  found   The findings.
   *  @param  sets    The merged vertices, every merge made.
   *  @param  strands The curves, the vertices merged elsewhere and those at one place.
   *  @return         Whether a cluster is twisted, and the leaders where split.
   */
  auto Untwisted(const Findings &found, Unions &sets, const Tangle &strands) -> Untwist;

  /**
   *  The pieces of every curve between its cuts in order, curve by curve, each curve's in order
   *  along it, their ends numbered as the arrangement's vertices; pieces whose ends merged dropped.
   *
   *  @tparam VisitT   The kind of what is done with a piece.
   *  @param  found    The findings.
   *  @param  sequence The cuts of the curves that have more than their ends'.
   *  @param  visit    What is done with each piece, given each vertex's number too.
   */
  template <typename VisitT>
  void Notched(const Findings &found, const Sequence &sequence, const VisitT &visit)
  {
    for (std::size_t curve = 0; curve < found.Count; ++curve)
    {
      const std::size_t            from  = sequence.Starts[curve];
      const std::size_t            until = sequence.Starts[curve + 1];
      const auto                   ended = Ended(found, curve);
      const std::span<const Notch> cuts =
          from < until ? std::span<const Notch>(sequence.Cuts).subspan(from, until - from)
                       : std::span<const Notch>(ended.first).first(ended.second);
      for (std::size_t index = 1; index < cuts.size(); ++index)
      {
        visit(cuts[index - 1], cuts[index]);
      }
    }
  }

  /**
   *  The pieces each shared stretch makes one: of its two curves, those inside it whose ends are
   *  the same two vertices and that the stretch maps one onto the other, so that two pieces of a
   *  circle between the same two vertices, its two halves, stay two.
   *
   *  @param  found  The findings, some stretch shared.
   *  @param  pieces The pieces.
   *  @return        The sets of pieces that are one.
   */
  auto Coincide(const Findings &found, std::span<const Piece> pieces) -> Unions;

  /**
   *  How the cuts of the curves are made pieces: the cuts in order, and each vertex's number.
   */
  struct Numbering
  {
    const Sequence          *Order = nullptr; /**< The cuts in order.    */
    const List<std::size_t> *Ids   = nullptr; /**< Each vertex's number. */
  };

  /**
   *  The pieces of the curves, each between two cuts whose vertices are two, curve by curve. Where
   *  chained, a run of lines joined at lone joints is one piece, from the first's start to the
   *  last's end: its lines, uncut, would be pieces in a row through joints no other piece reaches,
   *  which are no vertices. Its two ends are two vertices where the loops close, as a boolean's
   *  must: the start is the end of a curve before it, the one before or the last of its loop, and
   *  an end within the length of it would be in a pair with that curve, which is next to it only
   *  where it is the curve after the run, and the run's first would then be in a pair with it.
   *
   *  @tparam VisitT    The kind of what is done with a piece.
   *  @param  found     The findings.
   *  @param  numbering The cuts in order and the vertices' numbers.
   *  @param  visit     What is done with each piece.
   */
  template <typename VisitT>
  void Hewn(const Findings &found, const Numbering &numbering, const VisitT &visit)
  {
    const List<std::size_t> &ids    = *numbering.Ids;
    const auto               joined = [&found](std::size_t curve) -> bool
    {
      return curve < found.Strung.size() && found.Strung[curve] != 0;
    };
    std::size_t through = 0;
    Notched(found, *numbering.Order,
            [&](const Notch &before, const Notch &after) -> void
            {
              const std::size_t curve = before.Source;
              if (curve < through)
              {
                return;
              }
              std::size_t last = curve;
              while (joined(last))
              {
                ++last;
              }
              const bool        run   = last > curve;
              const std::size_t start = ids[before.Vertex];
              const std::size_t end   = ids[run ? found.Firsts[last + 1] - 1 : after.Vertex];
              through                 = run ? last + 1 : through;
              if (run || start != end)
              {
                visit(Piece{.Source  = curve,
                            .Stretch = {.Start = before.Along, .End = run ? 1.0 : after.Along},
                            .Ends    = {start, end},
                            .Lines   = last - curve + 1});
              }
            });
  }

  /**
   *  Whether two curves or more are lines: what a cluster needs to be twisted, so that one that
   *  none can be is not looked for, as among the arcs of a bulged polygon's corners, which merge by
   *  distance from points a few steps of a double apart.
   *
   *  @param  sources The lines and arcs.
   *  @return         Whether they are.
   */
  auto Twistable(std::span<const Curve> sources) -> bool;

  /**
   *  Each vertex's lone tell, as Arranged's Lone holds it: the two lines of a crossing, the lesser
   *  first, where no other point was numbered so; Joint where the points numbered so are two ends,
   *  of one curve and the next; Crowded elsewhere. Each crossing's lines are those of the two cuts
   *  at its point, and a crossing whose point is a curve's end is that end's twin, so that a vertex
   *  of one point that is no end is a crossing cut on both its lines.
   *
   *  @param found The findings, of lines alone and no stretch shared.
   *  @param ids   Each point's vertex.
   *  @param lone  Each vertex's tell, sized.
   */
  void Lonesome(const Findings &found, std::span<const std::size_t> ids,
                std::pmr::vector<std::array<std::uint32_t, 2>> &lone);

  /**
   *  Whether an arrangement's vertices are told lone: one made for a boolean, of lines alone, no
   *  stretch shared.
   *
   *  @param  inputs   What it is made of.
   *  @param  found    The findings.
   *  @param  arranged The arrangement, its curves laid.
   *  @return          Whether they are.
   */
  auto Solitary(const Inputs &inputs, const Findings &found, const Arranged &arranged) -> bool;
} // namespace Geoxy::Weave

#endif /* GEOXY_WEAVE_ORDER_HPP */
