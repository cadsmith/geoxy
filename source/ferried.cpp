/***************************************************************************************************
 *  @file      source/ferried.cpp
 *  @brief     The loops that bound the faces a boolean keeps, read with no face traced: the
 *             windings ferried along each run's pieces from one vertex to the next.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "arranged.hpp"
#include "scratch.hpp"
#include "traced.hpp"

#include "geoxy/line.hpp"
#include "geoxy/point.hpp"
#include "geoxy/region.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <optional>
#include <span>
#include <tuple>
#include <utility>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  One tally less another, operand by operand.
     *
     *  @param  lhs The tally.
     *  @param  rhs The one taken away.
     *  @return     The difference.
     */
    auto Without(const Tally &lhs, const Tally &rhs) -> Tally
    {
      return {lhs[0] - rhs[0], lhs[1] - rhs[1]};
    }

    /**
     *  Two tallies summed, operand by operand.
     *
     *  @param  lhs One tally.
     *  @param  rhs The other.
     *  @return     The sum.
     */
    auto Stacked(const Tally &lhs, const Tally &rhs) -> Tally
    {
      return {lhs[0] + rhs[0], lhs[1] + rhs[1]};
    }

    /**
     *  Unknown windings, each known as another's plus a difference: sets of them joined one tie at
     *  a time, each set led by its least, each member holding what it is more than its parent.
     */
    class Strata
    {
    public:
      /**
       *  A new unknown, a set of its own.
       *
       *  @return Its index.
       */
      auto Add() -> std::uint32_t
      {
        const auto index = static_cast<std::uint32_t>(m_Parents.size());
        m_Parents.push_back(index);
        m_Above.push_back({});
        return index;
      }

      /**
       *  How many unknowns there are.
       *
       *  @return The count.
       */
      [[nodiscard]] auto Size() const -> std::size_t { return m_Parents.size(); }

      /**
       *  An unknown's set's leader, and what the unknown is more than it; the path to it made
       *  straight.
       *
       *  @param  member The unknown.
       *  @return        The leader and the difference.
       */
      auto Root(std::uint32_t member) -> std::pair<std::uint32_t, Tally>
      {
        std::uint32_t root  = member;
        Tally         above = {};
        while (m_Parents[root] != root)
        {
          above = Stacked(above, m_Above[root]);
          root  = m_Parents[root];
        }
        Tally rest = above;
        while (m_Parents[member] != root)
        {
          const std::uint32_t parent = m_Parents[member];
          const Tally         own    = m_Above[member];
          m_Parents[member]          = root;
          m_Above[member]            = rest;
          rest                       = Without(rest, own);
          member                     = parent;
        }
        return {root, above};
      }

      /**
       *  That one unknown is another's plus a difference: their sets joined, the greater leader
       *  under the lesser.
       *
       *  @param  one        The unknown.
       *  @param  other      The other.
       *  @param  difference What the first is more than the second.
       *  @return            Whether that agrees with what was known.
       */
      auto Tie(std::uint32_t one, std::uint32_t other, const Tally &difference) -> bool
      {
        const auto [first, above]  = Root(one);
        const auto [second, below] = Root(other);
        const Tally gap            = Stacked(Without(difference, above), below);
        if (first == second)
        {
          return gap == Tally{};
        }
        const bool lower                  = first < second;
        m_Parents[lower ? second : first] = lower ? first : second;
        m_Above[lower ? second : first]   = lower ? Without({}, gap) : gap;
        return true;
      }

    private:
      List<std::uint32_t> m_Parents = {}; /**< Each unknown's parent.             */
      List<Tally>         m_Above   = {}; /**< What each is more than its parent. */
    };

    /**
     *  A piece's winding as ferried: its stretch's unknown, and the winding on its left more than
     *  that unknown.
     */
    struct Hauled
    {
      std::uint32_t Stratum = 0;  /**< The unknown.         */
      Tally         Left    = {}; /**< The winding over it. */
    };

    /**
     *  A junction at a vertex neither a lone crossing nor a lone joint: the vertex, the piece that
     *  ends there and the one that starts there, next in its run.
     */
    struct Snag
    {
      std::uint32_t Vertex = 0; /**< The vertex.    */
      std::uint32_t In     = 0; /**< The piece in.  */
      std::uint32_t Out    = 0; /**< The piece out. */
    };

    /**
     *  A vertex as the first run to pass it left it: at a lone crossing the unknown and winding
     *  left of the piece out; at any other vertex the unknown alone, telling it reached; or none.
     */
    struct Rendezvous
    {
      std::uint32_t Stratum = Arranged::Crowded; /**< The unknown, or Crowded for none. */
      Tally         Left    = {};                /**< The winding over it.              */
    };

    /**
     *  The farthest of some points one way along an axis: a joint of a run of lines, no vertex,
     *  with the run's edge and how many of its lines lie before it; or a lone joint's vertex, with
     *  the piece in and the piece out.
     */
    struct Outmost
    {
      Point       At    = {};    /**< The point.                          */
      std::size_t In    = 0;     /**< The run, or the piece in.           */
      std::size_t Out   = 0;     /**< The run, or the piece out.          */
      std::size_t Joint = 0;     /**< The run's lines before it, or none. */
      bool        Found = false; /**< Whether any was found.              */
    };

    /**
     *  The four ways along the axes the outer side is told from: least x, most x, least y, most y.
     */
    constexpr std::array<Vector, 4> Rimward = {
        Vector{.X = -1.0, .Y =  0.0},
        Vector{.X =  1.0, .Y =  0.0},
        Vector{.X =  0.0, .Y = -1.0},
        Vector{.X = 0.0, .Y = 1.0}};

    /**
     *  Whether one point lies farther than another one of the four ways.
     *
     *  @param  lhs The point.
     *  @param  rhs The other.
     *  @param  way The way's index in Rimward.
     *  @return     Whether it does.
     */
    auto Farther(Point lhs, Point rhs, std::size_t way) -> bool
    {
      const double one   = way < 2 ? lhs.X : lhs.Y;
      const double other = way < 2 ? rhs.X : rhs.Y;
      return way % 2 == 0 ? one < other : one > other;
    }

    /**
     *  Whether a way lies strictly between two others, counterclockwise from the first: as Order
     *  turns them, from the x-axis's negative end.
     *
     *  @param  from  The first.
     *  @param  until The last.
     *  @param  way   The way.
     *  @return       Whether it does.
     */
    auto Amid(Vector from, Vector until, Vector way) -> bool
    {
      const bool after  = Preceding(from, way) < 0;
      const bool before = Preceding(way, until) < 0;
      return Preceding(from, until) < 0 ? after && before : after || before;
    }

    /**
     *  What ferrying windings along an arrangement's runs works out: each piece's winding, each
     *  lone crossing's first passage, the snags, the unknowns, and whether a tie disagreed.
     */
    struct Ferrying
    {
      List<Hauled>           Hauls    = {}; /**< Each piece's winding.            */
      List<Rendezvous>       Crossed  = {}; /**< Each vertex's first passage.     */
      List<Snag>             Snags    = {}; /**< The junctions at other vertices. */
      Strata                 Unknowns = {}; /**< The unknowns.                    */
      std::array<Outmost, 4> Joints   = {}; /**< The farthest joints of runs.     */
      std::array<Outmost, 4> Corners  = {}; /**< The farthest lone joints.        */
      List<std::size_t>      Around   = {}; /**< The halves about each snagged vertex,
                                                 counterclockwise, vertex by vertex. */
      List<std::size_t> Fans = {};          /**< Where each snagged vertex's start in
                                                 Around, and past the last, its end. */
      List<std::uint32_t> Snagged = {};     /**< The snagged vertices, in order. */
      bool                Broken  = false;  /**< Whether a tie disagreed.        */
    };

    /**
     *  The unknown and the winding on the left of a half, or on its right.
     *
     *  @param  arranged The arrangement.
     *  @param  ferrying The windings ferried.
     *  @param  half     The half, and 1 for its right.
     *  @return          The unknown and the winding over it.
     */
    auto Flanked(const Arranged &arranged, const Ferrying &ferrying,
                 const std::array<std::size_t, 2> &half) -> std::pair<std::uint32_t, Tally>
    {
      const Hauled &haul  = ferrying.Hauls[half[0] / 2];
      const bool    right = (half[0] % 2) != half[1];
      return {haul.Stratum,
              right ? Without(haul.Left, arranged.Edges[half[0] / 2].Weighs) : haul.Left};
    }

    /**
     *  Two sides tied: the winding on one is that on the other.
     *
     *  @param ferrying The windings ferried.
     *  @param one      One side's unknown and winding over it.
     *  @param other    The other's.
     */
    void Tied(Ferrying &ferrying, const std::pair<std::uint32_t, Tally> &one,
              const std::pair<std::uint32_t, Tally> &other)
    {
      if (!ferrying.Unknowns.Tie(one.first, other.first, Without(other.second, one.second)))
      {
        ferrying.Broken = true;
      }
    }

    /**
     *  Whether a line's other at a lone crossing leaves it on its left going forward, as the exact
     *  sign of their ways' cross product tells. Not told where the ways lie within twice the angle
     *  of each other or of opposite, where Fanned might order the halves by number: worked out
     *  alike from either line, the lesser's way first.
     *
     *  @param  arranged The arrangement.
     *  @param  line     The line.
     *  @param  other    The other line.
     *  @return          Whether it does, or none.
     */
    auto Leftward(const Arranged &arranged, std::size_t line, std::size_t other)
        -> std::optional<bool>
    {
      const Line  &low    = std::get<Line>(arranged.Sources[std::min(line, other)]);
      const Line  &high   = std::get<Line>(arranged.Sources[std::max(line, other)]);
      const Vector first  = low.End - low.Start;
      const Vector second = high.End - high.Start;
      const double cross  = Cross(first, second);
      const double angle  = 2.0 * std::max(arranged.Within.Radians, LeastAngle);
      if (!(std::abs(cross) > angle * (Length(first) * Length(second))))
      {
        return std::nullopt;
      }
      return (cross > 0.0) == (line < other);
    }

    /**
     *  A line's passage of a lone crossing: the unknown and the winding left of its piece out,
     *  whether the other line leaves it on its left, and the two lines' tallies.
     */
    struct Passage
    {
      std::uint32_t Stratum  = 0;     /**< The unknown.                   */
      Tally         Left     = {};    /**< The winding over it.           */
      bool          Leftward = false; /**< Whether the other leaves left. */
      Tally         Own      = {};    /**< This line's tally.             */
      Tally         Other    = {};    /**< The other line's tally.        */
    };

    /**
     *  A line's passage of a lone crossing kept where it is the first of the two, else tied to the
     *  first's: the sector after this one's piece out, counterclockwise, is the one before the
     *  other's piece out where the other leaves on this one's left, its right side, else the one
     *  after the other's piece in, the left side of that piece.
     *
     *  @param ferrying The windings ferried.
     *  @param vertex   The crossing's vertex.
     *  @param passage  The passage.
     */
    void Passed(Ferrying &ferrying, std::size_t vertex, const Passage &passage)
    {
      Rendezvous &first = ferrying.Crossed[vertex];
      if (first.Stratum == Arranged::Crowded)
      {
        first = {.Stratum = passage.Stratum, .Left = passage.Left};
        return;
      }
      const Tally across =
          passage.Leftward ? Without(first.Left, passage.Other) : Stacked(first.Left, passage.Own);
      Tied(ferrying, {passage.Stratum, passage.Left}, {first.Stratum, across});
    }

    /**
     *  The joints of a run of lines kept where farther than those kept each way along the axes.
     *
     *  @param arranged The arrangement.
     *  @param edge     The run's edge.
     *  @param ferrying The windings ferried, their farthest joints kept.
     */
    void Jointed(const Arranged &arranged, std::size_t edge, Ferrying &ferrying)
    {
      for (std::size_t joint = 1; joint < arranged.Edges[edge].Lines; ++joint)
      {
        const Point point = Beaded(arranged, edge, joint);
        for (std::size_t way = 0; way < Rimward.size(); ++way)
        {
          Outmost &kept = ferrying.Joints.at(way);
          if (!kept.Found || Farther(point, kept.At, way))
          {
            kept = {.At = point, .In = edge, .Out = edge, .Joint = joint, .Found = true};
          }
        }
      }
    }

    /**
     *  A lone joint's vertex kept where farther than those kept each way along the axes, with the
     *  piece in and the piece out.
     *
     *  @param arranged The arrangement.
     *  @param pieces   The piece in and the piece out.
     *  @param ferrying The windings ferried, their farthest lone joints kept.
     */
    void Cornered(const Arranged &arranged, const std::array<std::size_t, 2> &pieces,
                  Ferrying &ferrying)
    {
      const Point point = arranged.Points[arranged.Edges[pieces[0]].Ends[1]];
      for (std::size_t way = 0; way < Rimward.size(); ++way)
      {
        Outmost &kept = ferrying.Corners.at(way);
        if (!kept.Found || Farther(point, kept.At, way))
        {
          kept = {.At = point, .In = pieces[0], .Out = pieces[1], .Joint = 0, .Found = true};
        }
      }
    }

    /**
     *  A run's passage through the junction of one of its pieces and the next: its first piece, the
     *  two pieces, the unknown the windings are over and the winding left of the piece in.
     */
    struct Leg
    {
      std::size_t   First   = 0;  /**< The run's first piece. */
      std::size_t   In      = 0;  /**< The piece in.          */
      std::size_t   Out     = 0;  /**< The piece out.         */
      std::uint32_t Stratum = 0;  /**< The unknown.           */
      Tally         Left    = {}; /**< The winding over it.   */
    };

    /**
     *  How a run passes a junction: at a lone joint, where one curve ends and the next starts, with
     *  no change, as false; at a lone crossing, whether the other line leaves on its left, where
     *  told; else not told, a snag.
     *
     *  @param  arranged The arrangement.
     *  @param  leg      The passage.
     *  @return          Whether the other line leaves on the left, or none.
     */
    auto Passing(const Arranged &arranged, const Leg &leg) -> std::optional<bool>
    {
      const Arranged::Edge               &piece = arranged.Edges[leg.In];
      const std::array<std::uint32_t, 2> &lone  = arranged.Lone[piece.Ends[1]];
      if (piece.Source != arranged.Edges[leg.Out].Source)
      {
        return lone[0] == Arranged::Joint ? std::optional<bool>{false} : std::nullopt;
      }
      if (lone[0] >= Arranged::Joint)
      {
        return std::nullopt;
      }
      return Leftward(arranged, piece.Source, lone[0] == piece.Source ? lone[1] : lone[0]);
    }

    /**
     *  A junction passed: a lone joint kept where farthest, the vertex marked reached, and the
     *  winding ferried on; at a lone crossing changed by the other line's tally and the passage
     *  kept or tied; at a snag, the snag kept and a new unknown begun; the last piece tied back to
     *  the first.
     *
     *  @param arranged The arrangement.
     *  @param leg      The passage, its unknown and winding moved on.
     *  @param ferrying The windings ferried.
     */
    void Stepped(const Arranged &arranged, Leg &leg, Ferrying &ferrying)
    {
      const Arranged::Edge     &piece    = arranged.Edges[leg.In];
      const std::size_t         vertex   = piece.Ends[1];
      const bool                through  = piece.Source == arranged.Edges[leg.Out].Source;
      const std::optional<bool> leftward = Passing(arranged, leg);
      if (!through && leftward)
      {
        Cornered(arranged, {leg.In, leg.Out}, ferrying);
      }
      if (!through || !leftward)
      {
        ferrying.Crossed[vertex].Stratum = leg.Stratum;
      }
      if (!leftward)
      {
        ferrying.Snags.push_back({.Vertex = static_cast<std::uint32_t>(vertex),
                                  .In     = static_cast<std::uint32_t>(leg.In),
                                  .Out    = static_cast<std::uint32_t>(leg.Out)});
        if (leg.Out != leg.First)
        {
          leg.Stratum             = ferrying.Unknowns.Add();
          leg.Left                = {};
          ferrying.Hauls[leg.Out] = {.Stratum = leg.Stratum, .Left = {}};
        }
        return;
      }
      if (through)
      {
        const std::array<std::uint32_t, 2> &lone = arranged.Lone[vertex];
        const Tally tally = arranged.Tallies[lone[0] == piece.Source ? lone[1] : lone[0]];
        leg.Left          = Without(leg.Left, *leftward ? tally : Without({}, tally));
        Passed(ferrying, vertex,
               {.Stratum  = leg.Stratum,
                .Left     = leg.Left,
                .Leftward = *leftward,
                .Own      = arranged.Tallies[piece.Source],
                .Other    = tally});
      }
      if (leg.Out == leg.First)
      {
        Tied(ferrying, {ferrying.Hauls[leg.First].Stratum, {}}, {leg.Stratum, leg.Left});
        return;
      }
      ferrying.Hauls[leg.Out] = {.Stratum = leg.Stratum, .Left = leg.Left};
    }

    /**
     *  The windings ferried along one run's pieces, from the first: a new unknown at the first and
     *  after each snag, each junction passed in turn.
     *
     *  @param  arranged The arrangement.
     *  @param  pieces   The run's first piece and one past its last.
     *  @param  ferrying The windings ferried.
     *  @return          Whether each piece ends where the next starts.
     */
    auto Haul(const Arranged &arranged, const std::array<std::size_t, 2> &pieces,
              Ferrying &ferrying) -> bool
    {
      const auto [first, last] = pieces;
      Leg leg                  = {.First = first, .Stratum = ferrying.Unknowns.Add()};
      ferrying.Hauls[first]    = {.Stratum = leg.Stratum, .Left = {}};
      for (leg.In = first; leg.In < last; ++leg.In)
      {
        leg.Out = leg.In + 1 < last ? leg.In + 1 : first;
        if (arranged.Edges[leg.In].Ends[1] != arranged.Edges[leg.Out].Ends[0])
        {
          return false;
        }
        Jointed(arranged, leg.In, ferrying);
        Stepped(arranged, leg, ferrying);
      }
      return true;
    }

    /**
     *  The halves about each snagged vertex put in order as Thronged orders them, each sector's two
     *  sides tied: the winding left of each half is that right of the next counterclockwise.
     *
     *  @param arranged The arrangement.
     *  @param ferrying The windings ferried, their snags gathered.
     */
    void Snagging(const Arranged &arranged, Ferrying &ferrying)
    {
      List<Snag> &snags = ferrying.Snags;
      std::ranges::sort(
          snags, [](const Snag &lhs, const Snag &rhs) -> bool
          { return lhs.Vertex < rhs.Vertex || (lhs.Vertex == rhs.Vertex && lhs.Out < rhs.Out); });
      List<Spoke> spokes;
      for (std::size_t start = 0; start < snags.size();)
      {
        const std::uint32_t vertex = snags[start].Vertex;
        const std::size_t   from   = ferrying.Around.size();
        for (; start < snags.size() && snags[start].Vertex == vertex; ++start)
        {
          ferrying.Around.push_back((2 * std::size_t{snags[start].In}) + 1);
          ferrying.Around.push_back(2 * std::size_t{snags[start].Out});
        }
        const std::span<std::size_t> fan = std::span<std::size_t>(ferrying.Around).subspan(from);
        Ordered(arranged, fan, spokes);
        for (std::size_t place = 0; place < fan.size(); ++place)
        {
          Tied(ferrying, Flanked(arranged, ferrying, {fan[place], 0}),
               Flanked(arranged, ferrying, {fan[(place + 1) % fan.size()], 1}));
        }
        ferrying.Fans.push_back(from);
        ferrying.Snagged.push_back(vertex);
      }
      ferrying.Fans.push_back(ferrying.Around.size());
    }

    /**
     *  The halves about a snagged vertex, counterclockwise, or none where it is not snagged.
     *
     *  @param  ferrying The windings ferried, their snags ordered.
     *  @param  vertex   The vertex.
     *  @return          The halves.
     */
    auto Snarled(const Ferrying &ferrying, std::size_t vertex) -> std::span<const std::size_t>
    {
      const auto found = std::ranges::lower_bound(ferrying.Snagged, vertex);
      if (found == ferrying.Snagged.end() || *found != vertex)
      {
        return {};
      }
      const auto slot = static_cast<std::size_t>(found - ferrying.Snagged.begin());
      return std::span<const std::size_t>(ferrying.Around)
          .subspan(ferrying.Fans[slot], ferrying.Fans[slot + 1] - ferrying.Fans[slot]);
    }

    /**
     *  The way a half of a line, or of a run of lines, leaves its vertex, as Spoked tells it.
     *
     *  @param  arranged The arrangement, of lines.
     *  @param  half     The half.
     *  @return          The way.
     */
    auto Outbound(const Arranged &arranged, std::size_t half) -> Vector
    {
      const Arranged::Edge &edge = arranged.Edges[half / 2];
      const bool            back = half % 2 == 1;
      const Line           &line =
          std::get<Line>(arranged.Sources[edge.Source + (back ? edge.Lines - 1 : 0)]);
      return back ? line.Start - line.End : line.End - line.Start;
    }

    /**
     *  The side of the sector about a vertex that one of the four ways points into, where just one
     *  does: the left of the half before it, counterclockwise.
     *
     *  @param  arranged The arrangement.
     *  @param  ferrying The windings ferried.
     *  @param  fan      The halves about the vertex, counterclockwise, and the way's index.
     *  @return          The side's unknown and winding over it, or none.
     */
    auto Outside(const Arranged &arranged, const Ferrying &ferrying,
                 const std::pair<std::span<const std::size_t>, std::size_t> &fan)
        -> std::optional<std::pair<std::uint32_t, Tally>>
    {
      const auto [halves, way] = fan;
      std::optional<std::pair<std::uint32_t, Tally>> side;
      std::size_t                                    found = 0;
      for (std::size_t place = 0; place < halves.size(); ++place)
      {
        const Vector from  = Outbound(arranged, halves[place]);
        const Vector until = Outbound(arranged, halves[(place + 1) % halves.size()]);
        if (Amid(from, until, Rimward.at(way)))
        {
          side = Flanked(arranged, ferrying, {halves[place], 0});
          ++found;
        }
      }
      return found == 1 ? side : std::nullopt;
    }

    /**
     *  The side the outer face lies on at the farthest point one way along an axis: a run's joint
     *  where one lies farther than every vertex, the left of the run where the way points between
     *  its line after and its line before, counterclockwise, else its right; else the farthest lone
     *  joint where none lies farther, the sector it points into; else a snagged vertex's.
     *
     *  @param  arranged The arrangement.
     *  @param  ferrying The windings ferried.
     *  @param  rim      The farthest vertex that way, and the way's index.
     *  @return          The side's unknown and winding over it, or none told.
     */
    auto Outer(const Arranged &arranged, const Ferrying &ferrying,
               const std::array<std::size_t, 2> &rim)
        -> std::optional<std::pair<std::uint32_t, Tally>>
    {
      const auto [vertex, way] = rim;
      const Point    farthest  = arranged.Points[vertex];
      const Outmost &joint     = ferrying.Joints.at(way);
      if (joint.Found && Farther(joint.At, farthest, way))
      {
        const std::size_t source = arranged.Edges[joint.In].Source + joint.Joint;
        const Line       &after  = std::get<Line>(arranged.Sources[source]);
        const Line       &before = std::get<Line>(arranged.Sources[source - 1]);
        const bool left = Amid(after.End - after.Start, before.Start - before.End, Rimward.at(way));
        return Flanked(arranged, ferrying, {2 * joint.In, left ? 0U : 1U});
      }
      const Outmost &corner = ferrying.Corners.at(way);
      if (corner.Found && !Farther(farthest, corner.At, way))
      {
        const std::array<std::size_t, 2> halves = {(2 * corner.In) + 1, 2 * corner.Out};
        return Outside(arranged, ferrying, {halves, way});
      }
      return Outside(arranged, ferrying, {Snarled(ferrying, vertex), way});
    }

    /**
     *  Each unknown's value, where the arrangement is one component and its outer face is told
     *  alike from its farthest points the four ways along the axes, that face's windings 0.
     *
     *  @param  arranged The arrangement.
     *  @param  ferrying The windings ferried, every tie made.
     *  @param  lifts    Each unknown's value, sized.
     *  @return          Whether told.
     */
    auto Anchoring(const Arranged &arranged, Ferrying &ferrying, List<Tally> &lifts) -> bool
    {
      const std::size_t   count = ferrying.Unknowns.Size();
      List<std::uint32_t> roots(count);
      std::size_t         components = 0;
      for (std::size_t index = 0; index < count; ++index)
      {
        std::tie(roots[index], lifts[index]) =
            ferrying.Unknowns.Root(static_cast<std::uint32_t>(index));
        components += static_cast<std::size_t>(roots[index] == index);
      }
      if (components != 1)
      {
        return false;
      }
      /**
       *  The farthest vertices each way that an edge reaches: each was passed by a run, as a lone
       *  crossing, a joint or a snag; a loop shrunk to a point leaves one no edge reaches.
       */
      std::array<std::size_t, 4> rims = {};
      rims.fill(arranged.Edges.front().Ends[0]);
      for (std::size_t vertex = 0; vertex < arranged.Points.size(); ++vertex)
      {
        if (ferrying.Crossed[vertex].Stratum == Arranged::Crowded)
        {
          continue;
        }
        for (std::size_t way = 0; way < rims.size(); ++way)
        {
          rims.at(way) = Farther(arranged.Points[vertex], arranged.Points[rims.at(way)], way)
                           ? vertex
                           : rims.at(way);
        }
      }
      std::optional<Tally> base;
      for (std::size_t way = 0; way < rims.size(); ++way)
      {
        const auto side = Outer(arranged, ferrying, {rims.at(way), way});
        const auto outer =
            side ? std::optional<Tally>{Without({}, Stacked(lifts[side->first], side->second))}
                 : std::nullopt;
        if (!outer || (base && *base != *outer))
        {
          return false;
        }
        base = outer;
      }
      for (Tally &lift : lifts)
      {
        lift = Stacked(base.value_or(Tally{}), lift);
      }
      return true;
    }

    /**
     *  The halves of the pieces whose sides the rule keeps one of, each keyed by the vertex it
     *  leaves, then the half, then whether it bounds what is kept, in that order; the halves about
     *  a vertex of more than two then put counterclockwise, as Thronged orders them all.
     *
     *  @param  arranged The arrangement.
     *  @param  ferrying The windings ferried, their snags ordered.
     *  @param  sides    Each unknown's value, and the rule.
     *  @return          The keys.
     */
    auto Bordering(const Arranged &arranged, const Ferrying &ferrying,
                   const std::pair<std::span<const Tally>, Keeps> &sides) -> List<std::uint64_t>
    {
      const auto [lifts, keeps] = sides;
      List<std::uint64_t> keys;
      for (std::size_t edge = 0; edge < arranged.Edges.size(); ++edge)
      {
        const Hauled         &haul   = ferrying.Hauls[edge];
        const Arranged::Edge &piece  = arranged.Edges[edge];
        const Tally           left   = Stacked(lifts[haul.Stratum], haul.Left);
        const Tally           right  = Without(left, piece.Weighs);
        const bool            inside = keeps(left[0] != 0, left[1] != 0);
        if (inside != keeps(right[0] != 0, right[1] != 0))
        {
          const std::uint64_t half = 2 * std::uint64_t{edge};
          keys.push_back((std::uint64_t{piece.Ends[0]} << 33U) | (half << 1U) | (inside ? 1U : 0U));
          keys.push_back(
              (std::uint64_t{piece.Ends[1]} << 33U) | ((half + 1) << 1U) | (inside ? 0U : 1U));
        }
      }
      std::ranges::sort(keys);
      List<std::size_t>   halves;
      List<std::uint64_t> group;
      for (std::size_t start = 0; start < keys.size();)
      {
        const std::uint64_t vertex = keys[start] >> 33U;
        std::size_t         stop   = start;
        halves.clear();
        for (; stop < keys.size() && keys[stop] >> 33U == vertex; ++stop)
        {
          halves.push_back((keys[stop] >> 1U) & 0xFFFF'FFFFU);
        }
        const std::span<const std::size_t> fan = Snarled(ferrying, vertex);
        if (halves.size() > 2 && fan.empty())
        {
          /**
           *  A lone crossing's four halves, their ways farther apart than the angle, as Ferried
           *  told them: in order by their ways alone, as Fanned leaves them.
           */
          std::ranges::sort(
              halves, [&arranged](std::size_t lhs, std::size_t rhs) -> bool
              { return Preceding(Outbound(arranged, lhs), Outbound(arranged, rhs)) < 0; });
        }
        group.clear();
        const auto own = std::span<const std::uint64_t>(keys).subspan(start, stop - start);
        for (const std::size_t half : fan.empty() ? std::span<const std::size_t>(halves) : fan)
        {
          std::ranges::copy_if(own, std::back_inserter(group), [half](std::uint64_t key) -> bool
                               { return ((key >> 1U) & 0xFFFF'FFFFU) == half; });
        }
        std::ranges::copy(group, keys.begin() + static_cast<std::ptrdiff_t>(start));
        start = stop;
      }
      return keys;
    }

    /**
     *  The loops the bounding halves make, each half followed by the first bounding one clockwise
     *  from its twin about the vertex it ends at, as Stitched follows them.
     *
     *  @param  arranged The arrangement.
     *  @param  keys     The keys of the boundary's halves, each vertex's counterclockwise.
     *  @return          The loops.
     */
    auto Threading(const Arranged &arranged, std::span<const std::uint64_t> keys) -> Region
    {
      List<std::uint32_t> places(2 * arranged.Edges.size());
      for (std::size_t place = 0; place < keys.size(); ++place)
      {
        places[(keys[place] >> 1U) & 0xFFFF'FFFFU] = static_cast<std::uint32_t>(place);
      }
      List<std::uint8_t> done(keys.size(), 0);
      List<std::size_t>  run;
      Region             region = {};
      for (std::size_t first = 0; first < keys.size(); ++first)
      {
        if ((keys[first] & 1U) == 0 || done[first] != 0)
        {
          continue;
        }
        run.clear();
        std::size_t place = first;
        for (bool round = false; !round; round = place == first)
        {
          done[place]              = 1;
          const std::uint64_t half = (keys[place] >> 1U) & 0xFFFF'FFFFU;
          run.push_back(half);
          /**
           *  The keys about a vertex alternate between halves that bound and halves that do not, as
           *  only the boundary's are kept: so the one just clockwise of the twin bounds.
           */
          const std::size_t   twin   = places[half ^ 1U];
          const std::uint64_t vertex = keys[twin] >> 33U;
          place                      = twin;
          if (twin > 0 && keys[twin - 1] >> 33U == vertex)
          {
            --place;
            continue;
          }
          while (place + 1 < keys.size() && keys[place + 1] >> 33U == vertex)
          {
            ++place;
          }
        }
        region.Loops.push_back(Sewn(arranged, run));
      }
      return region;
    }
  } // namespace

  /**
   *  Each run's pieces lie in a row among the edges, in its curves' order and along each, as Hewn
   *  makes them: so each run is walked once, and the windings are told of every piece before any is
   *  ordered about a vertex.
   */
  auto Ferried(const Arranged &arranged, Keeps keeps, Region &region) -> Ferry
  {
    if (arranged.Lone.empty())
    {
      return Ferry::Unfit;
    }
    Ferrying    ferrying = {.Hauls   = List<Hauled>(arranged.Edges.size()),
                            .Crossed = List<Rendezvous>(arranged.Points.size())};
    std::size_t cursor   = 0;
    for (std::size_t run = 1; run < arranged.Runs.size(); ++run)
    {
      const std::size_t first = cursor;
      while (cursor < arranged.Edges.size() && arranged.Edges[cursor].Source < arranged.Runs[run])
      {
        ++cursor;
      }
      if (first < cursor && !Haul(arranged, {first, cursor}, ferrying))
      {
        return Ferry::Unfit;
      }
    }
    Snagging(arranged, ferrying);
    if (ferrying.Broken)
    {
      return Ferry::Twisted;
    }
    List<Tally> lifts(ferrying.Unknowns.Size());
    if (!Anchoring(arranged, ferrying, lifts))
    {
      return Ferry::Unfit;
    }
    region = Threading(arranged, Bordering(arranged, ferrying, {lifts, keeps}));
    return Ferry::Done;
  }
} // namespace Geoxy
