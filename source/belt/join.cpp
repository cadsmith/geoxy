/***************************************************************************************************
 *  @file      source/belt/join.cpp
 *  @brief     The stretches kept joined end to start into the offset's polylines.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "belt/join.hpp"

#include "belt/carrier.hpp"
#include "belt/make.hpp"
#include "belt/near.hpp"
#include "belt/piece.hpp"
#include "belt/work.hpp"
#include "scratch.hpp"
#include "twofold.hpp"

#include "geoxy/offset.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <span>
#include <utility>
#include <vector>

namespace Geoxy::Belt
{
  namespace
  {
    /**
     *  The stretches as loops are made of them.
     */
    struct Joining
    {
      const Work *Curves = nullptr;       /**< The curves and the stretches, in order of
                                               their starts.                             */
      const Nodes      *Points = nullptr; /**< The points' numbers. */
      List<std::size_t> Starts = {};      /**< Where the stretches leaving each point
                                               start among them, and one past the last.  */
      List<std::uint8_t> Used = {};       /**< Which stretches are in a chain. */
      List<std::size_t>  Left = {};       /**< For each point the chain being made has left,
                                               one more than where in the chain; 0 for the
                                               others.                                    */
    };

    /**
     *  The stretch a chain goes on by from the end of one: of those leaving that point not yet in a
     *  chain, the first met turning clockwise from the way back.
     *
     *  @param  joining The stretches.
     *  @param  last    The chain's last stretch.
     *  @return         The stretch; none where none leaves.
     */
    auto Next(const Joining &joining, std::size_t last) -> std::size_t
    {
      const Work    &work = *joining.Curves;
      const Stretch &now  = work.Stretches[last];
      /**
       *  Mostly one stretch leaves, or none: no angle is worked out to choose it.
       */
      std::size_t only   = None;
      std::size_t unused = 0;
      for (std::size_t index = joining.Starts[now.End]; index < joining.Starts[now.End + 1];
           ++index)
      {
        if (joining.Used[index] == 0)
        {
          only = index;
          ++unused;
        }
      }
      if (unused < 2)
      {
        return only;
      }
      const Vector into = Heading(work.Curves[now.Curve], joining.Points->At(now.End));
      const double back = std::atan2(-into.Y, -into.X);
      std::size_t  next = None;
      double       best = std::numeric_limits<double>::infinity();
      for (std::size_t index = joining.Starts[now.End]; index < joining.Starts[now.End + 1];
           ++index)
      {
        const Stretch &leaving = work.Stretches[index];
        const Vector   out = Heading(work.Curves[leaving.Curve], joining.Points->At(leaving.Start));
        double         turn = std::fmod(back - std::atan2(out.Y, out.X), Whole);
        if (turn <= 1e-9)
        {
          turn += Whole;
        }
        /**
         *  One already in a chain is never chosen: it is put two turns behind the others.
         */
        turn += 2.0 * Whole * static_cast<double>(joining.Used[index]);
        if (turn < best)
        {
          best = turn;
          next = index;
        }
      }
      return next;
    }

    /**
     *  A line's carrier: the one it was drawn by, across a corner or an end, or its piece's moved.
     *
     *  @param  work  The curves.
     *  @param  curve The curve.
     *  @return       Its carrier.
     */
    auto Borne(const Work &work, std::size_t curve) -> Carrier
    {
      if (curve < work.Made.size() && work.Made[curve] != 0)
      {
        return work.WideCurves[curve];
      }
      const std::size_t source = work.Curves[curve].Source;
      return Shifted(Carried(work.Shape[source], work.Vertices[source].Bulge),
                     work.Shape[source].Left, work.Ways[curve]);
    }

    /**
     *  Whether two stretches in a row are one edge: along one line or one circle, running on.
     *
     *  @param  work   The curves and the reach.
     *  @param  first  One stretch.
     *  @param  second The next.
     *  @return        Whether they are.
     */
    auto Same(const Work &work, const Stretch &first, const Stretch &second) -> bool
    {
      const Piece &one   = work.Curves[first.Curve];
      const Piece &other = work.Curves[second.Curve];
      if (first.Curve == second.Curve || one.Round != other.Round)
      {
        return first.Curve == second.Curve;
      }
      if (one.Round)
      {
        /**
         *  One circle where center and radius are the same doubles and, for two pieces' circles,
         *  the same in two doubles a number: two circles nearer each other than doubles tell are
         *  two, and the edge turns from one to the other at a vertex.
         */
        if (one.Left != other.Left
            || !(one.Center == other.Center)
            || !std::equal_to<double>{}(one.Radius, other.Radius))
        {
          return false;
        }
        if (one.Source == other.Source)
        {
          return true;
        }
        /**
         *  Each as its piece's carrier moved: a circle about a vertex has no piece and is one with
         *  no piece's circle, so the last piece stands in for none, and tells them two.
         */
        const auto borne = [&work](std::size_t curve) -> std::array<double, 6>
        {
          const std::size_t source = std::min(work.Curves[curve].Source, work.Shape.size() - 1);
          const Carrier     wide = Shifted(Carried(work.Shape[source], work.Vertices[source].Bulge),
                                           work.Shape[source].Left, work.Ways[curve]);
          return {wide.X.High, wide.X.Low,       wide.Y.High,
                  wide.Y.Low,  wide.Radius.High, wide.Radius.Low};
        };
        return borne(first.Curve) == borne(second.Curve);
      }
      /**
       *  One line where the two lines, each as it was drawn in two doubles a number, run one way
       *  and the one passes through the other, as near as two doubles tell: two lines that only
       *  doubles take for one are two, and the edge turns from one to the other at a vertex. Two
       *  that turn by more than their rounded ends may make them are two with no more read.
       */
      const Vector way  = one.End - one.Start;
      const Vector next = other.End - other.Start;
      const double round =
          64.0
          * std::numeric_limits<double>::epsilon()
          * work.Scale
          * (std::abs(way.X) + std::abs(way.Y) + std::abs(next.X) + std::abs(next.Y));
      return Dot(way, next) > 0.0
          && std::abs(Cross(way, next)) <= round
          && Alike(work, Borne(work, first.Curve), Borne(work, second.Curve));
    }

    /**
     *  An edge of a loop: stretches in a row along one line or one circle.
     */
    struct Edge
    {
      std::size_t First = 0; /**< Its first stretch. */
      std::size_t Last  = 0; /**< Its last.          */
    };

    /**
     *  The bulge of an edge.
     *
     *  @param  joining The stretches.
     *  @param  edge    The edge.
     *  @return         0 for a line; for an arc that of its sweep from its start to its end.
     */
    auto Bulged(const Joining &joining, const Edge &edge) -> double
    {
      const Work    &work  = *joining.Curves;
      const Stretch &first = work.Stretches[edge.First];
      const Stretch &last  = work.Stretches[edge.Last];
      const Piece   &piece = work.Curves[first.Curve];
      if (!piece.Round)
      {
        return 0.0;
      }
      if (piece.Shallow)
      {
        /**
         *  A nearly straight arc's from its chord and its radius.
         */
        Piece stretch = piece;
        stretch.Start = joining.Points->At(first.Start);
        stretch.End   = joining.Points->At(last.End);
        return Bowed(stretch);
      }
      /**
       *  The tangent of a quarter of the turn from its start to its end about the center, with no
       *  angle worked out: the sine and cosine of half the turn, to one scale, are taken from the
       *  two ways' dot and cross products by whichever sum does not cancel, and halved again.
       */
      const Vector from   = joining.Points->At(first.Start) - piece.Center;
      const Vector onto   = joining.Points->At(last.End) - piece.Center;
      const double along  = Dot(from, onto);
      const double aside  = piece.Left ? Cross(from, onto) : -Cross(from, onto);
      const double size   = std::sqrt(Dot(from, from) * Dot(onto, onto));
      const double flip   = aside < 0.0 ? -1.0 : 1.0;
      const double sine   = along >= 0.0 ? flip * aside : size - along;
      const double cosine = along >= 0.0 ? flip * (size + along) : aside;
      const double whole  = std::sqrt((sine * sine) + (cosine * cosine));
      const double bulge  = cosine >= 0.0 ? sine / (whole + cosine) : (whole - cosine) / sine;
      return piece.Left ? bulge : -bulge;
    }

    /**
     *  A loop of stretches, joined end to start.
     */
    struct Chain
    {
      List<std::size_t> Links = {}; /**< The stretches, in order, the last ending where the
                                         first starts.                                     */
      bool Open = false;            /**< Whether it does not: a chain from an end to an end. */
    };

    /**
     *  The stretch a chain goes on by: the next from its last stretch, or, where none leaves the
     *  point that one ends at, from the one before with that one taken back off the chain, and so
     *  on back. A stretch taken back leads nowhere: a piece left where two curves touch within the
     *  reach, of no loop.
     *
     *  @param  joining The stretches.
     *  @param  links   The chain, which loses the stretches that lead nowhere.
     *  @return         The stretch; none where the whole chain led nowhere.
     */
    auto Onward(Joining &joining, List<std::size_t> &links) -> std::size_t
    {
      std::size_t next = None;
      while (!links.empty() && next == None)
      {
        next = Next(joining, links.back());
        if (next == None)
        {
          joining.Left[joining.Curves->Stretches[links.back()].Start] = 0;
          links.pop_back();
        }
      }
      return next;
    }

    /**
     *  The loops a chain from a stretch not yet in one comes to.
     *
     *  The chain goes on from stretch to stretch. Each time it comes back to a point it has left,
     *  the stretches since are a loop, taken off it, and it goes on from that point.
     *
     *  @param  joining The stretches; those of the chain are marked.
     *  @param  first   The stretch.
     *  @return         The loops, each closed.
     */
    auto Chained(Joining &joining, std::size_t first) -> List<Chain>
    {
      const Work       &work = *joining.Curves;
      List<Chain>       loops;
      List<std::size_t> links;
      links.reserve(work.Stretches.size());
      for (std::size_t link = first; link != None; link = Onward(joining, links))
      {
        const Stretch &now      = work.Stretches[link];
        joining.Used[link]      = 1;
        joining.Left[now.Start] = links.size() + 1;
        links.push_back(link);
        const std::size_t since = joining.Left[now.End];
        if (since != 0)
        {
          const auto from = links.begin() + static_cast<std::ptrdiff_t>(since - 1);
          loops.push_back({.Links = {from, links.end()}});
          for (const std::size_t taken : loops.back().Links)
          {
            joining.Left[work.Stretches[taken].Start] = 0;
          }
          links.erase(from, links.end());
        }
      }
      return loops;
    }

    /**
     *  A chain's edges: stretches in a row along one line or one circle made one, the last with the
     *  first where a closed chain starts in the middle of an edge.
     *
     *  @param  work  The curves and the stretches.
     *  @param  chain The chain.
     *  @return       The edges.
     */
    auto Edges(const Work &work, const Chain &chain) -> List<Edge>
    {
      List<Edge> edges;
      edges.reserve(chain.Links.size());
      for (const std::size_t link : chain.Links)
      {
        if (!edges.empty() && Same(work, work.Stretches[edges.back().Last], work.Stretches[link]))
        {
          edges.back().Last = link;
        }
        else
        {
          edges.push_back({.First = link, .Last = link});
        }
      }
      if (!chain.Open
          && edges.size() > 1
          && Same(work, work.Stretches[edges.back().Last], work.Stretches[edges.front().First]))
      {
        edges.front().First = edges.back().First;
        edges.pop_back();
      }
      return edges;
    }

    /**
     *  A polyline run the other way: each vertex takes the bulge of the edge that arrived at it,
     *  its sign changed.
     *
     *  @param  loop The polyline.
     *  @return      The same one from its last vertex to its first.
     */
    auto Flipped(const Polyline &loop) -> Polyline
    {
      const std::size_t count  = loop.Vertices.size();
      Polyline          turned = {.Vertices = {}, .Closed = loop.Closed};
      turned.Vertices.reserve(count);
      for (std::size_t index = 0; index < count; ++index)
      {
        const std::size_t from = count - 1 - index;
        const std::size_t edge = (from + count - 1) % count;
        const bool        last = !loop.Closed && from == 0;
        turned.Vertices.push_back({.Position = loop.Vertices[from].Position,
                                   .Bulge    = last ? 0.0 : -loop.Vertices[edge].Bulge});
      }
      return turned;
    }

    /**
     *  A chain as a polyline: a vertex where each edge starts, with the edge's bulge; a whole
     *  circle as two half turns from its rightmost point.
     *
     *  @param  joining The stretches.
     *  @param  chain   The chain.
     *  @return         The polyline, the band on its left.
     */
    auto Drawn(const Joining &joining, const Chain &chain) -> Polyline
    {
      const Work      &work  = *joining.Curves;
      const List<Edge> edges = Edges(work, chain);
      Polyline         loop  = {.Vertices = {}, .Closed = !chain.Open};
      const Piece     &first = work.Curves[work.Stretches[edges.front().First].Curve];
      if (edges.size() == 1 && !chain.Open)
      {
        const double bulge = first.Left ? 1.0 : -1.0;
        const Vector reach = {.X = first.Radius};
        loop.Vertices      = {{.Position = first.Center + reach, .Bulge = bulge},
                              {.Position = first.Center - reach, .Bulge = bulge}};
        return loop;
      }
      loop.Vertices.reserve(edges.size() + 1);
      for (const Edge &edge : edges)
      {
        loop.Vertices.push_back({.Position = joining.Points->At(work.Stretches[edge.First].Start),
                                 .Bulge    = Bulged(joining, edge)});
      }
      if (chain.Open)
      {
        loop.Vertices.push_back(
            {.Position = joining.Points->At(work.Stretches[edges.back().Last].End), .Bulge = 0.0});
      }
      return loop;
    }

    /**
     *  Whether a loop of the band's edge is an edge of the region asked: both ways, or of an open
     *  shape, every one; else one that lies on the side of what the shape fills that is asked,
     *  outside it outward and inside it inward. The band's edge never crosses the shape, so what
     *  the shape fills is the same all along a loop, and is read at the middle of its first
     *  stretch.
     *
     *  @param  joining The shape and the stretches.
     *  @param  chain   The loop.
     *  @return         Whether it is.
     */
    auto Sided(const Joining &joining, const Chain &chain) -> bool
    {
      const Work &work = *joining.Curves;
      if (!work.Closed || work.Asked.Sides == Side::Both)
      {
        return true;
      }
      const Stretch &first  = work.Stretches[chain.Links.front()];
      const Point    middle = Mid(work.Curves[first.Curve],
                                  {joining.Points->At(first.Start), joining.Points->At(first.End)},
                                  first.Onto - first.From);
      const int      windings = Wound(work.Shape, middle);
      const bool     filled = work.Asked.Fill == Fills::EvenOdd ? windings % 2 != 0 : windings != 0;
      return filled == (work.Asked.Sides == Side::Inward);
    }

    /**
     *  A chain from a stretch on: each stretch followed by the one that leaves where it ends, the
     *  first met turning clockwise where several do, till none leaves or the chain is back where it
     *  started.
     *
     *  @param  joining The stretches, those in the chain marked.
     *  @param  first   The stretch it starts with.
     *  @return         The chain, open unless it came back to its start.
     */
    auto Strung(Joining &joining, std::size_t first) -> Chain
    {
      const Work       &work  = *joining.Curves;
      const std::size_t start = work.Stretches[first].Start;
      Chain             chain = {.Links = {}, .Open = true};
      for (std::size_t link = first; link != None; link = Next(joining, chain.Links.back()))
      {
        joining.Used[link] = 1;
        chain.Links.push_back(link);
        if (work.Stretches[link].End == start)
        {
          chain.Open = false;
          break;
        }
      }
      return chain;
    }
  } // namespace

  auto Looped(const Work &work, const Nodes &points) -> std::vector<Polyline>
  {
    Joining joining = {.Curves = &work, .Points = &points};
    joining.Starts.assign(points.Count() + 1, 0);
    for (const Stretch &stretch : work.Stretches)
    {
      ++joining.Starts[stretch.Start + 1];
    }
    for (std::size_t node = 0; node < points.Count(); ++node)
    {
      joining.Starts[node + 1] += joining.Starts[node];
    }
    joining.Used.assign(work.Stretches.size(), 0);
    joining.Left.assign(points.Count(), 0);
    const bool            inward = work.Closed && work.Asked.Sides == Side::Inward;
    std::vector<Polyline> loops;
    for (std::size_t first = 0; first < work.Stretches.size(); ++first)
    {
      if (joining.Used[first] != 0)
      {
        continue;
      }
      for (const Chain &chain : Chained(joining, first))
      {
        if (!Sided(joining, chain))
        {
          continue;
        }
        Polyline loop = Drawn(joining, chain);
        loops.push_back(inward ? Flipped(loop) : std::move(loop));
      }
    }
    return loops;
  }

  void Flanking(Work &work)
  {
    std::ranges::sort(
        work.Stretches, [](const Stretch &first, const Stretch &second) -> bool
        { return std::pair{first.Start, first.End} < std::pair{second.Start, second.End}; });
    List<Stretch> sided;
    for (std::size_t start = 0; start < work.Stretches.size();)
    {
      std::size_t     end  = start;
      bool            side = true;
      const std::pair ends = {work.Stretches[start].Start, work.Stretches[start].End};
      for (; end < work.Stretches.size()
             && std::pair{work.Stretches[end].Start, work.Stretches[end].End} == ends;
           ++end)
      {
        side = side && work.Curves[work.Stretches[end].Curve].Flank == work.Flank;
      }
      if (side)
      {
        sided.push_back(work.Stretches[start]);
      }
      start = end;
    }
    work.Stretches = std::move(sided);
  }

  auto Sidelong(const Work &work, const Nodes &points) -> std::vector<Polyline>
  {
    Joining joining = {.Curves = &work, .Points = &points};
    joining.Starts.assign(points.Count() + 1, 0);
    List<std::uint32_t> arrive(points.Count(), 0);
    for (const Stretch &stretch : work.Stretches)
    {
      ++joining.Starts[stretch.Start + 1];
      ++arrive[stretch.End];
    }
    for (std::size_t node = 0; node < points.Count(); ++node)
    {
      joining.Starts[node + 1] += joining.Starts[node];
    }
    joining.Used.assign(work.Stretches.size(), 0);
    std::vector<Polyline> chains;
    /**
     *  First the chains that start where nothing arrives, then whatever is left.
     */
    for (const bool ends : {true, false})
    {
      for (std::size_t first = 0; first < work.Stretches.size(); ++first)
      {
        if (joining.Used[first] != 0 || (ends && arrive[work.Stretches[first].Start] != 0))
        {
          continue;
        }
        const Polyline chain = Drawn(joining, Strung(joining, first));
        /**
         *  The band lies on every curve's left: on the shape's left the curves run against it.
         */
        chains.push_back(work.Flank < 0 ? Flipped(chain) : chain);
      }
    }
    return chains;
  }
} // namespace Geoxy::Belt
