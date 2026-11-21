/***************************************************************************************************
 *  @file      source/belt/keep.cpp
 *  @brief     The stretches of a band's curves that lie on its edge, kept.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "belt/keep.hpp"

#include "belt/carrier.hpp"
#include "belt/make.hpp"
#include "belt/near.hpp"
#include "belt/piece.hpp"
#include "belt/work.hpp"
#include "scratch.hpp"
#include "twofold.hpp"

#include "geoxy/box.hpp"
#include "geoxy/offset.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <span>
#include <tuple>
#include <utility>
#include <vector>

namespace Geoxy::Belt
{
  namespace
  {
    /**
     *  How a stretch's middle stands to the other parts of the shape.
     */
    struct Measure
    {
      int  Stands = 1;     /**< Below 0 nearer than the distance, above 0 farther, 0 just it. */
      bool Behind = false; /**< Whether a part just the distance away lies on the curve's right:
                                the band is then on both sides of the stretch.                    */
    };

    /**
     *  How a stretch whose middle lies the distance from another part of the shape, as near as
     *  doubles tell, stands to that part in two doubles a number.
     *
     *  @param  work   The shape and the curves.
     *  @param  curve  The stretch's curve.
     *  @param  middle The stretch's middle.
     *  @return        How it stands, and whether a part just the distance away lies past it.
     */
    auto Measured(Work &work, std::size_t curve, Point middle) -> Measure
    {
      const Piece  &piece    = work.Curves[curve];
      const double  distance = work.Asked.Distance;
      const Carrier wide     = Widened(work, curve);
      const double  slack    = 1e-24 * distance * distance;
      const Vector  ahead    = Heading(piece, middle);
      const Vector  right    = {.X = ahead.Y, .Y = -ahead.X};
      double        least    = std::numeric_limits<double>::infinity();
      bool          behind   = false;
      for (const std::size_t other : Around(work, middle))
      {
        if (Beyond(work.Pads[other], middle) || Spared(work, {&piece, middle}, other))
        {
          continue;
        }
        const Piece &shape = work.Shape[other];
        const bool   past  = Past(work, other, middle);
        double       clear = 0.0;
        Vector       away  = {};
        if (past)
        {
          /**
           *  Past a flat end the piece holds only up to its other end's vertex: how far that is, in
           *  two doubles too.
           */
          const Point  hub = other == 0 ? shape.End : shape.Start;
          const double far = Hub(work, {&piece, middle}, other);
          clear = far < std::numeric_limits<double>::infinity()
                    ? Clear({.On = &wide, .Near = middle, .Distance = distance}, hub).High
                    : far;
          away  = middle - hub;
        }
        else
        {
          clear =
              Clear({.On = &wide, .Near = middle, .Distance = distance}, shape, Broad(work, other))
                  .High;
          away = Away(shape, middle);
        }
        least = std::min(least, clear);
        /**
         *  Just the distance from this part, the band it makes lies past the curve where the part
         *  itself does: on the curve's right.
         */
        if (std::abs(clear) <= slack && Dot(away, right) < 0.0)
        {
          behind = true;
        }
      }
      return {.Stands = least < -slack ? -1 : static_cast<int>(least > slack), .Behind = behind};
    }

    /**
     *  How a curve's point that lies on a mitre's edge, as near as doubles tell, stands to that
     *  mitre in two doubles a number: read against each curve of the mitre's edge it lies on.
     *
     *  @param  work  The curves and the wedges.
     *  @param  pair  The mitre's number, and the number of the curve the point is of.
     *  @param  point The point.
     *  @return       Below 0 inside the mitre, 0 on its edge as near as two doubles tell or by its
     *                way from the corner, above 0 outside it.
     */
    auto Rimmed(Work &work, const std::pair<std::size_t, std::size_t> &pair, Point point) -> int
    {
      const double  tie    = work.Scale * 1e-9;
      const Carrier wide   = Widened(work, pair.second);
      bool          inside = true;
      bool          read   = false;
      for (std::size_t index = 0; index < work.Wedges[pair.first].Rim.size(); ++index)
      {
        const std::size_t edge = work.Wedges[pair.first].Edges[index];
        if (edge == None || Apart(work.Wedges[pair.first].Rim[index], point) > tie)
        {
          continue;
        }
        const Carrier         rim   = Widened(work, edge);
        const Piece          &curve = work.Curves[edge];
        const Twofold::Double value = Beside({.On = &wide, .Near = point}, rim);
        /**
         *  The band lies on the curve's left: within its circle where it runs counterclockwise.
         */
        const Vector way  = {.X = rim.WayX.High, .Y = rim.WayY.High};
        const bool   flip = curve.Round ? curve.Left : Dot(curve.End - curve.Start, way) < 0.0;
        const double left = flip ? -value.High : value.High;
        const double size = (curve.Round ? curve.Radius : Long(way)) * work.Scale;
        if (left < -1e-24 * size)
        {
          return 1;
        }
        if (!(left > 1e-24 * size))
        {
          inside = false;
        }
        read = true;
      }
      return read && inside ? -1 : 0;
    }

    /**
     *  How a point of a circle about a corner stands to that corner's own mitre: inside it just
     *  where the circle runs between the ways to the mitre's two ends, however thin the mitre. That
     *  is read from the ways, not measured.
     *
     *  @param  work  The reach and the distance.
     *  @param  wedge The wedge.
     *  @param  from  The curve the point is of, and the point.
     *  @return       Below 0 inside the mitre, above 0 outside it; 0 where the curve is no circle
     *                about the wedge's corner, or the wedge no mitre.
     */
    auto Spanned(const Work &work, const Wedge &wedge, const std::pair<const Piece *, Point> &from)
        -> int
    {
      const Piece *curve = from.first;
      if (!curve->Round || curve->Source != None || !wedge.Coned)
      {
        return 0;
      }
      const Vector apart = curve->Center - wedge.About;
      if (Dot(apart, apart) > 0.0)
      {
        return 0;
      }
      const Vector off   = from.second - wedge.About;
      const double slack = work.Reach * work.Asked.Distance;
      return Cross(wedge.Span[0], off) > slack && Cross(off, wedge.Span[1]) > slack ? -1 : 1;
    }

    /**
     *  How a point stands to one mitre or square: by its box, then by how near its edge it lies; on
     *  its edge as near as doubles tell, a curve's point is read against that edge again in two
     *  doubles a number.
     *
     *  @param  work  The curves and the wedges.
     *  @param  pair  The wedge's number, and the number of the curve the point is of.
     *  @param  point The point.
     *  @return       Below 0 inside it, 0 on its edge, above 0 outside it.
     */
    auto Stood(Work &work, const std::pair<std::size_t, std::size_t> &pair, Point point) -> int
    {
      const Wedge &wedge = work.Wedges[pair.first];
      const double tie   = work.Scale * 1e-9;
      const Vector pad   = {.X = tie, .Y = tie};
      if (Beyond({.Min = wedge.Hull.Min - pad, .Max = wedge.Hull.Max + pad}, point))
      {
        return 1;
      }
      double edge = std::numeric_limits<double>::infinity();
      for (const Piece &piece : wedge.Rim)
      {
        edge = std::min(edge, Apart(piece, point));
      }
      if (edge > tie)
      {
        return Wound(wedge.Rim, point) != 0 ? -1 : 1;
      }
      return Rimmed(work, pair, point);
    }

    /**
     *  How a curve's point stands to the mitres and squares, the curve's own left out; and, where
     *  it lies on the edge of one, whether that one lies on the curve's right there: read from
     *  which way the nearest piece of its edge runs, the mitre being on one side of it all along.
     *
     *  @param  work  The wedges.
     *  @param  point The point.
     *  @param  from  The curve the point is of.
     *  @return       Below 0 inside one of them, 0 on the edge of one and inside none, above 0
     *                outside them all; and whether one whose edge it is on lies past it.
     */
    auto Wedging(Work &work, Point point, std::size_t from) -> Measure
    {
      const Piece &curve  = work.Curves[from];
      const Vector ahead  = Heading(curve, point);
      Measure      stands = {};
      for (std::size_t index = 0; index < work.Wedges.size(); ++index)
      {
        const Wedge &wedge   = work.Wedges[index];
        const int    spanned = Spanned(work, wedge, {&curve, point});
        if (spanned < 0)
        {
          return {.Stands = -1};
        }
        if (spanned > 0 || index == curve.Bounds)
        {
          continue;
        }
        const int stood = Stood(work, {index, from}, point);
        if (stood < 0)
        {
          return {.Stands = -1};
        }
        if (stood == 0)
        {
          const Piece *nearest = wedge.Rim.data();
          for (const Piece &piece : wedge.Rim)
          {
            nearest = Apart(piece, point) < Apart(*nearest, point) ? &piece : nearest;
          }
          const bool past = (Dot(Heading(*nearest, point), ahead) < 0.0) == wedge.Left;
          stands.Stands   = 0;
          if (past)
          {
            stands.Behind = true;
          }
        }
      }
      return stands;
    }

    /**
     *  Whether a stretch of a curve lies on the region's edge: its middle lies outside the band, no
     *  nearer the shape than the distance and in no mitre; where it lies just on the band's edge
     *  elsewhere too, the band is not on its far side as well; and for a closed shape asked inward
     *  or outward its middle lies on the side of what the shape fills that is asked.
     *
     *  @param  work   The shape and the curves.
     *  @param  curve  The curve.
     *  @param  middle The stretch's middle.
     *  @return        Whether it does.
     */
    auto Edging(Work &work, std::size_t curve, Point middle) -> bool
    {
      const Piece &piece    = work.Curves[curve];
      const double distance = work.Asked.Distance;
      const double tie      = work.Scale * 1e-9;
      const double nearest  = Nearest(work, piece, middle);
      /**
       *  A line across a corner or an end is drawn from rounded points: where it lies the distance
       *  from the shape as near as doubles tell, two doubles tell no more, and it is taken to lie
       *  just the distance from it.
       */
      Measure measured = {.Stands = nearest < distance - tie ? -1 : 1};
      if (nearest >= distance - tie && nearest <= distance + tie)
      {
        measured = Measured(work, curve, middle);
      }
      const Measure wedged = Wedging(work, middle, curve);
      /**
       *  Just on the band's edge by another part of it, the stretch is an edge of the region only
       *  if the band is not on its far side too: that part, or that mitre, lies on its right.
       */
      return measured.Stands >= 0 && wedged.Stands >= 0 && !measured.Behind && !wedged.Behind;
    }

    /**
     *  The cuts put in order: curve after curve, and along each curve. Each is dealt to its curve's
     *  place among them, counted first, and slid back past those of its curve that come after it: a
     *  curve has few.
     *
     *  @param work The curves and their cuts.
     */
    void Ordered(Work &work)
    {
      List<std::size_t> next(work.Curves.size() + 1, 0);
      for (const Cutting &cut : work.Cuts)
      {
        ++next[cut.Curve + 1];
      }
      for (std::size_t curve = 0; curve < work.Curves.size(); ++curve)
      {
        next[curve + 1] += next[curve];
      }
      const List<std::size_t> starts = next;
      /**
       *  A copy to deal them into: every place of it is written again.
       */
      List<Cutting> dealt(work.Cuts);
      for (const Cutting &cut : work.Cuts)
      {
        std::size_t place = next[cut.Curve]++;
        for (; place > starts[cut.Curve]
               && std::pair{cut.Along, cut.Node}
                      < std::pair{dealt[place - 1].Along, dealt[place - 1].Node};
             --place)
        {
          dealt[place] = dealt[place - 1];
        }
        dealt[place] = cut;
      }
      work.Cuts = std::move(dealt);
    }

    /**
     *  Whether the cut of a number and the next are the two ends of a stretch: on one curve, at two
     *  points.
     *
     *  @param  work  The cuts, in order.
     *  @param  index The first cut's number, not the last cut's.
     *  @return       Whether they are.
     */
    auto Spanning(const Work &work, std::size_t index) -> bool
    {
      const Cutting &from = work.Cuts[index];
      const Cutting &onto = work.Cuts[index + 1];
      return from.Curve == onto.Curve && from.Node != onto.Node;
    }

    /**
     *  What is known of the stretches before any is measured.
     */
    struct Runs
    {
      List<std::size_t> Before = {}; /**< For the cut that starts a stretch: the cut that starts
                                          the stretch it runs on from, or none.                */
      List<std::uint8_t> Known = {}; /**< For the cut that starts a stretch: 1 on the region's
                                          edge, 2 not, 0 not yet told.                         */
    };

    /**
     *  Which way one curve turns from another where they cross, read in two doubles: the cross of
     *  their headings at the point, a line's its way and a circle's square to the way from its
     *  center, each from its carrier.
     *
     *  @param  work   The curves.
     *  @param  curves The two curves.
     *  @param  point  Where they cross, as doubles hold it.
     *  @return        Above 0 where the second heads to the first's left, below 0 to its right; 0
     *                 where the cross is no more than the point's rounding may make it.
     */
    auto Veering(Work &work, std::array<std::size_t, 2> curves, Point point) -> double
    {
      const auto heading = [&work, point](std::size_t curve) -> std::array<Twofold::Double, 2>
      {
        const Carrier wide  = Widened(work, curve);
        const Piece  &piece = work.Curves[curve];
        if (!wide.Round)
        {
          const Vector way  = {.X = wide.WayX.High, .Y = wide.WayY.High};
          const double sign = Dot(piece.End - piece.Start, way) < 0.0 ? -1.0 : 1.0;
          return {Twofold::Multiply(wide.WayX, Twofold::Lifted(sign)),
                  Twofold::Multiply(wide.WayY, Twofold::Lifted(sign))};
        }
        const Twofold::Double outX = Twofold::Difference(Twofold::Lifted(point.X), wide.X);
        const Twofold::Double outY = Twofold::Difference(Twofold::Lifted(point.Y), wide.Y);
        return piece.Left ? std::array{Twofold::Negated(outY), outX}
                          : std::array{outY, Twofold::Negated(outX)};
      };
      const std::array<Twofold::Double, 2> here  = heading(curves[0]);
      const std::array<Twofold::Double, 2> there = heading(curves[1]);
      const double turn  = Twofold::Difference(Twofold::Multiply(here[0], there[1]),
                                               Twofold::Multiply(here[1], there[0]))
                               .High;
      const double slack = 4.0
                         * std::numeric_limits<double>::epsilon()
                         * work.Scale
                         * (std::abs(here[0].High)
                            + std::abs(here[1].High)
                            + std::abs(there[0].High)
                            + std::abs(there[1].High));
      return std::abs(turn) > slack ? turn : 0.0;
    }

    /**
     *  Two curves that cross at a point nothing else is cut at, each within itself: the band lies
     *  on the left of each, so the stretch of each that runs on the other's left is within it and
     *  no part of the edge; and the one that comes in on the other's right is on the edge just
     *  where the one that goes out on the other's right is, the edge turning from the one to the
     *  other there. Two that cross within a thousandth of one heading are left to be measured.
     *
     *  @param work  The curves and their cuts, in order.
     *  @param runs  What is known, added to.
     *  @param cross The two cuts' numbers, and the point.
     */
    void Handed(Work &work, Runs &runs, const std::pair<std::array<std::size_t, 2>, Point> &cross)
    {
      const std::size_t one   = cross.first[0];
      const std::size_t other = cross.first[1];
      const Vector      here  = Heading(work.Curves[work.Cuts[one].Curve], cross.second);
      const Vector      there = Heading(work.Curves[work.Cuts[other].Curve], cross.second);
      double            turn  = Cross(here, there);
      if (!(std::abs(turn) > 1e-3 * Long(here) * Long(there)))
      {
        /**
         *  A shallow crossing's turn is read in two doubles from the curves' carriers, and taken
         *  where it is more than the point's own rounding may make it: two curves that all but
         *  touch cross twice a hair apart, and the stretch of each between the two is told by the
         *  crossings themselves, where no measure of its middle tells it from the band's edge.
         */
        turn = Veering(work, {work.Cuts[one].Curve, work.Cuts[other].Curve}, cross.second);
        if (!(turn > 0.0) && !(turn < 0.0))
        {
          return;
        }
      }
      /**
       *  The curve that heads to the other's left: it goes out within the band and comes in on the
       *  edge's side; the other comes in within the band and goes out on the edge's side.
       */
      const std::size_t lefts  = turn > 0.0 ? other : one;
      const std::size_t rights = turn > 0.0 ? one : other;
      runs.Known[lefts]        = 2;
      runs.Known[rights - 1]   = 2;
      runs.Before[rights]      = lefts - 1;
    }

    /**
     *  What is known of the stretches before any is measured. A stretch runs on from the one that
     *  ends where it starts, at a point one curve ends at, another starts at, and nothing else is
     *  cut at: along such a run no other curve is crossed, so the whole of it lies on the region's
     *  edge or none of it does. And where two curves cross, two of the four stretches are told at
     *  once and one runs on from another.
     *
     *  @param  work   The cuts, in order.
     *  @param  points The points' numbers.
     *  @return        What is known.
     */
    auto Led(Work &work, const Nodes &points) -> Runs
    {
      /**
       *  At each point: how many cuts in all, how many counted, how many of those within their
       *  curves, the first two, and the cuts that end a curve and start one there.
       */
      struct Meeting
      {
        std::uint32_t              All     = 0;
        std::uint32_t              Cuts    = 0;
        std::uint32_t              Inner   = 0;
        std::array<std::size_t, 2> Held    = {None, None};
        std::size_t                Ended   = None;
        std::size_t                Started = None;
      };

      List<Meeting> met(points.Count());
      /**
       *  A curve cut within itself at the point its end is: the end is that cut over again, and is
       *  not counted; the cut beside an end's is its own curve's, each curve having both its ends.
       *  Two pieces in a row that turn by a hair cross a hair before their ends, and the point is a
       *  crossing, which hands the edge on, not a joint.
       */
      const auto twice = [&work](std::size_t end, std::size_t inner) -> bool
      {
        const Cutting &cut   = work.Cuts[end];
        const Cutting &other = work.Cuts[inner];
        return other.Node == cut.Node && other.Along > 0.0 && other.Along < 1.0;
      };
      for (std::size_t index = 0; index < work.Cuts.size(); ++index)
      {
        const Cutting &cut  = work.Cuts[index];
        Meeting       &here = met[cut.Node];
        ++here.All;
        if ((!(cut.Along < 1.0) && twice(index, index - 1))
            || (!(cut.Along > 0.0) && twice(index, index + 1)))
        {
          continue;
        }
        if (here.Cuts < 2)
        {
          here.Held.at(here.Cuts) = index;
        }
        ++here.Cuts;
        if (!(cut.Along < 1.0))
        {
          here.Ended = index - 1;
        }
        else if (!(cut.Along > 0.0))
        {
          here.Started = index;
        }
        else
        {
          ++here.Inner;
        }
      }
      Runs runs = {.Before = List<std::size_t>(work.Cuts.size(), None),
                   .Known  = List<std::uint8_t>(work.Cuts.size(), 0)};
      for (std::size_t node = 0; node < met.size(); ++node)
      {
        const Meeting &here = met[node];
        if (here.Cuts != 2)
        {
          continue;
        }
        /**
         *  Two curves cross, one of them perhaps at its very end, and nothing else is there.
         */
        if (here.Inner == 2 && here.All <= 3)
        {
          Handed(work, runs, {here.Held, points.At(node)});
        }
        else if (here.Started != None)
        {
          runs.Before[here.Started] = here.Ended;
        }
      }
      return runs;
    }
  } // namespace

  void Kept(Work &work, const Nodes &points)
  {
    const Vector pad = {.X = work.Pad, .Y = work.Pad};
    work.Pads.reserve(work.Shape.size());
    for (const Box &box : work.Hulls)
    {
      work.Pads.push_back({.Min = box.Min - pad, .Max = box.Max + pad});
    }
    if (work.Shape.size() <= Handful)
    {
      work.Everyone.resize(work.Shape.size());
      std::ranges::iota(work.Everyone, std::size_t{0});
    }
    else
    {
      work.Near.Lay(work.All, work.Pads);
    }
    Ordered(work);
    Runs                     runs   = Led(work, points);
    const List<std::size_t> &before = runs.Before;
    List<std::uint8_t>      &known  = runs.Known;
    for (std::size_t index = 0; index + 1 < work.Cuts.size(); ++index)
    {
      if (!Spanning(work, index))
      {
        continue;
      }
      /**
       *  Back along the stretches this one runs on from, to one already told or one that runs on
       *  from none, or round to this one again: that one is measured, or its answer taken, for them
       *  all. Most run on from none, and are their own.
       */
      const auto back = [&](std::size_t stretch) -> bool
      {
        return known[stretch] == 0 && before[stretch] != None && before[stretch] != index;
      };
      std::size_t base = index;
      while (back(base))
      {
        base = before[base];
      }
      if (known[base] == 0)
      {
        /**
         *  The longest of them is the one measured: a stretch a hair long, as between the two
         *  points two curves that all but touch cross at, has a middle no measure tells from the
         *  band's edge, and its answer, taken for the long stretches it runs on into, drops them.
         */
        const auto span = [&work, &points](std::size_t stretch) -> double
        {
          const Vector across =
              points.At(work.Cuts[stretch + 1].Node) - points.At(work.Cuts[stretch].Node);
          return Dot(across, across);
        };
        std::size_t longest = index;
        for (std::size_t stretch = index; stretch != base;)
        {
          stretch = before[stretch];
          longest = span(stretch) > span(longest) ? stretch : longest;
        }
        const Cutting &from = work.Cuts[longest];
        const Cutting &onto = work.Cuts[longest + 1];
        const Point    middle =
            Mid(work.Curves[from.Curve], {points.At(from.Node), points.At(onto.Node)},
                onto.Along - from.Along);
        known[base] = Edging(work, from.Curve, middle) ? 1 : 2;
      }
      for (std::size_t stretch = index; stretch != base; stretch = before[stretch])
      {
        known[stretch] = known[base];
      }
      if (known[base] == 1)
      {
        work.Stretches.push_back({.Curve = work.Cuts[index].Curve,
                                  .From  = work.Cuts[index].Along,
                                  .Onto  = work.Cuts[index + 1].Along,
                                  .Start = work.Cuts[index].Node,
                                  .End   = work.Cuts[index + 1].Node});
      }
    }
  }

  void Once(Work &work)
  {
    std::ranges::sort(work.Stretches,
                      [&work](const Stretch &first, const Stretch &second) -> bool
                      {
                        return std::tuple{first.Start, first.End, work.Curves[first.Curve].Round}
                             < std::tuple{second.Start, second.End,
                                          work.Curves[second.Curve].Round};
                      });
    const auto twice =
        std::ranges::unique(work.Stretches, [](const Stretch &first, const Stretch &second) -> bool
                            { return first.Start == second.Start && first.End == second.End; });
    work.Stretches.erase(twice.begin(), twice.end());
  }
} // namespace Geoxy::Belt
