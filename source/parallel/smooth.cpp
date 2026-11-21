/***************************************************************************************************
 *  @file      source/parallel/smooth.cpp
 *  @brief     A curve that is neither line nor arc laid with pairs of arcs tangent to its offset.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "parallel/smooth.hpp"

#include "kernel/piece.hpp"
#include "parallel/bend.hpp"
#include "real.hpp"
#include "scratch.hpp"
#include "trig.hpp"

#include "geoxy/box.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>
#include <variant>

namespace Geoxy::Parallel
{
  namespace
  {
    /**
     *  A point of a run, its offset, and the way the run goes there.
     */
    struct Probe
    {
      Point  At     = {}; /**< The run's point.              */
      Point  Offset = {}; /**< Its offset.                   */
      Vector Way    = {}; /**< The run's unit tangent there. */
    };

    /**
     *  A run's point at a parameter, and its offset.
     *
     *  Where the run stalls, its derivative 0 as at a spline's pole repeated, its way is read a
     *  little inward, where it moves again the way it leaves or comes in.
     *
     *  @param  run      The run.
     *  @param  along    Where, from 0 to 1.
     *  @param  distance The distance.
     *  @return          The point, its offset and the way.
     */
    auto Probed(const Run &run, double along, double distance) -> Probe
    {
      const Kernel::Sample sample = Sampled(run, along);
      Vector               way    = sample.Way;
      if (Zero(Largest(way)))
      {
        way = Sampled(run, along < 0.5 ? along + 0x1p-30 : along - 0x1p-30).Way;
      }
      const Vector unit = way / Length(way);
      return {.At     = sample.At,
              .Offset = sample.At + (Vector{.X = unit.Y, .Y = -unit.X} * distance),
              .Way    = unit};
    }

    /**
     *  A stretch of a run, between two probes.
     */
    struct Tract
    {
      double From  = 0.0; /**< Where it starts.        */
      double To    = 1.0; /**< Where it ends.          */
      Probe  Start = {};  /**< The probe at its start. */
      Probe  End   = {};  /**< The probe at its end.   */
    };

    /**
     *  A segment from one point to another as a moved piece.
     *
     *  @param  start Where it starts.
     *  @param  end   Where it ends.
     *  @param  bulge Its bulge.
     *  @return       The piece.
     */
    auto Piece(Point start, Point end, double bulge) -> Moved
    {
      const std::array<Vertex, 2> pair = {Vertex{.Position = start, .Bulge = bulge},
                                          Vertex{.Position = end}};
      const Bend                  bend = Bent(pair, 0);
      return {.Start  = start,
              .End    = end,
              .Center = bend.Center,
              .Radius = bend.Radius,
              .Sweep  = bend.Sweep};
    }

    /**
     *  The signed angle from one way to another.
     *
     *  @param  from The way turned from.
     *  @param  onto The way turned to.
     *  @return      The angle, counterclockwise positive.
     */
    auto Turn(Vector from, Vector onto) -> double
    {
      return Trig::Atan2(Cross(from, onto), Dot(from, onto));
    }

    /**
     *  Two arcs tangent to each other, from one point along one way to another along another.
     */
    struct Biarc
    {
      Point  Joint = {};  /**< Where the first ends and the second starts. */
      double First = 0.0; /**< The first's bulge.                          */
      double Last  = 0.0; /**< The second's bulge.                         */
    };

    /**
     *  The pair of arcs from one offset point to another, tangent to the offset at both, their
     *  tangents' lengths to the joint equal.
     *
     *  @param  stretch The stretch, its probes' offsets the ends.
     *  @param  side    1 where the offset runs the run's way, -1 where it runs back.
     *  @return         The pair; a straight one where both ways and the chord are one.
     */
    auto Paired(const Tract &stretch, double side) -> Biarc
    {
      const Point  start  = stretch.Start.Offset;
      const Point  end    = stretch.End.Offset;
      const Vector lead   = stretch.Start.Way * side;
      const Vector trail  = stretch.End.Way * side;
      const Vector chord  = end - start;
      const double square = Dot(chord, chord);
      const double turned = 2.0 * (1.0 - Dot(lead, trail));
      double       reach  = 0.0;
      if (turned > 0x1p-40)
      {
        const double lean = Dot(chord, lead + trail);
        reach             = (std::sqrt((lean * lean) + (turned * square)) - lean) / turned;
      }
      else
      {
        reach = square / (4.0 * Dot(chord, lead));
      }
      const Point joint = Midpoint(start + (lead * reach), end - (trail * reach));
      return {.Joint = joint,
              .First = Trig::Tan(Turn(lead, joint - start) / 2.0),
              .Last  = Trig::Tan(Turn(end - joint, trail) / 2.0)};
    }

    /**
     *  A pair of arcs laid along a stretch, and the offset's points inside it.
     */
    struct Fitting
    {
      const Tract         *Stretch = nullptr; /**< The stretch. */
      const Biarc         *Pair    = nullptr; /**< The pair.    */
      std::array<Point, 3> Inner   = {};      /**< The offset a quarter, a half and three
                                                   quarters of the way.                */
    };

    /**
     *  Whether a pair of arcs fits a stretch of the offset: neither turns more than a quarter turn,
     *  and the offset's points inside the stretch lie within an aim of one of them.
     *
     *  @param  fitting The stretch, the pair and the points.
     *  @param  aim     How near the points must lie.
     *  @return         Whether it fits.
     */
    auto Fits(const Fitting &fitting, double aim) -> bool
    {
      const std::array<Vertex, 3> arcs = {
          Vertex{.Position = fitting.Stretch->Start.Offset, .Bulge = fitting.Pair->First},
          Vertex{.Position = fitting.Pair->Joint, .Bulge = fitting.Pair->Last},
          Vertex{.Position = fitting.Stretch->End.Offset},
      };
      const Bend first = Bent(std::span{arcs}.first(2), 0);
      const Bend last  = Bent(std::span{arcs}.last(2), 0);
      double     worst = 0.0;
      for (const Point &point : fitting.Inner)
      {
        worst = std::max(worst, std::min(Apart(first, point), Apart(last, point)));
      }
      return std::max(std::abs(first.Sweep), std::abs(last.Sweep)) <= std::numbers::pi / 2.0
          && worst                                                 <= aim;
    }

    /**
     *  What laying a run reads and where its arcs go.
     */
    struct Paving
    {
      const Parallel::Run    *Source = nullptr; /**< The run.                           */
      const Parallel::Laying *Asked  = nullptr; /**< The distance, accuracy and budget. */
      Layout                 *Into   = nullptr; /**< Where the arcs go.                 */
    };

    /**
     *  A stretch laid with a pair of arcs, or halved and each half laid, until the pair lies within
     *  the accuracy of the offset; a chord where it is halved as deep as it goes.
     *
     *  @param laying  The run, the distance and where the arcs go.
     *  @param stretch The stretch.
     */
    void Lay(const Paving &laying, const Tract &stretch)
    {
      Layout &into = *laying.Into;
      if (into.Spent || into.Pieces.size() + 2 > laying.Asked->Budget)
      {
        into.Spent = true;
        return;
      }
      const double               distance = laying.Asked->Distance;
      const double               aim      = laying.Asked->Accuracy / 2.0;
      const double               quarter  = (stretch.To - stretch.From) / 4.0;
      const double               middle   = stretch.From + (2.0 * quarter);
      const Probe                half     = Probed(*laying.Source, middle, distance);
      const std::array<Point, 3> inner    = {
          Probed(*laying.Source, stretch.From + quarter, distance).Offset,
          half.Offset,
          Probed(*laying.Source, stretch.To - quarter, distance).Offset,
      };
      Box spread = Grown(Grown(Box{}, stretch.Start.Offset), stretch.End.Offset);
      for (const Point &point : inner)
      {
        spread = Grown(spread, point);
      }
      if (Largest(Size(spread)) <= laying.Asked->Accuracy)
      {
        into.Pieces.push_back(Piece(stretch.Start.Offset, stretch.End.Offset, 0.0));
        return;
      }
      const Vector chord   = stretch.End.Offset - stretch.Start.Offset;
      const double ahead   = Dot(chord, stretch.Start.Way);
      const double after   = Dot(chord, stretch.End.Way);
      const bool   forward = (ahead < 0.0) == (after < 0.0) && !Zero(ahead) && !Zero(after);
      const Biarc  pair    = forward ? Paired(stretch, ahead < 0.0 ? -1.0 : 1.0) : Biarc{};
      const bool fits = forward && Fits({.Stretch = &stretch, .Pair = &pair, .Inner = inner}, aim);
      if (fits)
      {
        into.Pieces.push_back(Piece(stretch.Start.Offset, pair.Joint, pair.First));
        into.Pieces.push_back(Piece(pair.Joint, stretch.End.Offset, pair.Last));
        return;
      }
      Lay(laying, {.From = stretch.From, .To = middle, .Start = stretch.Start, .End = half});
      Lay(laying, {.From = middle, .To = stretch.To, .Start = half, .End = stretch.End});
    }
  } // namespace

  auto Sampled(const Run &run, double along) -> Kernel::Sample
  {
    if (const Oval *oval = std::get_if<Oval>(&run))
    {
      const double angle  = oval->Start + (along * oval->Sweep);
      const Vector way    = Trig::Way(angle);
      const double cosine = way.X;
      const double sine   = way.Y;
      return {.At  = oval->Map.Origin + (oval->Map.X * cosine) + (oval->Map.Y * sine),
              .Way = ((oval->Map.Y * cosine) - (oval->Map.X * sine)) * oval->Sweep};
    }
    return Kernel::Sampled(std::get<Kernel::Piece>(run), along);
  }

  auto Laid(std::span<const Run> runs, const Parallel::Laying &laying) -> Layout
  {
    Layout layout;
    for (const Run &run : runs)
    {
      Lay({.Source = &run, .Asked = &laying, .Into = &layout},
          {.From  = 0.0,
           .To    = 1.0,
           .Start = Probed(run, 0.0, laying.Distance),
           .End   = Probed(run, 1.0, laying.Distance)});
    }
    return layout;
  }
} // namespace Geoxy::Parallel
