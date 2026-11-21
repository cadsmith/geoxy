/***************************************************************************************************
 *  @file      source/belt/make.cpp
 *  @brief     The curves a band's edge can run along, made: each piece moved to its sides, and the
 *             carriers they lie on.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "belt/make.hpp"

#include "belt/piece.hpp"
#include "belt/work.hpp"
#include "scratch.hpp"

#include "geoxy/box.hpp"
#include "geoxy/offset.hpp"
#include "geoxy/point.hpp"
#include "geoxy/polyline.hpp"
#include "geoxy/segment.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace Geoxy::Belt
{
  namespace
  {
    /**
     *  The bulge under which an arc is measured from its ends, its radius and the way it bends, not
     *  from its center: there the center lies more than 2,500 chords off, and what is measured from
     *  it is coarser than a ten million millionth of the chord.
     */
    constexpr double Shoal = 1e-4;
  } // namespace

  void Shaped(Work &work)
  {
    const std::size_t count = work.Vertices.size();
    /**
     *  Room for all that is made from the pieces, taken once: two moved curves and a circle a
     *  piece, a mitre's three more, and a cut at each end of every curve with as many again.
     */
    work.Shape.reserve(count);
    work.Rights.reserve(count);
    work.Lefts.reserve(count);
    work.Curves.reserve((6 * count) + 8);
    work.Ways.reserve((6 * count) + 8);
    work.Cuts.reserve((24 * count) + 32);
    work.Feet.reserve((2 * count) + 2);
    work.Wedges.reserve(count + 2);
    work.Stretches.reserve((12 * count) + 16);
    for (std::size_t index = 0; index < (work.Closed ? count : count - 1); ++index)
    {
      const Point  start = work.Vertices[index].Position;
      const Point  end   = work.Vertices[(index + 1) % count].Position;
      const double bulge = work.Vertices[index].Bulge;
      if (!(std::abs(bulge) >= Straight))
      {
        work.Shape.push_back({.Start = start, .End = end, .Source = work.Shape.size()});
        continue;
      }
      const double factor  = (1.0 - (bulge * bulge)) / (4.0 * bulge);
      const Vector chord   = end - start;
      const Point  center  = {.X = ((start.X + end.X) / 2.0) - (chord.Y * factor),
                              .Y = ((start.Y + end.Y) / 2.0) + (chord.X * factor)};
      const bool   shallow = std::abs(bulge) < Shoal;
      work.Shape.push_back(
          Arced({.Start  = start,
                 .End    = end,
                 .Center = center,
                 .Radius = shallow ? Long(chord) * (1.0 + (bulge * bulge)) / (4.0 * std::abs(bulge))
                                   : Long(start - center),
                 .Left   = bulge > 0.0,
                 .Source = work.Shape.size(),
                 .Shallow = shallow}));
    }
  }

  void Moved(Work &work)
  {
    const double distance = work.Asked.Distance;
    for (const Piece &piece : work.Shape)
    {
      Piece right        = Sided(piece, -distance);
      right.Spared[0]    = piece.Source;
      right.Flank        = 1;
      const bool righted = work.Keeps >= 0 && (!right.Round || right.Radius > 0.0);
      work.Rights.push_back(righted ? Put(work, right, -distance) : None);
      Piece left        = Sided(piece, distance);
      left.Spared[0]    = piece.Source;
      left.Flank        = -1;
      const bool lefted = work.Keeps <= 0 && (!left.Round || left.Radius > 0.0);
      work.Lefts.push_back(lefted ? Put(work, Reversed(left), distance) : None);
    }
  }

  void Hinged(Work &work, const Hinge &hinge)
  {
    const double distance = work.Asked.Distance;
    Piece        curve    = Arced({.Start  = hinge.About + (hinge.Sides[0] * distance),
                                   .End    = hinge.About + (hinge.Sides[1] * distance),
                                   .Center = hinge.About,
                                   .Radius = distance,
                                   .Left   = true});
    curve.Joined          = hinge.Joined;
    const auto lined      = [&work](std::size_t beside, bool line) -> std::size_t
    {
      return beside != None && work.Shape[beside].Round != line ? beside : None;
    };
    curve.Spared = {lined(hinge.Beside[0], true), lined(hinge.Beside[1], true)};
    curve.Beside = {lined(hinge.Beside[0], false), lined(hinge.Beside[1], false)};
    curve.Pivot  = hinge.About;
    curve.Flank  = hinge.Flank;
    /**
     *  A turn so slight that the circle's two ends are one point leaves no circle to put.
     */
    const Vector span = curve.End - curve.Start;
    if (!(Dot(span, span) > 0.0))
    {
      return;
    }
    const std::size_t put = Put(work, curve, 0.0);
    /**
     *  Each of its ends is the vertex moved the distance square to the piece that ends there,
     *  whether or not that piece's moved curve is put: which pieces, and at which end of each, is
     *  kept.
     */
    work.Feet.push_back({static_cast<std::uint32_t>(put),
                         static_cast<std::uint32_t>((2 * hinge.Ended[0]) + hinge.Tails[0]),
                         static_cast<std::uint32_t>((2 * hinge.Ended[1]) + hinge.Tails[1])});
  }

  void Capped(Work &work)
  {
    const Vector start = Aimed(Heading(work.Shape.front(), work.Shape.front().Start));
    const Vector end   = Aimed(Heading(work.Shape.back(), work.Shape.back().End));
    Hinged(work,
           {.About  = work.Shape.front().Start,
            .Sides  = {Vector{.X = -start.Y, .Y = start.X}, Vector{.X = start.Y, .Y = -start.X}},
            .Beside = {0, None},
            .Joined = {work.Lefts.front(), work.Rights.front()},
            .Ended  = {0, 0},
            .Tails  = {0, 0}});
    Hinged(work, {.About  = work.Shape.back().End,
                  .Sides  = {Vector{.X = end.Y, .Y = -end.X}, Vector{.X = -end.Y, .Y = end.X}},
                  .Beside = {work.Shape.size() - 1, None},
                  .Joined = {work.Rights.back(), work.Lefts.back()},
                  .Ended  = {work.Shape.size() - 1, work.Shape.size() - 1},
                  .Tails  = {1, 1}});
  }

  void Sized(Work &work)
  {
    Box all;
    work.Hulls.reserve(work.Shape.size());
    for (const Piece &piece : work.Shape)
    {
      work.Hulls.push_back(Hulled(piece));
      all = Joined(all, work.Hulls.back());
    }
    work.Scale       = std::max({work.Asked.Distance, std::abs(all.Min.X), std::abs(all.Min.Y),
                                 std::abs(all.Max.X), std::abs(all.Max.Y)});
    work.Reach       = work.Scale * 1e-14;
    work.Pad         = work.Asked.Distance + (work.Scale * 1e-6);
    const Vector pad = {.X = work.Pad, .Y = work.Pad};
    work.All         = {.Min = all.Min - pad, .Max = all.Max + pad};
  }
} // namespace Geoxy::Belt
