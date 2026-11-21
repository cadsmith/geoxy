/***************************************************************************************************
 *  @file      source/belt/cut.cpp
 *  @brief     The curves of a band's edge cut wherever two of them meet, each point numbered once.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "belt/cut.hpp"

#include "belt/carrier.hpp"
#include "belt/corner.hpp"
#include "belt/make.hpp"
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
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

namespace Geoxy::Belt
{
  namespace
  {
    /**
     *  How many curves at the most are each tried against every other, with no grid.
     */
    constexpr std::size_t Few = 96;

    /**
     *  A curve cut at each end of another that falls within it, the two along one line or circle.
     *
     *  @param marks The curves and the points.
     *  @param curve The curve cut.
     *  @param other The other curve.
     */
    void Overlapped(const Marks &marks, std::size_t curve, std::size_t other)
    {
      Work &work = *marks.Curves;
      for (const bool ends : {false, true})
      {
        const Point end = ends ? work.Curves[other].End : work.Curves[other].Start;
        if (Within(work, curve, end))
        {
          work.Cuts.push_back({.Curve = curve,
                               .Along = Along(work.Curves[curve], end),
                               .Node  = Numbered(marks, end, Sourced(other, ends))});
        }
      }
    }

    /**
     *  Whether two curves are two pieces in a row moved to one side, where the shape runs from the
     *  one into the other with no turn: they touch at the end they share alone.
     *
     *  @param  work   The curves.
     *  @param  first  One curve.
     *  @param  second The other.
     *  @return        Whether they are.
     */
    auto Smooth(const Work &work, std::size_t first, std::size_t second) -> bool
    {
      const Piece &one   = work.Curves[first];
      const Piece &other = work.Curves[second];
      if (one.Source                  == None
          || other.Source             == None
          || one.Bounds               != None
          || other.Bounds             != None
          || (work.Ways[first] > 0.0) != (work.Ways[second] > 0.0))
      {
        return false;
      }
      const std::size_t count = work.Shape.size();
      const std::size_t low   = std::min(one.Source, other.Source);
      const std::size_t high  = std::max(one.Source, other.Source);
      if (high == low + 1)
      {
        return work.Even[high] != 0;
      }
      return work.Closed && low == 0 && high == count - 1 && work.Even[0] != 0;
    }

    /**
     *  Two curves cut where they meet: at each point their lines or circles share that falls within
     *  both, read again in two doubles where they come near touching; and, where they run along one
     *  line or circle, each at the other's ends. A circle about a vertex and a curve it joins touch
     *  at the end they share alone, and are not met; nor are two pieces in a row moved to one side
     *  where the shape runs on with no turn.
     *
     *  @param marks  The curves and the points.
     *  @param first  One curve.
     *  @param second The other.
     */
    void Crossed(const Marks &marks, std::size_t first, std::size_t second)
    {
      Work        &work  = *marks.Curves;
      const Piece &one   = work.Curves[first];
      const Piece &other = work.Curves[second];
      /**
       *  The circles come after the moved pieces among the curves, so of two it is the second.
       */
      /**
       *  Nor are two edges of one mitre or square end, which are a loop that meets itself only
       *  where one edge runs into the next: that point is theirs already, and worked out again from
       *  their lines or circles it comes a few steps of a double off, as another.
       */
      if (other.Joined[0] == first
          || other.Joined[1] == first
          || Smooth(work, first, second)
          || (one.Bounds != None && one.Bounds == other.Bounds))
      {
        return;
      }
      /**
       *  Two that doubles take for one line or circle are read again in two doubles: where those
       *  tell them two, they cross as any two do.
       */
      const bool near  = one.Round == other.Round && Together(work, one, other);
      const bool alike = near && Alike(work, Carrier{Widened(work, first)}, Widened(work, second));
      Hits       hits  = Met(one, other);
      if (hits.Close || (near && !alike))
      {
        const Carrier wide = Widened(work, first);
        hits               = Finer(hits, Met(wide, Widened(work, second)));
      }
      for (const Point point : std::span{hits.At}.first(hits.Count))
      {
        const std::pair<double, bool> here  = Lying(work, first, point);
        const std::pair<double, bool> there = here.second ? Lying(work, second, point) : here;
        if (there.second)
        {
          const std::size_t node = Numbered(marks, point,
                                            {.Curve = static_cast<std::uint32_t>(first),
                                             .Other = static_cast<std::uint32_t>(second)});
          work.Cuts.push_back({.Curve = first, .Along = here.first, .Node = node});
          work.Cuts.push_back({.Curve = second, .Along = there.first, .Node = node});
        }
      }
      if (alike)
      {
        Overlapped(marks, first, second);
        Overlapped(marks, second, first);
      }
    }

    /**
     *  Every two of a few curves whose boxes come together cut where they meet, each tried against
     *  every later one.
     *
     *  @param marks The curves, the points and the curves' boxes.
     *  @param boxes The boxes.
     */
    void Paired(const Marks &marks, std::span<const Box> boxes)
    {
      const Work       &work  = *marks.Curves;
      const std::size_t count = boxes.size();
      /**
       *  Each box's four sides in rows of their own, grown by the reach, so that one curve's box is
       *  tried against all the later ones' in one pass that asks nothing between: which of them it
       *  comes together with is marked, and those are met after.
       */
      List<double>            sides(4 * count);
      const std::span<double> left  = std::span{sides}.subspan(0, count);
      const std::span<double> right = std::span{sides}.subspan(count, count);
      const std::span<double> low   = std::span{sides}.subspan(2 * count, count);
      const std::span<double> high  = std::span{sides}.subspan(3 * count, count);
      for (std::size_t index = 0; index < count; ++index)
      {
        left[index]  = boxes[index].Min.X - work.Reach;
        right[index] = boxes[index].Max.X + work.Reach;
        low[index]   = boxes[index].Min.Y - work.Reach;
        high[index]  = boxes[index].Max.Y + work.Reach;
      }
      List<std::uint32_t> near(count);
      /**
       *  Read through their iterators, which ask nothing of each place, so that the pass is one the
       *  compiler can take several boxes at a time.
       */
      const auto lefts  = left.begin();
      const auto rights = right.begin();
      const auto lows   = low.begin();
      const auto highs  = high.begin();
      const auto nears  = near.begin();
      const auto ends   = static_cast<std::ptrdiff_t>(count);
      for (std::ptrdiff_t first = 0; first < ends; ++first)
      {
        const Box box = boxes[static_cast<std::size_t>(first)];
        for (std::ptrdiff_t second = first + 1; second < ends; ++second)
        {
          nears[second] = static_cast<std::uint32_t>(lefts[second] <= box.Max.X)
                        & static_cast<std::uint32_t>(rights[second] >= box.Min.X)
                        & static_cast<std::uint32_t>(lows[second] <= box.Max.Y)
                        & static_cast<std::uint32_t>(highs[second] >= box.Min.Y);
        }
        for (std::ptrdiff_t second = first + 1; second < ends; ++second)
        {
          if (nears[second] != 0)
          {
            Crossed(marks, static_cast<std::size_t>(first), static_cast<std::size_t>(second));
          }
        }
      }
    }
  } // namespace

  auto Ending(const Work &work, std::size_t curve, bool ends) -> Recipe
  {
    const Piece      &piece = work.Curves[curve];
    const std::size_t side  = ends ? 1 : 0;
    if (curve < work.Told.size() && work.Told[curve].at(side).Kind != Making::Unknown)
    {
      return work.Told[curve].at(side);
    }
    if (piece.Source == None)
    {
      /**
       *  A circle about a vertex: found among those put, which are in order.
       */
      const auto feet = std::ranges::equal_range(
          work.Feet, curve, {},
          [](const std::array<std::uint32_t, 3> &foot) -> std::size_t { return foot[0]; });
      if (feet.empty())
      {
        return {.Kind = Making::Given};
      }
      const std::uint32_t ended = feet.front().at(side + 1);
      return {.Kind = Making::Foot, .First = ended / 2, .Second = ended % 2, .About = piece.Center};
    }
    const Piece &shape   = work.Shape[piece.Source];
    const bool   forward = work.Rights[piece.Source] == curve;
    return {.Kind   = Making::Foot,
            .First  = piece.Source,
            .Second = forward == ends ? 1U : 0U,
            .About  = forward == ends ? shape.End : shape.Start};
  }

  auto Remade(Work &work, const Recipe &by, Point near) -> Point
  {
    if (by.Kind == Making::Given || by.Kind == Making::Held)
    {
      return by.Kind == Making::Held ? by.About : near;
    }
    if (by.Kind == Making::Foot)
    {
      /**
       *  Where the circle of the distance about the vertex meets the line through the vertex square
       *  to the piece's heading there, with no root and no division: the point doubles found is
       *  stepped once onto both, how far it lies off the circle and how far its way from the vertex
       *  leans from square read in two doubles a number.
       */
      const std::array<Twofold::Double, 2> heading = Headway(work, by.First, by.Second == 1);
      const Twofold::Two                   outX    = Twofold::TwoSum(near.X, -by.About.X);
      const Twofold::Two                   outY    = Twofold::TwoSum(near.Y, -by.About.Y);
      const Twofold::Double                offX    = {.High = outX.High, .Low = outX.Low};
      const Twofold::Double                offY    = {.High = outY.High, .Low = outY.Low};
      const Twofold::Two reach = Twofold::TwoProduct(work.Asked.Distance, work.Asked.Distance);
      const double       far   = Twofold::Difference(Twofold::Add(Twofold::Multiply(offX, offX),
                                                                  Twofold::Multiply(offY, offY)),
                                                     {.High = reach.High, .Low = reach.Low})
                                     .High;
      const double       lean =
          Twofold::Add(Twofold::Multiply(heading[0], offX), Twofold::Multiply(heading[1], offY))
              .High;
      const Vector out  = {.X = 2.0 * offX.High, .Y = 2.0 * offY.High};
      const Vector way  = {.X = heading[0].High, .Y = heading[1].High};
      const double turn = Cross(out, way);
      return near
           + Vector{.X = ((lean * out.Y) - (far * way.Y)) / turn,
                    .Y = ((far * way.X) - (lean * out.X)) / turn};
    }
    const Carrier &second = Widened(work, by.Second);
    return Sharpened(Widened(work, by.First), second, near);
  }

  void Met(const Marks &marks)
  {
    Work             &work  = *marks.Curves;
    const std::size_t count = work.Curves.size();
    List<Box>         boxes;
    boxes.reserve(count);
    for (std::size_t index = 0; index < count; ++index)
    {
      boxes.push_back(Hulled(work.Curves[index]));
    }
    if (count <= Few)
    {
      Paired(marks, boxes);
      return;
    }
    Stalls grid;
    grid.Lay(work.All, boxes);
    const auto paired = [&](std::size_t row, std::size_t column) -> void
    {
      const std::span<const std::size_t> held = grid.Held(row, column);
      for (std::size_t one = 0; one < held.size(); ++one)
      {
        for (std::size_t other = one + 1; other < held.size(); ++other)
        {
          const Box &first  = boxes[held[one]];
          const Box &second = boxes[held[other]];
          if (grid.First({held[one], held[other]}, {row, column})
              && !Parted(first, second, work.Reach))
          {
            Crossed(marks, held[one], held[other]);
          }
        }
      }
    };
    for (std::size_t row = 0; row < grid.High(); ++row)
    {
      for (std::size_t column = 0; column < grid.Wide(); ++column)
      {
        paired(row, column);
      }
    }
  }
} // namespace Geoxy::Belt
