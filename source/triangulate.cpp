/***************************************************************************************************
 *  @file      source/triangulate.cpp
 *  @brief     A region cut into triangles: holes bridged, ears clipped with exact turns along a Z
 *             curve, as earcut does, and edges flipped to Delaunay's by the exact in-circle test.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "geoxy/triangulate.hpp"

#include "exact/circle.hpp"
#include "real.hpp"

#include "geoxy/orient.hpp"
#include "geoxy/point.hpp"
#include "geoxy/vector.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <utility>
#include <vector>

namespace Geoxy
{
  namespace
  {
    /**
     *  No node.
     */
    constexpr std::uint32_t None = std::numeric_limits<std::uint32_t>::max();

    /**
     *  How many points a region has before its corners are looked up along a Z curve.
     */
    constexpr std::size_t Hashed = 80;

    /**
     *  A corner of a ring: its point, its index among the loops' points, its neighbours on the ring
     *  and along the Z curve, and its place on that curve.
     */
    struct Node
    {
      Point         At    = {};   /**< The point.                           */
      std::uint32_t Index = 0;    /**< Its index among the loops' points.   */
      std::uint32_t Prev  = None; /**< The corner before it on the ring.    */
      std::uint32_t Next  = None; /**< The corner after it.                 */
      std::uint32_t Below = None; /**< The corner before it on the Z curve. */
      std::uint32_t Above = None; /**< The corner after it on the Z curve.  */
      std::uint32_t Z     = 0;    /**< Its place on the Z curve.            */
    };

    /**
     *  The rings being cut: their corners, how the Z curve is laid over them, and the triangles
     *  cut.
     */
    struct Net
    {
      std::vector<Node> Nodes  = {};    /**< Every corner, removed ones too.              */
      Point             Least  = {};    /**< The corner of the box the Z curve starts at. */
      double            Scale  = 0.0;   /**< Z curve steps to a unit of length.           */
      bool              Hashed = false; /**< Whether ears are sought along it.            */
      Geoxy::Mesh       Mesh   = {};    /**< The triangles cut, and whether all.          */
    };

    /**
     *  The way three corners turn.
     *
     *  @param  net     The rings.
     *  @param  corners The three corners.
     *  @return         Their orientation, exactly.
     */
    auto Handed(const Net &net, std::array<std::uint32_t, 3> corners) -> Orientation
    {
      return Orient(net.Nodes[corners[0]].At, net.Nodes[corners[1]].At, net.Nodes[corners[2]].At);
    }

    /**
     *  Whether a corner turns counterclockwise, as its ring runs.
     *
     *  @param  net  The rings.
     *  @param  node The corner.
     *  @return      Whether it does: a corner that turns the other way, or none, is reflex.
     */
    auto Convex(const Net &net, std::uint32_t node) -> bool
    {
      const Node &here = net.Nodes[node];
      return Handed(net, {here.Prev, node, here.Next}) == Orientation::Counterclockwise;
    }

    /**
     *  A corner added after another on a ring.
     *
     *  @param  net   The rings.
     *  @param  point The corner's point and its index.
     *  @param  last  The corner it follows, or none to start a ring.
     *  @return       The corner.
     */
    auto Added(Net &net, std::pair<Point, std::uint32_t> point, std::uint32_t last) -> std::uint32_t
    {
      const auto node = static_cast<std::uint32_t>(net.Nodes.size());
      net.Nodes.push_back({.At = point.first, .Index = point.second});
      Node &added = net.Nodes.back();
      if (last == None)
      {
        added.Prev = node;
        added.Next = node;
        return node;
      }
      Node &before                = net.Nodes[last];
      added.Next                  = before.Next;
      added.Prev                  = last;
      net.Nodes[before.Next].Prev = node;
      before.Next                 = node;
      return node;
    }

    /**
     *  A corner taken off its ring, and off the Z curve.
     *
     *  @param net  The rings.
     *  @param node The corner.
     */
    void Removed(Net &net, std::uint32_t node)
    {
      const Node &gone          = net.Nodes[node];
      net.Nodes[gone.Next].Prev = gone.Prev;
      net.Nodes[gone.Prev].Next = gone.Next;
      if (gone.Below != None)
      {
        net.Nodes[gone.Below].Above = gone.Above;
      }
      if (gone.Above != None)
      {
        net.Nodes[gone.Above].Below = gone.Below;
      }
    }

    /**
     *  A loop made a ring, turning as asked.
     *
     *  @param  net   The rings.
     *  @param  loop  The loop's points.
     *  @param  first The index of its first point among all the loops', and whether the ring is to
     *                run counterclockwise.
     *  @return       A corner of the ring; none for a loop of no point.
     */
    auto Netted(Net &net, std::span<const Point> loop, std::pair<std::uint32_t, bool> first)
        -> std::uint32_t
    {
      double area = 0.0;
      for (std::size_t index = 1; index + 1 < loop.size(); ++index)
      {
        area += Cross(loop[index] - loop[0], loop[index + 1] - loop[0]);
      }
      const bool    forward = (area > 0.0) == first.second;
      std::uint32_t last    = None;
      for (std::size_t step = 0; step < loop.size(); ++step)
      {
        const std::size_t index = forward ? step : loop.size() - 1 - step;
        last = Added(net, {loop[index], first.first + static_cast<std::uint32_t>(index)}, last);
      }
      return last;
    }

    /**
     *  A ring with its corners that turn neither way taken off: points the same as the next, and
     *  points on a line with their neighbours, which a triangle would have no area at.
     *
     *  @param  net   The rings.
     *  @param  start A corner of the ring.
     *  @return       A corner still on it; none where fewer than three are left.
     */
    auto Filtered(Net &net, std::uint32_t start) -> std::uint32_t
    {
      std::uint32_t node = start;
      std::uint32_t stop = start;
      while (net.Nodes[node].Next != net.Nodes[node].Prev)
      {
        const Node &here = net.Nodes[node];
        if (Handed(net, {here.Prev, node, here.Next}) == Orientation::Collinear)
        {
          Removed(net, node);
          node = here.Prev;
          stop = node;
          continue;
        }
        node = here.Next;
        if (node == stop)
        {
          return node;
        }
      }
      return None;
    }

    /**
     *  The leftmost corner of a ring, the lowest of those.
     *
     *  @param  net   The rings.
     *  @param  start A corner of the ring.
     *  @return       The corner.
     */
    auto Leftmost(const Net &net, std::uint32_t start) -> std::uint32_t
    {
      std::uint32_t best = start;
      std::uint32_t node = net.Nodes[start].Next;
      while (node != start)
      {
        const Point here = net.Nodes[node].At;
        const Point most = net.Nodes[best].At;
        if (std::pair{here.X, here.Y} < std::pair{most.X, most.Y})
        {
          best = node;
        }
        node = net.Nodes[node].Next;
      }
      return best;
    }

    /**
     *  Whether a point is within a triangle that turns counterclockwise, or on its sides.
     *
     *  @param  corners The triangle.
     *  @param  point   The point.
     *  @return         Whether it is.
     */
    auto Contained(const std::array<Point, 3> &corners, Point point) -> bool
    {
      return Orient(corners[0], corners[1], point) != Orientation::Clockwise
          && Orient(corners[1], corners[2], point) != Orientation::Clockwise
          && Orient(corners[2], corners[0], point) != Orientation::Clockwise;
    }

    /**
     *  Whether the diagonal from one corner to a point leaves it inside the ring.
     *
     *  @param  net   The rings.
     *  @param  node  The corner.
     *  @param  point The point.
     *  @return       Whether the diagonal starts into the ring's inside.
     */
    auto Visible(const Net &net, std::uint32_t node, Point point) -> bool
    {
      const Node &here   = net.Nodes[node];
      const Point before = net.Nodes[here.Prev].At;
      const Point after  = net.Nodes[here.Next].At;
      if (Convex(net, node))
      {
        return Orient(here.At, point, after)  != Orientation::Counterclockwise
            && Orient(here.At, before, point) != Orientation::Counterclockwise;
      }
      return Orient(here.At, point, before) == Orientation::Clockwise
          || Orient(here.At, after, point)  == Orientation::Clockwise;
    }

    /**
     *  The corner of the outline a hole is bridged to, as earcut finds it: the nearest edge the ray
     *  to the hole's left meets, its end the further left, unless a reflex corner within the
     *  triangle that end makes with the ray sees the hole at a smaller angle.
     *
     *  @param  net   The rings.
     *  @param  hole  The hole's leftmost corner.
     *  @param  outer A corner of the outline.
     *  @return       The corner; none where the ray meets no edge.
     */
    auto Bridge(const Net &net, std::uint32_t hole, std::uint32_t outer) -> std::uint32_t
    {
      const Point   here  = net.Nodes[hole].At;
      double        reach = -std::numeric_limits<double>::infinity();
      std::uint32_t best  = None;
      std::uint32_t node  = outer;
      for (bool first = true; first || node != outer; first = false)
      {
        const Node &edge = net.Nodes[node];
        const Point from = edge.At;
        const Point onto = net.Nodes[edge.Next].At;
        if (onto == here)
        {
          return edge.Next;
        }
        if (here.Y <= from.Y && onto.Y <= here.Y && onto.Y < from.Y)
        {
          const double cross = from.X + ((here.Y - from.Y) * (onto.X - from.X) / (onto.Y - from.Y));
          if (cross <= here.X && cross > reach)
          {
            reach = cross;
            best  = from.X < onto.X ? node : edge.Next;
          }
        }
        node = edge.Next;
      }
      if (best == None)
      {
        return None;
      }
      const Point                corner = net.Nodes[best].At;
      const std::array<Point, 3> seen =
          here.Y < corner.Y ? std::array{here, corner, Point{.X = reach, .Y = here.Y}}
                            : std::array{Point{.X = reach, .Y = here.Y}, corner, here};
      double        least = std::numeric_limits<double>::infinity();
      std::uint32_t found = best;
      node                = best;
      for (bool first = true; first || node != best; first = false)
      {
        const Point point = net.Nodes[node].At;
        if (point.X < here.X
            && corner.X <= point.X
            && Contained(seen, point)
            && Visible(net, node, here))
        {
          const double slope = std::abs(here.Y - point.Y) / (here.X - point.X);
          if (std::pair{slope, -point.X} < std::pair{least, -net.Nodes[found].At.X})
          {
            least = slope;
            found = node;
          }
        }
        node = net.Nodes[node].Next;
      }
      return found;
    }

    /**
     *  A ring split along a diagonal into two, each corner of it doubled: the bridge to a hole.
     *
     *  @param net  The rings.
     *  @param from The diagonal's first corner.
     *  @param into Its second.
     */
    void Bridging(Net &net, std::uint32_t from, std::uint32_t into)
    {
      const std::uint32_t twin   = Added(net, {net.Nodes[from].At, net.Nodes[from].Index}, None);
      const std::uint32_t other  = Added(net, {net.Nodes[into].At, net.Nodes[into].Index}, None);
      const std::uint32_t after  = net.Nodes[from].Next;
      const std::uint32_t before = net.Nodes[into].Prev;
      net.Nodes[from].Next       = into;
      net.Nodes[into].Prev       = from;
      net.Nodes[twin].Next       = after;
      net.Nodes[after].Prev      = twin;
      net.Nodes[other].Next      = twin;
      net.Nodes[twin].Prev       = other;
      net.Nodes[before].Next     = other;
      net.Nodes[other].Prev      = before;
    }

    /**
     *  A point's place along the Z curve over the box: its coordinates scaled to 15 bits each and
     *  their bits interleaved.
     *
     *  @param  net   The rings, with the box.
     *  @param  point The point.
     *  @return       Its place.
     */
    auto Zed(const Net &net, Point point) -> std::uint32_t
    {
      auto x = static_cast<std::uint32_t>((point.X - net.Least.X) * net.Scale);
      auto y = static_cast<std::uint32_t>((point.Y - net.Least.Y) * net.Scale);
      for (const auto &[shift, mask] : {std::pair{8U, 0x00FF00FFU}, std::pair{4U, 0x0F0F0F0FU},
                                        std::pair{2U, 0x33333333U}, std::pair{1U, 0x55555555U}})
      {
        x = (x | (x << shift)) & mask;
        y = (y | (y << shift)) & mask;
      }
      return x | (y << 1U);
    }

    /**
     *  A ring's corners laid along the Z curve, in order.
     *
     *  @param net   The rings.
     *  @param start A corner of the ring.
     */
    void Zordered(Net &net, std::uint32_t start)
    {
      std::vector<std::uint32_t> order;
      std::uint32_t              node = start;
      for (bool first = true; first || node != start; first = false)
      {
        net.Nodes[node].Z = Zed(net, net.Nodes[node].At);
        order.push_back(node);
        node = net.Nodes[node].Next;
      }
      std::ranges::sort(
          order, [&net](std::uint32_t lhs, std::uint32_t rhs) -> bool
          { return std::pair{net.Nodes[lhs].Z, lhs} < std::pair{net.Nodes[rhs].Z, rhs}; });
      for (std::size_t index = 0; index < order.size(); ++index)
      {
        Node &here = net.Nodes[order[index]];
        here.Below = index == 0 ? None : order[index - 1];
        here.Above = index + 1 == order.size() ? None : order[index + 1];
      }
    }

    /**
     *  Whether a corner keeps an ear from being cut: a reflex corner within the ear or on its
     *  sides, neither of the ear's other corners nor at its first corner's point.
     *
     *  @param  net  The rings.
     *  @param  ear  The ear's three corners.
     *  @param  node The corner.
     *  @return      Whether it does.
     */
    auto Blocks(const Net &net, const std::array<std::uint32_t, 3> &ear, std::uint32_t node) -> bool
    {
      const Point point = net.Nodes[node].At;
      const Point first = net.Nodes[ear[0]].At;
      return node != ear[0]
          && node != ear[2]
          && !(point == first)
          && Contained({first, net.Nodes[ear[1]].At, net.Nodes[ear[2]].At}, point)
          && !Convex(net, node);
    }

    /**
     *  Whether a point lies within a box, its sides too.
     *
     *  @param  point The point.
     *  @param  box   The box's least and most corners.
     *  @return       Whether it does.
     */
    auto Boxed(Point point, std::pair<Point, Point> box) -> bool
    {
      return box.first.X <= point.X
          && point.X     <= box.second.X
          && box.first.Y <= point.Y
          && point.Y     <= box.second.Y;
    }

    /**
     *  Whether a corner is an ear: it turns counterclockwise and no reflex corner lies within the
     *  triangle it makes with its neighbours, every other corner asked, or along the Z curve those
     *  within the triangle's box.
     *
     *  @param  net  The rings.
     *  @param  node The corner.
     *  @return      Whether it is.
     */
    auto Ear(const Net &net, std::uint32_t node) -> bool
    {
      const Node                        &here = net.Nodes[node];
      const std::array<std::uint32_t, 3> ear  = {here.Prev, node, here.Next};
      if (!Convex(net, node))
      {
        return false;
      }
      if (!net.Hashed)
      {
        for (std::uint32_t other = net.Nodes[here.Next].Next; other != here.Prev;
             other               = net.Nodes[other].Next)
        {
          if (Blocks(net, ear, other))
          {
            return false;
          }
        }
        return true;
      }
      const Point         first  = net.Nodes[ear[0]].At;
      const Point         second = here.At;
      const Point         third  = net.Nodes[ear[2]].At;
      const Point         least  = {.X = std::min({first.X, second.X, third.X}),
                                    .Y = std::min({first.Y, second.Y, third.Y})};
      const Point         most   = {.X = std::max({first.X, second.X, third.X}),
                                    .Y = std::max({first.Y, second.Y, third.Y})};
      const std::uint32_t low    = Zed(net, least);
      const std::uint32_t high   = Zed(net, most);
      for (std::uint32_t other = here.Below; other != None && net.Nodes[other].Z >= low;
           other               = net.Nodes[other].Below)
      {
        if (Boxed(net.Nodes[other].At, {least, most}) && Blocks(net, ear, other))
        {
          return false;
        }
      }
      for (std::uint32_t other = here.Above; other != None && net.Nodes[other].Z <= high;
           other               = net.Nodes[other].Above)
      {
        if (Boxed(net.Nodes[other].At, {least, most}) && Blocks(net, ear, other))
        {
          return false;
        }
      }
      return true;
    }

    /**
     *  Whether a corner on a line with its neighbours lies between them, the ring running straight
     *  on through it: a corner whose triangles would leave it in the middle of another's side if it
     *  were taken off, where one that turns back, or is one of its neighbours, may go.
     *
     *  @param  net  The rings.
     *  @param  node The corner, on a line with its neighbours.
     *  @return      Whether it lies strictly between them.
     */
    auto Flanked(const Net &net, std::uint32_t node) -> bool
    {
      const Node  &here    = net.Nodes[node];
      const Point  before  = net.Nodes[here.Prev].At;
      const Point  after   = net.Nodes[here.Next].At;
      const bool   upright = Exactly(before.X, after.X);
      const double low     = upright ? std::min(before.Y, after.Y) : std::min(before.X, after.X);
      const double high    = upright ? std::max(before.Y, after.Y) : std::max(before.X, after.X);
      const double along   = upright ? here.At.Y : here.At.X;
      return low < along && along < high;
    }

    /**
     *  A ring cut into ears, one after another, corners that turn neither way taken off as they
     *  come; a ring left where no ear is found all round, which only loops that cross leave.
     *
     *  @param net   The rings, the triangles cut added to its mesh.
     *  @param start A corner of the ring.
     */
    void Eared(Net &net, std::uint32_t start)
    {
      std::uint32_t node = start;
      std::uint32_t stop = start;
      while (net.Nodes[node].Prev != net.Nodes[node].Next)
      {
        const Node &here = net.Nodes[node];
        if (Handed(net, {here.Prev, node, here.Next}) == Orientation::Collinear
            && !Flanked(net, node))
        {
          Removed(net, node);
          node = here.Next;
          stop = node;
          continue;
        }
        if (Ear(net, node))
        {
          net.Mesh.Triangles.push_back(
              {net.Nodes[here.Prev].Index, here.Index, net.Nodes[here.Next].Index});
          Removed(net, node);
          node = net.Nodes[here.Next].Next;
          stop = node;
          continue;
        }
        node = here.Next;
        if (node == stop)
        {
          net.Mesh.Whole = false;
          return;
        }
      }
    }

    /**
     *  A triangle's corner, its index taken round.
     *
     *  @param  triangle The triangle.
     *  @param  corner   The corner's index, from 0 on, taken modulo 3.
     *  @return          The corner.
     */
    auto Vertexed(const std::array<std::uint32_t, 3> &triangle, std::uint32_t corner)
        -> std::uint32_t
    {
      return std::span{triangle}[corner % 3];
    }

    /**
     *  The triangles across each triangle's sides: for side k, from corner k to k + 1, the triangle
     *  that has it the other way round, or none on the region's boundary.
     *
     *  @param  triangles The triangles.
     *  @return           The neighbours, three per triangle.
     */
    auto Neighbours(std::span<const std::array<std::uint32_t, 3>> triangles)
        -> std::vector<std::array<std::uint32_t, 3>>
    {
      std::vector<std::pair<std::uint64_t, std::uint32_t>> sides;
      for (std::size_t triangle = 0; triangle < triangles.size(); ++triangle)
      {
        for (std::uint32_t side = 0; side < 3; ++side)
        {
          const std::uint64_t tail = Vertexed(triangles[triangle], side);
          const std::uint64_t head = Vertexed(triangles[triangle], side + 1);
          sides.emplace_back((tail << 32U) | head, static_cast<std::uint32_t>(triangle));
        }
      }
      std::ranges::sort(sides);
      sides.emplace_back(std::numeric_limits<std::uint64_t>::max(), None);
      std::vector<std::array<std::uint32_t, 3>> across(triangles.size(), {None, None, None});
      for (std::size_t triangle = 0; triangle < triangles.size(); ++triangle)
      {
        for (std::uint32_t side = 0; side < 3; ++side)
        {
          const std::uint64_t tail = Vertexed(triangles[triangle], side);
          const std::uint64_t head = Vertexed(triangles[triangle], side + 1);
          const auto          found =
              std::ranges::lower_bound(sides, std::pair{(head << 32U) | tail, std::uint32_t{0}});
          if (found->first == ((head << 32U) | tail))
          {
            std::span{across[triangle]}[side] = found->second;
          }
        }
      }
      return across;
    }

    /**
     *  The side of a triangle that runs from one corner to another.
     *
     *  @param  triangle The triangle.
     *  @param  from     The side's first corner.
     *  @return          Its index.
     */
    auto Side(const std::array<std::uint32_t, 3> &triangle, std::uint32_t from) -> std::uint32_t
    {
      if (triangle[0] == from)
      {
        return 0U;
      }
      return triangle[1] == from ? 1U : 2U;
    }

    /**
     *  Triangles flipped, side by side, until no corner lies inside the circle of the triangle
     *  across from it: Lawson's flips, each decided by the exact in-circle test, ending at the
     *  constrained Delaunay triangulation, the region's boundary kept.
     *
     *  @param points    Every loop's points, in order.
     *  @param triangles The triangles, flipped in place.
     */
    void Legalized(std::span<const Point>                     points,
                   std::vector<std::array<std::uint32_t, 3>> &triangles)
    {
      std::vector<std::array<std::uint32_t, 3>>            across = Neighbours(triangles);
      std::vector<std::pair<std::uint32_t, std::uint32_t>> stack;
      for (std::uint32_t triangle = 0; triangle < triangles.size(); ++triangle)
      {
        for (std::uint32_t side = 0; side < 3; ++side)
        {
          stack.emplace_back(triangle, side);
        }
      }
      while (!stack.empty())
      {
        const auto [first, side] = stack.back();
        stack.pop_back();
        const std::uint32_t second = Vertexed(across[first], side);
        if (second == None)
        {
          continue;
        }
        std::array<std::uint32_t, 3> &near = triangles[first];
        std::array<std::uint32_t, 3> &far  = triangles[second];
        const std::uint32_t           a    = Vertexed(near, side);
        const std::uint32_t           b    = Vertexed(near, side + 1);
        const std::uint32_t           c    = Vertexed(near, side + 2);
        const std::uint32_t           back = Side(far, b);
        const std::uint32_t           d    = Vertexed(far, back + 2);
        if (Exact::InCircle({.First = points[a], .Second = points[b], .Third = points[c]},
                            points[d])
            <= 0)
        {
          continue;
        }
        const std::uint32_t besideA = Vertexed(across[second], back + 1);
        const std::uint32_t besideD = Vertexed(across[second], back + 2);
        const std::uint32_t besideB = Vertexed(across[first], side + 1);
        const std::uint32_t besideC = Vertexed(across[first], side + 2);
        near                        = {a, d, c};
        far                         = {d, b, c};
        across[first]               = {besideA, second, besideC};
        across[second]              = {besideD, besideB, first};
        if (besideA != None)
        {
          std::span{across[besideA]}[Side(triangles[besideA], d)] = first;
        }
        if (besideB != None)
        {
          std::span{across[besideB]}[Side(triangles[besideB], c)] = second;
        }
        stack.emplace_back(first, 0);
        stack.emplace_back(first, 2);
        stack.emplace_back(second, 0);
        stack.emplace_back(second, 1);
      }
    }
  } // namespace

  auto Triangulate(std::span<const std::vector<Point>> loops, Fill fill) -> Mesh
  {
    Net                        net;
    std::vector<Point>         points;
    std::uint32_t              outer = None;
    std::vector<std::uint32_t> holes;
    for (std::size_t loop = 0; loop < loops.size(); ++loop)
    {
      const auto          first = static_cast<std::uint32_t>(points.size());
      const std::uint32_t ring  = Netted(net, loops[loop], {first, loop == 0});
      points.insert(points.end(), loops[loop].begin(), loops[loop].end());
      const std::uint32_t kept = ring == None ? None : Filtered(net, ring);
      if (loop == 0)
      {
        outer = kept;
      }
      else if (kept != None)
      {
        holes.push_back(Leftmost(net, kept));
      }
    }
    if (outer == None)
    {
      return {};
    }
    std::ranges::sort(holes,
                      [&net](std::uint32_t lhs, std::uint32_t rhs) -> bool
                      {
                        const Point left  = net.Nodes[lhs].At;
                        const Point right = net.Nodes[rhs].At;
                        return std::pair{left.X, left.Y} < std::pair{right.X, right.Y};
                      });
    for (const std::uint32_t hole : holes)
    {
      const std::uint32_t bridge = Bridge(net, hole, outer);
      if (bridge == None)
      {
        net.Mesh.Whole = false;
        continue;
      }
      Bridging(net, bridge, hole);
      outer = bridge;
    }
    if (points.size() > Hashed)
    {
      Point least = points.front();
      Point most  = points.front();
      for (const Point point : points)
      {
        least = {.X = std::min(least.X, point.X), .Y = std::min(least.Y, point.Y)};
        most  = {.X = std::max(most.X, point.X), .Y = std::max(most.Y, point.Y)};
      }
      const double size = std::max(most.X - least.X, most.Y - least.Y);
      net.Least         = least;
      net.Scale         = 32767.0 / size;
      net.Hashed        = true;
      Zordered(net, outer);
    }
    Eared(net, outer);
    if (fill == Fill::Delaunay)
    {
      Legalized(points, net.Mesh.Triangles);
    }
    return std::move(net.Mesh);
  }
} // namespace Geoxy
