# Offsetting

An offset is the curve at a fixed distance from another. geoxy gives it two ways: to one side of a curve, as the OFFSET command does, or as a region grown or shrunk by a distance.

## One side of a curve {#side}

`Offset` takes the distance signed: positive is to the right of the way the curve runs, which is outside a loop that runs counterclockwise. `Gaps` says how the gap that opens at an outer corner is closed: `Extend` runs each piece on to a sharp corner, `Fillet` rounds it, `Chamfer` cuts it straight.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
  const Geoxy::Polyline  path      = {.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                                   {.Position = {.X = 10, .Y = 0}},
                                                   {.Position = {.X = 10, .Y = 6}}},
                                      .Closed   = false};

  for (const Geoxy::Gap gaps : {Geoxy::Gap::Extend, Geoxy::Gap::Fillet, Geoxy::Gap::Chamfer})
  {
    const auto offset = Geoxy::Offset(path, {.Distance = 1, .Gaps = gaps}, tolerance);
    for (const Geoxy::Polyline &piece : *offset)
    {
      std::print("{} vertices:", piece.Vertices.size());
      for (const Geoxy::Vertex &vertex : piece.Vertices)
      {
        std::print(" ({:.4g}, {:.4g})", vertex.Position.X, vertex.Position.Y);
      }
      std::println("");
    }
  }
}
```

```text [Output]
3 vertices: (0, -1) (11, -1) (11, 6)
4 vertices: (0, -1) (10, -1) (11, 0) (11, 6)
4 vertices: (0, -1) (10.41, -1) (11, -0.4142) (11, 6)
```

An offset comes back as a list: an inward offset can pinch a shape into several pieces, or leave nothing at all.

## A shape grown or shrunk {#buffer}

`Buffer` treats a closed polyline as the area it fills. `Outward` is that area with everything within the distance of it, `Inward` the area less everything within the distance of its edge, and `Both` the buffer alone, both its edges. The distance is always positive; the side says the direction, whichever way the polyline runs.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
  const Geoxy::Polyline  plate     = {.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                                   {.Position = {.X = 10, .Y = 0}},
                                                   {.Position = {.X = 10, .Y = 6}},
                                                   {.Position = {.X = 0, .Y = 6}}},
                                      .Closed   = true};

  const auto area = [](const std::vector<Geoxy::Polyline> &loops)
  {
    double sum = 0;
    for (const Geoxy::Polyline &loop : loops)
    {
      sum += Geoxy::Area(loop);
    }
    return sum;
  };

  const auto grown  = Geoxy::Buffer(plate, {.Distance = 1, .Sides = Geoxy::Side::Outward}, tolerance);
  const auto shrunk = Geoxy::Buffer(plate, {.Distance = 1, .Sides = Geoxy::Side::Inward}, tolerance);
  const auto gone   = Geoxy::Buffer(plate, {.Distance = 4, .Sides = Geoxy::Side::Inward}, tolerance);

  std::println("the plate:      {}", Geoxy::Area(plate));
  std::println("grown by 1:     {:.6f} in {} loop", area(*grown), grown->size());
  std::println("shrunk by 1:    {} in {} loop", area(*shrunk), shrunk->size());
  std::println("shrunk by 4:    {} loops", gone->size());
}
```

```text [Output]
the plate:      60
grown by 1:     95.141593 in 1 loop
shrunk by 1:    32 in 1 loop
shrunk by 4:    0 loops
```

The loops come back with the region on their left: counterclockwise round it, clockwise round a hole.

## Open shapes, ends and self-crossing {#open}

An open polyline has no inside, so `Buffer` gives the outline of everything within the distance of it. `Ends` says how it closes round the two ends: `Round`, `Flat` or `Square`. For a closed polyline that crosses itself, `Fill` says what it fills: `NonZero` every point it winds about, `EvenOdd` those it winds about an odd number of times.

## Other curves {#others}

A line, a ray, an endless line and an arc offset to their own kind: `Offset(line, distance)` is a `Line`. An ellipse and a spline have no exact offset of their own kind, so theirs come back as polylines of arcs within `Offsetting::Approximation`.

## Which to use {#which}

Use `Offset` when the user picks a side of a curve. Use `Buffer` when the question is about area: a wall thickness, a clearance, a tool path round a pocket.
