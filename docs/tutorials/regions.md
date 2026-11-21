# Regions

A `Region` is an area: a list of `Loop`s, each a list of curves that close. Outer loops run counterclockwise and holes clockwise.

## The area around a pick {#boundary}

HATCH and BOUNDARY start from a click: which closed area is the point in? An `Arrangement` cuts a set of curves where they meet, and `Boundary` gives the face around a point, with its islands as holes.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  // Four walls that overrun their corners, and a column standing in the room.
  const std::array<Geoxy::Curve, 5> curves = {
      Geoxy::Line{.Start = {.X = -1, .Y = 0}, .End = {.X = 11, .Y = 0}},
      Geoxy::Line{.Start = {.X = 10, .Y = -1}, .End = {.X = 10, .Y = 7}},
      Geoxy::Line{.Start = {.X = 11, .Y = 6}, .End = {.X = -1, .Y = 6}},
      Geoxy::Line{.Start = {.X = 0, .Y = 7}, .End = {.X = 0, .Y = -1}},
      Geoxy::Circle({.X = 5, .Y = 3}, 1),
  };
  const Geoxy::Tolerance   tolerance = {.Length = 1e-9, .Radians = 1e-9};
  const Geoxy::Arrangement arrangement{curves, {.Within = tolerance, .Gap = 0}};

  const auto room = Geoxy::Boundary(arrangement, {.X = 2, .Y = 2});
  std::println("{} loops", room->Loops.size());
  std::println("area {:.6f}", Geoxy::Area(*room, Geoxy::Accuracy{}));
}
```

```text [Output]
2 loops
area 56.858407
```

`Closing::Gap` closes small gaps between loose ends, as HPGAPTOL does.

## Is a point inside? {#locate}

`Locate(region, point, tolerance)` answers `Inside`, `Outside` or `On` the edge.

## Union, difference, intersection {#booleans}

`Union`, `Difference` and `Intersection` combine two regions of lines, arcs and polylines. Arcs stay arcs: nothing is turned into chords.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
  const Geoxy::Accuracy  best      = {};

  const auto square = [](double left, double low, double size) -> Geoxy::Region
  {
    const Geoxy::Polyline outline = {.Vertices = {{.Position = {.X = left, .Y = low}},
                                                  {.Position = {.X = left + size, .Y = low}},
                                                  {.Position = {.X = left + size, .Y = low + size}},
                                                  {.Position = {.X = left, .Y = low + size}}},
                                     .Closed   = true};
    return {.Loops = {Geoxy::Loop{.Edges = {Geoxy::Curve{outline}}}}};
  };
  const Geoxy::Region first  = square(0, 0, 4);
  const Geoxy::Region second = square(2, 2, 4);

  std::println("union        {}", Geoxy::Area(*Geoxy::Union(first, second, tolerance), best));
  std::println("difference   {}", Geoxy::Area(*Geoxy::Difference(first, second, tolerance), best));
  std::println("intersection {}", Geoxy::Area(*Geoxy::Intersection(first, second, tolerance), best));
}
```

```text [Output]
union        28
difference   12
intersection 4
```

## Cleaning a loop that crosses itself {#filled}

`Filled` takes a region whose loops may cross themselves or each other and gives the same area as loops that never cross: a bow-tie becomes its two triangles.

## Selecting by a window {#select}

`Inside(curve, box)` says whether a curve lies wholly within a box, and `Intersects(curve, box)` whether it touches the box at all: a window selection and a crossing selection. Both also take a region or a polygon in place of the box.
