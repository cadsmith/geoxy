# Indexing a drawing

Asking every curve of a drawing a question takes time in proportion to the drawing. geoxy has three indexes that ask only the curves that could matter, and give the same answers as asking them all.

| Index | Holds | Good for |
| :--- | :--- | :--- |
| `Tree` | curves, added and removed one at a time | picking, snapping, ray shots while editing |
| `Cloud` | points, added and removed one at a time | endpoint snaps, joining loose ends |
| `Grid` | a fixed set of curves | every intersection of a drawing at once |

## A tree of curves {#tree}

A `Tree` is built from a list of curves or grown with `Insert`. Each curve is told by an index: its place in the list, or what `Insert` returned.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  Geoxy::Tree tree;
  for (int row = 0; row < 100; ++row)
  {
    for (int column = 0; column < 100; ++column)
    {
      const Geoxy::Point corner = {.X = 10.0 * column, .Y = 10.0 * row};
      tree.Insert(Geoxy::Line{.Start = corner, .End = {.X = corner.X + 6, .Y = corner.Y + 3}});
    }
  }
  std::println("{} curves", tree.Size());

  // The three curves nearest the cursor, within a pick aperture of 5 units.
  const Geoxy::Point cursor = {.X = 503, .Y = 498};
  for (const Geoxy::Tree::Close &close : tree.Nearest(cursor, 3, 5.0))
  {
    std::println("curve {} is {:.4f} away, nearest at ({:.4f}, {:.4f})", close.Index, close.Distance,
                 close.Nearest.At.X, close.Nearest.At.Y);
  }

  // A crossing window.
  const Geoxy::Box window = {.Min = {.X = 100, .Y = 100}, .Max = {.X = 125, .Y = 125}};
  std::println("{} curves cross the window, {} lie wholly inside", tree.Crossing(window).size(),
               tree.Window(window).size());
}
```

```text [Output]
10000 curves
curve 5050 is 3.1305 away, nearest at (501.6000, 500.8000)
9 curves cross the window, 6 lie wholly inside
```

`Nearest` takes how many curves to return at most and how far away at most; pass `Tree::Anywhere` for no limit. `Cast(ray, tolerance)` gives every place a ray hits the curves, in order along the ray: the first boundary for TRIM and EXTEND, or the edges a hatch line crosses.

`Erase` removes a curve. Its index is then free and is given to the next curve inserted, so do not keep an index after erasing it.

## A cloud of points {#cloud}

A `Cloud` does the same for points: `Within(box)` for those in a box, `Nearest(point, count, within)` for those closest to a point.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const std::array<Geoxy::Point, 5> ends = {Geoxy::Point{.X = 0, .Y = 0}, Geoxy::Point{.X = 4, .Y = 0},
                                            Geoxy::Point{.X = 4, .Y = 3}, Geoxy::Point{.X = 0, .Y = 3},
                                            Geoxy::Point{.X = 9, .Y = 9}};
  const Geoxy::Cloud cloud{ends};

  for (const Geoxy::Cloud::Close &close : cloud.Nearest({.X = 3.6, .Y = 2.9}, 2, Geoxy::Cloud::Anywhere))
  {
    std::println("point {} at ({}, {}) is {:.4f} away", close.Index, cloud[close.Index].X, cloud[close.Index].Y,
                 close.Distance);
  }
}
```

```text [Output]
point 2 at (4, 3) is 0.4123 away
point 1 at (4, 0) is 2.9275 away
```

## Every intersection at once {#grid}

A `Grid` takes all the curves and a tolerance, prepares each curve once, and `IntersectAll` reports every pair that meets. `IntersectWith` meets one more curve against them all.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  // A 30 by 30 lattice of lines: 30 across and 30 upward.
  std::vector<Geoxy::Curve> curves;
  for (int index = 0; index < 30; ++index)
  {
    curves.emplace_back(Geoxy::Line{.Start = {.X = 0, .Y = index + 0.5}, .End = {.X = 30, .Y = index + 0.5}});
    curves.emplace_back(Geoxy::Line{.Start = {.X = index + 0.5, .Y = 0}, .End = {.X = index + 0.5, .Y = 30}});
  }
  const Geoxy::Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
  const Geoxy::Grid      grid{curves, tolerance};

  std::size_t pairs  = 0;
  std::size_t points = 0;
  Geoxy::IntersectAll(grid, tolerance, [&](const Geoxy::Grid::Pair &, const Geoxy::Marking &marking)
                 {
                   ++pairs;
                   points += marking.Marks.size();
                 });
  std::println("{} pairs meet in {} points", pairs, points);
}
```

```text [Output]
900 pairs meet in 900 points
```

The marking handed to the sink lasts only for that call: copy what you keep. `IntersectAll` also runs on several threads when given a `Workers`, with the same pairs in the same order.

## Which to use {#which}

Use a `Tree` for what changes as the user edits and for anything asked at the cursor. Use a `Grid` when a command needs all intersections of a set of curves that will not change during it. A program can hold both over the same curves.
