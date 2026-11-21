# Drawing curves

A screen draws straight lines and triangles. These calls turn curves into them, as finely as the view needs and no finer.

## Chords for the screen {#flatten}

`Flatten` gives the points of a chain of chords that stays within an `Accuracy` of the curve. Pass the size of a pixel in drawing units and a curve gets just enough points to look right at the current zoom.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Curve arc = Geoxy::Arc{.Center = {.X = 0, .Y = 0}, .Radius = 50, .Start = Geoxy::Degrees(0), .Sweep = Geoxy::Degrees(90)};

  for (const double pixel : {1.0, 0.1, 0.01})
  {
    const std::vector<Geoxy::Point> points = Geoxy::Flatten(arc, Geoxy::Accuracy{.Error = pixel});
    std::println("a pixel of {:>4}: {} points", pixel, points.size());
  }
}
```

```text [Output]
a pixel of    1: 5 points
a pixel of  0.1: 14 points
a pixel of 0.01: 41 points
```

## Only what is in view {#clip}

`Clip(curve, box)` gives the stretches of a curve's parameter that lie inside a box, so that a curve far larger than the view, or one without end, is drawn only where it shows. Take each stretch with `Part` and flatten that.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Curve endless = Geoxy::Xline{.Through = {.X = 0, .Y = 0}, .Way = {.X = 1, .Y = 1}};
  const Geoxy::Box   view    = {.Min = {.X = 2, .Y = 0}, .Max = {.X = 10, .Y = 6}};

  for (const Geoxy::Range &stretch : Geoxy::Clip(endless, view))
  {
    const Geoxy::Point from = Geoxy::At(endless, stretch.Start);
    const Geoxy::Point to   = Geoxy::At(endless, stretch.End);
    std::println("visible from ({}, {}) to ({}, {})", from.X, from.Y, to.X, to.Y);
  }
}
```

```text [Output]
visible from (2, 2) to (6, 6)
```

## Triangles to fill {#triangulate}

`Triangulate` cuts a polygon with holes into triangles for a solid fill. It takes the outline first and then the holes, each a list of points, and returns the triangles as indices into those points counted through all the loops in order.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const std::vector<std::vector<Geoxy::Point>> loops = {
      {{.X = 0, .Y = 0}, {.X = 10, .Y = 0}, {.X = 10, .Y = 10}, {.X = 0, .Y = 10}},
      {{.X = 4, .Y = 4}, {.X = 6, .Y = 4}, {.X = 6, .Y = 6}, {.X = 4, .Y = 6}},
  };
  const Geoxy::Mesh mesh = Geoxy::Triangulate(loops, Geoxy::Fill::Ears);
  std::println("{} triangles, {}", mesh.Triangles.size(), mesh.Whole ? "covering it all" : "incomplete");
}
```

```text [Output]
8 triangles, covering it all
```

`Fill::Delaunay` gives better shaped triangles for the same area.

## Linetypes {#dashes}

`Dashes` lays a dash pattern along a curve by length. A `Pattern` holds the lengths of its elements as a DXF LTYPE does: positive for a dash, negative for a gap, 0 for a dot. Each dash comes back as a curve of the same kind, so a dashed arc is made of arcs.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Curve   line   = Geoxy::Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 10, .Y = 0}};
  const Geoxy::Pattern dashed = {.Elements = {2, -1}};

  for (const Geoxy::Curve &dash : Geoxy::Dashes(line, dashed, Geoxy::Accuracy{}))
  {
    const Geoxy::Line &piece = std::get<Geoxy::Line>(dash);
    std::println("{} to {}", piece.Start.X, piece.End.X);
  }
}
```

```text [Output]
0 to 2
3 to 5
6 to 8
9 to 10
```

## Wide polylines {#outline}

`Outline` gives the filled shape of a polyline drawn with a width, one closed polyline per segment, with the width able to change along each segment as LWPOLYLINE allows.
