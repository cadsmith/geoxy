# Moving and scaling

Move, rotate, scale and mirror are one operation: a `Matrix` applied to a curve by `Transformed`.

## Making a matrix {#matrix}

| Call | Gives |
| :--- | :--- |
| `Translation(vector)` | a move |
| `Rotation(angle)`, `Rotation(angle, about)` | a turn about the origin or a point |
| `Scaling(factor)`, `Scaling(alongX, alongY)` | a scale, even or not |
| `Mirror(first, second)` | a reflection across the line through two points |

Matrices multiply, the right one applied first: `Translation(v) * Rotation(a)` turns, then moves. `Inverse` undoes one where it can be undone.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Curve line = Geoxy::Line{.Start = {.X = 1, .Y = 0}, .End = {.X = 3, .Y = 0}};

  const Geoxy::Matrix turn  = Geoxy::Rotation(Geoxy::Degrees(90));
  const Geoxy::Matrix move  = Geoxy::Translation({.X = 10, .Y = 0});
  const Geoxy::Matrix both  = move * turn;

  const Geoxy::Line turned = std::get<Geoxy::Line>(Geoxy::Transformed(line, both));
  std::println("({}, {}) to ({}, {})", turned.Start.X, turned.Start.Y, turned.End.X, turned.End.Y);
}
```

```text [Output]
(10, 1) to (10, 3)
```

A quarter turn made from `Degrees(90)` is exact: no small error creeps into the coordinates.

## What a curve becomes {#becomes}

Most curves keep their kind. Two do not always: a circle scaled unevenly is no longer a circle, so `Transformed` of an `Arc` returns a `Curve`, which holds an `Arc` where the map keeps circles round and an `Ellipse` where it does not. The same holds for a polyline with arcs.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Arc circle = Geoxy::Circle({.X = 0, .Y = 0}, 1);

  const Geoxy::Curve even   = Geoxy::Transformed(circle, Geoxy::Scaling(3));
  const Geoxy::Curve uneven = Geoxy::Transformed(circle, Geoxy::Scaling(3, 1));

  std::println("scaled evenly:   {}", std::holds_alternative<Geoxy::Arc>(even) ? "an arc" : "an ellipse");
  std::println("scaled unevenly: {}", std::holds_alternative<Geoxy::Arc>(uneven) ? "an arc" : "an ellipse");
}
```

```text [Output]
scaled evenly:   an arc
scaled unevenly: an ellipse
```

## Copies and arrays {#arrays}

A rectangular or polar array is a loop over matrices: apply `Translation` or `Rotation` with a growing step to the same curve.

## Curves drawn in another plane {#space}

DXF and DWG store some entities in their own plane, given by a normal. `Space::Facing(normal)` gives that plane's frame and `Space::Plan(frame, elevation)` the matrix that carries its coordinates into the drawing's.
