# Snapping

Object snaps are a handful of questions about one curve and one point. Each has its own call.

| Snap | Call |
| :--- | :--- |
| Nearest | `Nearest(curve, point)` |
| Perpendicular | `Perpendiculars(curve, point)` |
| Tangent | `Tangents(curve, point)` |
| Endpoint | `Ends(curve)` |
| Midpoint | `PlaceAt(curve, Length(curve) / 2)` |
| Center | the `Center` of an `Arc` |
| Intersection | `Intersect(first, second, tolerance)` |

## Nearest, perpendicular, tangent {#nearest}

`Nearest` gives the one place on the curve closest to the point. `Perpendiculars` gives every place where the line from the point meets the curve squarely, and `Tangents` every place where a line from the point touches it.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Curve circle = Geoxy::Circle({.X = 0, .Y = 0}, 5);
  const Geoxy::Point cursor = {.X = 0, .Y = 13};

  const Geoxy::Place nearest = Geoxy::Nearest(circle, cursor);
  std::println("nearest: ({}, {})", nearest.At.X, nearest.At.Y);

  for (const Geoxy::Place &foot : Geoxy::Perpendiculars(circle, cursor))
  {
    std::println("perpendicular: ({}, {})", foot.At.X, foot.At.Y);
  }
  for (const Geoxy::Place &touch : Geoxy::Tangents(circle, cursor))
  {
    std::println("tangent: ({:.6f}, {:.6f})", touch.At.X, touch.At.Y);
  }
}
```

```text [Output]
nearest: (0, 5)
perpendicular: (0, 5)
perpendicular: (0, -5)
tangent: (4.615385, 1.923077)
tangent: (-4.615385, 1.923077)
```

## Ends and the middle {#ends}

`Ends` gives a curve's start and end places; a closed curve or an endless line has none, so each is optional. The middle by length is `PlaceAt` at half the length, which is not the middle of the parameter on a curve whose speed varies.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Curve     path = Geoxy::Polyline{.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                                              {.Position = {.X = 8, .Y = 0}},
                                                              {.Position = {.X = 8, .Y = 2}}},
                                                .Closed   = false};
  const Geoxy::Accuracy  best = {};
  const Geoxy::Endpoints ends = Geoxy::Ends(path);
  std::println("from ({}, {}) to ({}, {})", ends.Start->At.X, ends.Start->At.Y, ends.End->At.X, ends.End->At.Y);

  const Geoxy::Place middle = Geoxy::PlaceAt(path, Geoxy::Length(path, best) / 2, best);
  std::println("middle by length: ({}, {})", middle.At.X, middle.At.Y);
}
```

```text [Output]
from (0, 0) to (8, 2)
middle by length: (5, 0)
```

## Circles and lines by construction {#construct}

`Circle(first, second)` is the circle on a diameter, `Circle(first, second, third)` the one through three points. `TangentCircles` gives the circles tangent to two lines or arcs with a radius, or to three; `TangentLines` the lines from a point tangent to a circle, or tangent to two circles.

## Among many curves {#many}

Asking every curve of a drawing for its nearest point on each mouse move is too slow. A `Tree` answers the nearest curves to the cursor directly; see [Indexing a drawing](/tutorials/indexing/).
