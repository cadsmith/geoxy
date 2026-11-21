# Curves

A curve in geoxy is a small value you fill in yourself. There are nine kinds, and `Curve` is a `std::variant` that holds any one of them.

## The nine kinds {#kinds}

| Kind | Holds | Parameter runs |
| :--- | :--- | :--- |
| `Line` | `Start`, `End` | 0 at the start to 1 at the end |
| `Ray` | `Start`, `Way` | 0 at the start, without end |
| `Xline` | `Through`, `Way` | without end both ways, 0 at `Through` |
| `Arc` | `Center`, `Radius`, `Start`, `Sweep` | 0 to 1 across the sweep |
| `Ellipse` | `Map`, `Start`, `Sweep` | 0 to 1 across the sweep |
| `Polyline` | `Vertices`, `Closed` | one unit a segment |
| `Spline` | `Degree`, `Knots`, `Poles`, `Weights` | its knots' span |
| `Bezier::Quadratic` | `Start`, `Control`, `End` | 0 to 1 |
| `Bezier::Cubic` | `Start`, `Out`, `In`, `End` | 0 to 1 |

A circle is an arc swept a whole turn: `Circle(center, radius)` makes one. An ellipse is the unit circle carried by a `Matrix`, so its axes need not be square to each other. A spline with `Weights` is rational (a NURBS curve); with none, every weight is 1.

## Angles {#angles}

An `Angle` is made by `Degrees` or `Radians` and remembers which, so that 90 degrees stays exactly 90 and a quarter turn's sine is exactly 1:

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Arc quarter = {.Center = {.X = 0, .Y = 0},
                              .Radius = 2,
                              .Start  = Geoxy::Degrees(0),
                              .Sweep  = Geoxy::Degrees(90)};
  const Geoxy::Point end = Geoxy::At(quarter, 1.0);
  std::println("the arc ends at ({}, {})", end.X, end.Y);
  std::println("its sweep is {} radians", quarter.Sweep.InRadians());
}
```

```text [Output]
the arc ends at (0, 2)
its sweep is 1.5707963267948966 radians
```

A positive sweep runs counterclockwise, a negative one clockwise.

## Polylines and bulges {#polylines}

A polyline is a list of vertices. Each vertex has a `Bulge` for the segment that leaves it: 0 for a straight segment, otherwise the tangent of a quarter of the arc's sweep, positive when the arc bends to the left. A bulge of 1 is a half circle. This is the bulge DXF and DWG files hold, so a polyline read from a drawing goes in unchanged.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  // A slot: two straight sides joined by two half circles.
  const Geoxy::Polyline slot = {.Vertices = {{.Position = {.X = 0, .Y = 0}, .Bulge = 0},
                                             {.Position = {.X = 6, .Y = 0}, .Bulge = 1},
                                             {.Position = {.X = 6, .Y = 2}, .Bulge = 0},
                                             {.Position = {.X = 0, .Y = 2}, .Bulge = 1}},
                                .Closed   = true};
  std::println("{} segments", Geoxy::Count(slot));
  std::println("length {}", Geoxy::Length(slot));
  std::println("area   {}", Geoxy::Area(slot));
}
```

```text [Output]
4 segments
length 18.283185307179586
area   15.141592653589793
```

`Count` and `Nth` read a polyline's segments one at a time, each as a `Segment`: a `Line` or an `Arc`.

## Places on a curve {#places}

A place on a curve is told by its parameter, called `along` everywhere. `At(curve, along)` gives the point there, `Derivative` the velocity, `Tangent` and `Normal` the unit directions, `Curvature` how sharply it bends. `Domain(curve)` gives the range the parameter runs over, and `Ends(curve)` its two end places, where it has them.

Many answers come back as a `Place`: the point and the parameter together, so that the next call can go on from it.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Curve circle = Geoxy::Circle({.X = 0, .Y = 0}, 5);

  const Geoxy::Range domain = Geoxy::Domain(circle);
  std::println("the parameter runs from {} to {}", domain.Start, domain.End);

  const Geoxy::Place place = Geoxy::Nearest(circle, {.X = 0, .Y = 9});
  std::println("nearest (0, 9): ({}, {}) at {}", place.At.X, place.At.Y, place.Along);

  const auto tangent = Geoxy::Tangent(circle, place.Along);
  std::println("the circle runs ({}, {}) there", tangent->X, tangent->Y);
}
```

```text [Output]
the parameter runs from 0 to 1
nearest (0, 9): (0, 5) at 0.25
the circle runs (-1, 0) there
```

## One kind or any kind {#variant}

Every operation has an overload for `Curve`, and most have one for each kind as well. Calling the kind's own overload skips the variant and sometimes returns a narrower type: `Offset` of a `Line` is a `Line`, of a `Curve` a list of curves. To act on what a `Curve` holds, use `std::visit` or `std::get_if` as with any variant.
