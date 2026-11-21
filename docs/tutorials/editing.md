# Editing curves

The editing commands of a CAD program are short compositions of a few calls. Each returns a new curve; nothing is changed in place.

## Trim, extend and break {#trim}

All three work on a place along the curve. Find the place with `Intersect` or `Nearest`, then cut there.

- `Trim(curve, along, keep)` cuts at a place and keeps the part toward the start or the end.
- `Extend(curve, along)` runs a curve on to a place past one of its ends.
- `Split(curve, alongs)` cuts at several places and gives every part.
- `Part(curve, range)` gives the stretch between two places.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
  const Geoxy::Curve     line      = Geoxy::Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 10, .Y = 0}};
  const Geoxy::Curve     edge      = Geoxy::Line{.Start = {.X = 4, .Y = -5}, .End = {.X = 4, .Y = 5}};
  const Geoxy::Curve     far       = Geoxy::Line{.Start = {.X = 14, .Y = -5}, .End = {.X = 14, .Y = 5}};

  const auto show = [](const char *what, const Geoxy::Curve &curve)
  {
    const Geoxy::Line &got = std::get<Geoxy::Line>(curve);
    std::println("{}: ({}, {}) to ({}, {})", what, got.Start.X, got.Start.Y, got.End.X, got.End.Y);
  };

  // Trim: cut the line where the edge crosses it and keep its start.
  const double cut = Geoxy::Intersect(line, edge, tolerance).Hits.front().First.Along;
  show("trimmed", *Geoxy::Trim(line, cut, Geoxy::Keep::Start));

  // Extend: the far edge is past the line's end, so meet it with the line run on without end.
  const Geoxy::Line &segment = std::get<Geoxy::Line>(line);
  const Geoxy::Curve endless = Geoxy::Xline{.Through = segment.Start, .Way = segment.End - segment.Start};
  const double       reach   = Geoxy::Intersect(endless, far, tolerance).Hits.front().First.Along;
  show("extended", *Geoxy::Extend(line, reach));

  // Break: cut at two places and count the parts.
  const std::array places = {0.25, 0.75};
  std::println("broken into {} parts", Geoxy::Split(line, places).size());
}
```

```text [Output]
trimmed: (0, 0) to (4, 0)
extended: (0, 0) to (14, 0)
broken into 3 parts
```

An `Xline` through a line's start with the line's own direction has the same parameter as the line, which is what makes the place found on it usable on the line.

## Join {#join}

`Join` joins curves that touch end to end into one curve: a polyline when they are lines, arcs and polylines, a spline otherwise. It turns each curve round as needed.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};
  const std::array<Geoxy::Curve, 3> pieces = {
      Geoxy::Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 4, .Y = 0}},
      Geoxy::Arc{.Center = {.X = 4, .Y = 1}, .Radius = 1, .Start = Geoxy::Degrees(-90), .Sweep = Geoxy::Degrees(180)},
      Geoxy::Line{.Start = {.X = 0, .Y = 2}, .End = {.X = 4, .Y = 2}},
  };
  const auto joined = Geoxy::Join(pieces, tolerance);
  const Geoxy::Polyline &path = std::get<Geoxy::Polyline>(*joined);
  for (const Geoxy::Vertex &vertex : path.Vertices)
  {
    std::println("({}, {}) bulge {:g}", vertex.Position.X, vertex.Position.Y, vertex.Bulge + 0.0);
  }
}
```

```text [Output]
(0, 0) bulge 0
(4, 0) bulge 1
(4, 2) bulge 0
(0, 2) bulge 0
```

## Fillet and chamfer {#fillet}

`Fillet` joins two curves with an arc of a radius tangent to both, and cuts or runs each on to where the arc touches it. The two picks say which side of each curve to keep, as the points clicked do in a CAD program. `Chamfer` does the same with a straight edge.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Curve across = Geoxy::Line{.Start = {.X = 0, .Y = 0}, .End = {.X = 10, .Y = 0}};
  const Geoxy::Curve upward = Geoxy::Line{.Start = {.X = 10, .Y = 0}, .End = {.X = 10, .Y = 8}};

  const auto corner = Geoxy::Fillet(across, upward, {.Radius = 2, .Picks = {Geoxy::Point{.X = 5, .Y = 0}, Geoxy::Point{.X = 10, .Y = 4}}});
  if (!corner)
  {
    return 1;
  }
  std::println("arc about ({}, {}) of radius {}", corner->Round->Center.X, corner->Round->Center.Y, corner->Round->Radius);
  std::println("touching at ({}, {}) and ({}, {})", corner->Touches[0].X, corner->Touches[0].Y, corner->Touches[1].X, corner->Touches[1].Y);
}
```

```text [Output]
arc about (8, 2) of radius 2
touching at (8, 0) and (10, 2)
```

`Fillet(polyline, radius, tolerance)` and `Chamfer(polyline, bevel, tolerance)` round or cut every corner of a polyline at once, and list the corners where it did not fit.

## Reverse {#reverse}

`Reversed` gives the same curve run the other way.
