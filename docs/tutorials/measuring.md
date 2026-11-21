# Measuring

Lengths, areas, extents and centroids, each a function of a curve.

## Length {#length}

`Length` is exact for lines, arcs and polylines. For an ellipse, a spline or a Bézier curve it takes an `Accuracy`. `Length(curve, range)` measures a stretch of the parameter, and `PlaceAt(curve, distance)` goes the other way: the place a distance along the curve from its start.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Arc half = {.Center = {.X = 0, .Y = 0}, .Radius = 1, .Start = Geoxy::Degrees(0), .Sweep = Geoxy::Degrees(180)};
  std::println("half a unit circle: {}", Geoxy::Length(half));

  const Geoxy::Bezier::Cubic bend = {.Start = {.X = 0, .Y = 0}, .Out = {.X = 0, .Y = 4}, .In = {.X = 6, .Y = 4}, .End = {.X = 6, .Y = 0}};
  std::println("a cubic:            {}", Geoxy::Length(bend, Geoxy::Accuracy{}));

  const Geoxy::Place third = Geoxy::PlaceAt(half, Geoxy::Length(half) / 3);
  std::println("a third of the way: ({:.6f}, {:.6f})", third.At.X, third.At.Y);
}
```

```text [Output]
half a unit circle: 3.141592653589793
a cubic:            9.51494597792503
a third of the way: (0.500000, 0.866025)
```

## Area and centroid {#area}

`Area` of a closed curve is the area it encloses, positive when it runs counterclockwise and negative when clockwise. For an open curve it is the area between the curve and the chord that closes it. `Centroid` is the center of that area.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Polyline plate = {.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                              {.Position = {.X = 6, .Y = 0}},
                                              {.Position = {.X = 6, .Y = 2}},
                                              {.Position = {.X = 2, .Y = 2}},
                                              {.Position = {.X = 2, .Y = 5}},
                                              {.Position = {.X = 0, .Y = 5}}},
                                 .Closed   = true};
  std::println("area {}", Geoxy::Area(plate));
  const auto center = Geoxy::Centroid(plate);
  std::println("centroid ({}, {})", center->X, center->Y);

  const auto moments = Geoxy::Inertia(Geoxy::Curve{plate}, Geoxy::Accuracy{});
  std::println("second moments about the centroid: {} and {}", moments->AboutX, moments->AboutY);
}
```

```text [Output]
area 18
centroid (2.3333333333333335, 1.8333333333333333)
second moments about the centroid: 33.5 and 54
```

`Inertia` gives the second moments MASSPROP reports, and `Axes` their principal directions.

## Extents {#extents}

`Extents` is the smallest axis-aligned `Box` that holds a curve: its true extremes, not the box of its control points. Boxes combine with `Joined`, grow with `Grown` and `Padded`, and are tested with `Contains` and `Overlaps`.

```cpp [main.cpp]
#include "geoxy/geoxy.hpp"

#include <print>

auto main() -> int
{
  const Geoxy::Arc arc = {.Center = {.X = 0, .Y = 0}, .Radius = 10, .Start = Geoxy::Degrees(45), .Sweep = Geoxy::Degrees(90)};
  const Geoxy::Box box = Geoxy::Extents(arc);
  std::println("x from {:.6f} to {:.6f}", box.Min.X, box.Max.X);
  std::println("y from {:.6f} to {:.6f}", box.Min.Y, box.Max.Y);
}
```

```text [Output]
x from -7.071068 to 7.071068
y from 7.071068 to 10.000000
```

## Equal parts {#divide}

`Divide(curve, parts, accuracy)` gives the places that cut a curve into equal lengths, and `Measure(curve, length, accuracy)` the places a fixed length apart, as the DIVIDE and MEASURE commands place their points.
