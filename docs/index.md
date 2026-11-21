---
layout: home

hero:
  name: geoxy
  text: Exact and fast 2D geometry for CAD
  tagline: The curves a drawing is made of, and what drawing, editing and snapping need of them. C++23, no dependencies, no exceptions, and answers right to the last digit a double holds.
  actions:
    - theme: brand
      text: Get started
      link: /guide/start/
    - theme: alt
      text: See the tutorials
      link: /tutorials/intersecting/

features:
  - title: Where curves meet
    details: Every crossing, touch and shared stretch of any two curves, from lines to rational splines. Between lines and arcs each point is the nearest double to the true one.
    link: /tutorials/intersecting/
    linkText: Intersect
  - title: Offsets that hold up
    details: One side of a curve as OFFSET gives it, or a shape grown and shrunk as a region, with round, sharp or chamfered corners.
    link: /tutorials/offsetting/
    linkText: Offset and Buffer
  - title: The editing commands
    details: Trim, extend, break, join, fillet and chamfer, each a short call that returns a new curve and leaves the old one alone.
    link: /tutorials/editing/
    linkText: Edit
  - title: Regions
    details: The area a click lies in, with its islands, and unions, differences and intersections that keep arcs as arcs.
    link: /tutorials/regions/
    linkText: Regions
  - title: Indexes for large drawings
    details: A tree of curves and a cloud of points for picking and snapping while editing, and a grid for every intersection at once.
    link: /tutorials/indexing/
    linkText: Tree, Cloud and Grid
  - title: Nothing hidden
    details: No global state, no exceptions, nothing beyond the standard library. A call that can fail says why in its return type.
    link: /guide/errors/
    linkText: Errors
---

## Measured, not claimed

Every figure here comes from a run against an independent reference: exact arithmetic, a 120 digit offset program, or asking every curve one by one.

<div class="numbers">
  <div class="number"><div class="figure">0</div><div class="what">wrong offsets</div><div class="how">of 37,888 shapes judged against the reference; 88% identical to the last bit.</div></div>
  <div class="number"><div class="figure">7.0 ms</div><div class="what">every crossing of 100,000 curves</div><div class="how">167,701 points, each one the nearest double to the truth.</div></div>
  <div class="number"><div class="figure">2 µs</div><div class="what">nearest curve to the cursor</div><div class="how">among 100,000, against about 450 µs asking each in turn.</div></div>
  <div class="number"><div class="figure">100%</div><div class="what">of lines and branches tested</div><div class="how">no change is kept unless the tests still run all of the library.</div></div>
</div>

## A few lines each

Curves are plain values, operations are free functions, and a result is checked before it is used.

::: code-group

```cpp [Intersect]
const Geoxy::Curve     line      = Geoxy::Line{.Start = {.X = -10, .Y = 3}, .End = {.X = 10, .Y = 3}};
const Geoxy::Curve     circle    = Geoxy::Circle({.X = 0, .Y = 0}, 5);
const Geoxy::Tolerance tolerance = {.Length = 1e-9, .Radians = 1e-9};

for (const Geoxy::Hit &hit : Geoxy::Intersect(line, circle, tolerance).Hits)
{
  std::println("({}, {})", hit.At.X, hit.At.Y); // (-4, 3) and (4, 3)
}
```

```cpp [Offset]
const Geoxy::Polyline plate = {.Vertices = {{.Position = {.X = 0, .Y = 0}},
                                            {.Position = {.X = 10, .Y = 0}},
                                            {.Position = {.X = 10, .Y = 6}},
                                            {.Position = {.X = 0, .Y = 6}}},
                               .Closed   = true};

// The plate grown by 1, its corners rounded: one loop of area 95.14.
const auto grown = Geoxy::Buffer(plate, {.Distance = 1, .Sides = Geoxy::Side::Outward}, tolerance);
```

```cpp [Snap]
Geoxy::Tree tree{curves};

// The three curves nearest the cursor, within the pick aperture.
for (const Geoxy::Tree::Close &close : tree.Nearest(cursor, 3, aperture))
{
  std::println("curve {} at ({}, {})", close.Index, close.Nearest.At.X, close.Nearest.At.Y);
}
```

```cpp [Fillet]
const auto corner = Geoxy::Fillet(across, upward, {.Radius = 2, .Picks = {pickOne, pickTwo}});
if (corner)
{
  draw(corner->First, *corner->Round, corner->Second);
}
```

```cmake [CMake]
add_subdirectory(geoxy)
target_link_libraries(app PRIVATE geoxy::geoxy)
```

:::

## Part of a family

<div class="family">
  <div class="here"><b>geoxy</b><span>The geometry: curves and what is worked out from them.</span></div>
  <a href="https://github.com/cadsmith/dxfio"><b>dxfio</b><span>Reads and writes DXF drawings.</span></a>
  <div><b>dwgio</b><span>Reads and writes DWG drawings.</span></div>
</div>
