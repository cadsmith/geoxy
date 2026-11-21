# Limits

What geoxy does not do, and where its answers have bounds.

## Two dimensions {#plane}

geoxy is a library of the plane. `Space` holds only what is needed to bring a curve drawn in a tilted plane of a DXF or DWG file into the drawing's own coordinates.

## Geometry only {#scope}

There are no layers, blocks, styles, dimensions or files here. Reading and writing drawings is the work of dxfio and dwgio. Hatch patterns, text and constraints belong to the application: geoxy gives the regions, the clipping and the fills they are built from.

## Regions are lines and arcs {#regions}

`Arrangement`, `Filled`, `Union`, `Difference` and `Intersection` take lines, arcs of up to a turn and polylines. A region bounded by an ellipse or a spline is not combined yet; such a call returns `Degenerate`.

## Offsets {#offsets}

`Buffer` takes one polyline. A shape with holes is offset loop by loop and the results combined by the caller. Offsets of ellipses and splines are approximations by arcs, within the accuracy asked.

## Work budgets {#budgets}

A few operations stop after a fixed amount of work so that no input can stall a program. `Flatten` lays at most `MostChords` chords, `Divide` gives at most `MostPlaces` places, `Dashes` repeats a pattern at most `MostPatterns` times before drawing the curve whole, and a spline's degree is at most `MaxDegree`. An intersection that ran out of budget says so in `Spent`.

## Indexes {#indexes}

An index erased from a `Tree` or a `Cloud` is reused by the next insertion. A `Grid` is made for a fixed set of curves; one curve can be replaced with `Moved`, but curves are not added or removed.

## Threads {#threads}

geoxy keeps no global state that a caller sets. A `Tree`, a `Cloud` or a `Grid` must not be changed while another thread reads it. `IntersectAll` with a `Workers` is the one call that runs on several threads itself.
