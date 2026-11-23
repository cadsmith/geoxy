add_library(geoxy
  aberth.cpp
  abut.cpp
  angle.cpp
  area.cpp
  arranged.cpp
  arrangement.cpp
  at.cpp
  axis.cpp
  belt/belt.cpp
  belt/carrier.cpp
  belt/corner.cpp
  belt/cut.cpp
  belt/end.cpp
  belt/join.cpp
  belt/keep.cpp
  belt/make.cpp
  belt/near.cpp
  boughs.cpp
  box.cpp
  boxes.cpp
  chain.cpp
  chamfer.cpp
  clip.cpp
  cloud.cpp
  construct.cpp
  contain.cpp
  crossing.cpp
  curvature.cpp
  dash.cpp
  derivative.cpp
  direction.cpp
  divide.cpp
  domain.cpp
  eased.cpp
  ends.cpp
  expansion.cpp
  extents.cpp
  ferried.cpp
  fillet.cpp
  fit.cpp
  flatten.cpp
  graze.cpp
  grid.cpp
  heading.cpp
  hodograph.cpp
  inertia.cpp
  intersect.cpp
  join.cpp
  knots.cpp
  lap.cpp
  length.cpp
  loop.cpp
  matrix.cpp
  nearest.cpp
  offset.cpp
  orient.cpp
  outline.cpp
  part.cpp
  piecewise.cpp
  polish.cpp
  precise.cpp
  prepared.cpp
  rational.cpp
  region.cpp
  reverse.cpp
  scratch.cpp
  secular.cpp
  segment.cpp
  sound.cpp
  spline.cpp
  stitched.cpp
  tangency.cpp
  tangents.cpp
  tolerance.cpp
  touch.cpp
  traced.cpp
  transform.cpp
  tree.cpp
  triangulate.cpp
  vector.cpp
  version.cpp
  winding.cpp
  workers.cpp
)

add_subdirectory(corner)
add_subdirectory(crossings)
add_subdirectory(exact)
add_subdirectory(kernel)
add_subdirectory(lengths)
add_subdirectory(met)
add_subdirectory(near)
add_subdirectory(parallel)
add_subdirectory(space)
add_subdirectory(swift)
add_subdirectory(trigo)
add_subdirectory(weave)

add_library(geoxy::geoxy ALIAS geoxy)

target_include_directories(geoxy
  PRIVATE
    ${PROJECT_SOURCE_DIR}/source
)

target_compile_features(geoxy PUBLIC cxx_std_23)

find_package(Threads REQUIRED)

target_link_libraries(geoxy
  PRIVATE
    $<BUILD_INTERFACE:geoxy_build>
    Threads::Threads
)

target_compile_definitions(geoxy
  PRIVATE
    GEOXY_VERSION="${PROJECT_VERSION}"
    $<$<BOOL:${GEOXY_ENABLE_TRIGO}>:GEOXY_TRIGO>
)

set_target_properties(geoxy PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN YES
  POSITION_INDEPENDENT_CODE ON
)

get_target_property(geoxy_type geoxy TYPE)
if (geoxy_type STREQUAL "STATIC_LIBRARY")
  target_compile_definitions(geoxy PUBLIC GEOXY_STATIC)
endif()

if (GEOXY_ENABLE_TIDY)
  find_program(GEOXY_CLANG_TIDY clang-tidy REQUIRED
    HINTS /opt/homebrew/opt/llvm/bin /usr/local/opt/llvm/bin)
  set_target_properties(geoxy PROPERTIES CXX_CLANG_TIDY ${GEOXY_CLANG_TIDY})
endif()
