set(GEOXY_HEADERS
  geoxy/accuracy.hpp
  geoxy/angle.hpp
  geoxy/api.hpp
  geoxy/arc.hpp
  geoxy/area.hpp
  geoxy/arrangement.hpp
  geoxy/at.hpp
  geoxy/bezier/cubic.hpp
  geoxy/bezier/quadratic.hpp
  geoxy/box.hpp
  geoxy/chamfer.hpp
  geoxy/clip.hpp
  geoxy/cloud.hpp
  geoxy/construct.hpp
  geoxy/contain.hpp
  geoxy/curvature.hpp
  geoxy/curve.hpp
  geoxy/dash.hpp
  geoxy/derivative.hpp
  geoxy/divide.hpp
  geoxy/domain.hpp
  geoxy/ellipse.hpp
  geoxy/ends.hpp
  geoxy/error.hpp
  geoxy/extents.hpp
  geoxy/fillet.hpp
  geoxy/fit.hpp
  geoxy/flatten.hpp
  geoxy/geoxy.hpp
  geoxy/grid.hpp
  geoxy/inertia.hpp
  geoxy/intersect.hpp
  geoxy/join.hpp
  geoxy/length.hpp
  geoxy/line.hpp
  geoxy/loop.hpp
  geoxy/matrix.hpp
  geoxy/nearest.hpp
  geoxy/offset.hpp
  geoxy/orient.hpp
  geoxy/outline.hpp
  geoxy/part.hpp
  geoxy/place.hpp
  geoxy/point.hpp
  geoxy/polyline.hpp
  geoxy/prepared.hpp
  geoxy/ray.hpp
  geoxy/region.hpp
  geoxy/reverse.hpp
  geoxy/segment.hpp
  geoxy/space/frame.hpp
  geoxy/space/vector.hpp
  geoxy/spline.hpp
  geoxy/tangents.hpp
  geoxy/tolerance.hpp
  geoxy/transform.hpp
  geoxy/tree.hpp
  geoxy/triangulate.hpp
  geoxy/vector.hpp
  geoxy/version.hpp
  geoxy/workers.hpp
  geoxy/xline.hpp
)

target_sources(geoxy
  PUBLIC
    FILE_SET HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES ${GEOXY_HEADERS}
)

file(GLOB_RECURSE found RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  CONFIGURE_DEPENDS geoxy/*.hpp)

foreach(header IN LISTS found)
  if (NOT header IN_LIST GEOXY_HEADERS)
    message(FATAL_ERROR "${header} is not listed in include/CMakeLists.txt")
  endif()
endforeach()
