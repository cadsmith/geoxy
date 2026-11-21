add_executable(geoxy_bench EXCLUDE_FROM_ALL
  exact/batch.cpp
  exact/circle.cpp
  exact/orient.cpp
  fill.cpp
  grid.cpp
  intersect.cpp
  measures.cpp
  offset.cpp
  region.cpp
)

target_link_libraries(geoxy_bench
  PRIVATE
    geoxy_internal
    geoxy_predicates
    benchmark::benchmark_main
)

