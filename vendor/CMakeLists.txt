set(FETCHCONTENT_QUIET OFF)

include(FetchContent)

if (BUILD_TESTING AND PROJECT_IS_TOP_LEVEL)
  add_subdirectory(doctest)
  add_subdirectory(predicates)
  find_program(GEOXY_ROCQ rocq)
  if (GEOXY_ROCQ)
    add_subdirectory(flocq)
  endif()
endif()

if (GEOXY_ENABLE_BENCH AND PROJECT_IS_TOP_LEVEL)
  add_subdirectory(benchmark)
endif()
