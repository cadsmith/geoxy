include(${PROJECT_SOURCE_DIR}/cmake/cases.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/stamp.cmake)

add_library(geoxy_internal INTERFACE)

target_link_libraries(geoxy_internal
  INTERFACE
    geoxy::geoxy
    geoxy_build
)

target_include_directories(geoxy_internal
  INTERFACE
    ${PROJECT_SOURCE_DIR}/source
)

add_library(geoxy_doctest OBJECT
  main.cpp
)

target_link_libraries(geoxy_doctest
  PUBLIC
    doctest::doctest
    geoxy_build
)

target_compile_definitions(geoxy_doctest
  PUBLIC
    DOCTEST_CONFIG_USE_STD_HEADERS
)

add_library(geoxy_support STATIC
  support/text.cpp
)

target_link_libraries(geoxy_support
  PUBLIC
    geoxy_internal
)

target_include_directories(geoxy_support
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(geoxy_tests
  fuzz/draw.cpp
  fuzz/oracle.cpp
)

target_link_libraries(geoxy_tests
  PRIVATE
    geoxy_internal
    geoxy_doctest
    geoxy_predicates
    geoxy_support
)

target_compile_definitions(geoxy_tests
  PRIVATE
    GEOXY_FUZZ_DIR="${CMAKE_CURRENT_SOURCE_DIR}/fuzz"
    GEOXY_TESTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  target_compile_options(geoxy_tests
    PRIVATE
      "SHELL:-Xclang -fno-pch-timestamp"
  )
endif()

target_precompile_headers(geoxy_tests
  PRIVATE
    <algorithm>
    <array>
    <bit>
    <cmath>
    <cstddef>
    <cstdint>
    <functional>
    <limits>
    <numbers>
    <optional>
    <span>
    <string>
    <variant>
    <vector>
    <doctest/doctest.h>
)

add_subdirectory(api)
add_subdirectory(geoxy)
add_subdirectory(proof)
add_subdirectory(suites)

if (GEOXY_ENABLE_BENCH)
  add_subdirectory(bench)
endif()

if (GEOXY_ENABLE_FUZZ)
  add_subdirectory(fuzz)
endif()

geoxy_cases(geoxy_tests ${CMAKE_CURRENT_SOURCE_DIR}/geoxy)

if (GEOXY_ENABLE_COVERAGE)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    set(LLVM_PROFDATA xcrun llvm-profdata)
    set(LLVM_COV xcrun llvm-cov)
  elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    get_filename_component(COMPILER_DIR ${CMAKE_CXX_COMPILER} DIRECTORY)
    find_program(LLVM_PROFDATA llvm-profdata HINTS ${COMPILER_DIR})
    find_program(LLVM_COV llvm-cov HINTS ${COMPILER_DIR})
  endif()

  if (LLVM_PROFDATA AND LLVM_COV)
    set(COVERAGE_SOURCES
      ${PROJECT_SOURCE_DIR}/include
      ${PROJECT_SOURCE_DIR}/source
    )
    set(COVERAGE_RAW ${CMAKE_CURRENT_BINARY_DIR}/profiles)
    add_custom_target(coverage
      COMMAND
        ${CMAKE_COMMAND} -E rm -rf ${COVERAGE_RAW}
      COMMAND
        ${CMAKE_COMMAND} -E env
        LLVM_PROFILE_FILE=${COVERAGE_RAW}/geoxy-%8m.profraw
        ${CMAKE_CTEST_COMMAND} --test-dir ${CMAKE_BINARY_DIR}
        --tests-regex ^geoxy/ --no-tests=error
        --parallel 8 --output-on-failure
      COMMAND
        ${CMAKE_COMMAND}
        "-DCOVERAGE_TOOL=${LLVM_PROFDATA}"
        -DCOVERAGE_RAW=${COVERAGE_RAW}
        -DCOVERAGE_PROFILE=geoxy_tests.profdata
        -P ${CMAKE_CURRENT_SOURCE_DIR}/coverage.cmake
      COMMAND
        ${LLVM_COV} report $<TARGET_FILE:geoxy_tests>
        -instr-profile=geoxy_tests.profdata
        -show-branch-summary
        ${COVERAGE_SOURCES}
      COMMAND
        ${LLVM_COV} show $<TARGET_FILE:geoxy_tests>
        -instr-profile=geoxy_tests.profdata
        -format=html -output-dir=coverage_html
        ${COVERAGE_SOURCES}
      COMMAND
        ${CMAKE_COMMAND}
        "-DCOVERAGE_TOOL=${LLVM_COV}"
        -DCOVERAGE_BINARY=$<TARGET_FILE:geoxy_tests>
        -DCOVERAGE_PROFILE=geoxy_tests.profdata
        "-DCOVERAGE_SOURCES=${COVERAGE_SOURCES}"
        -P ${CMAKE_CURRENT_SOURCE_DIR}/coverage.cmake
      WORKING_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}
      DEPENDS
        geoxy_tests
        geoxy_tests_cases
      VERBATIM
    )
  endif()
endif()
