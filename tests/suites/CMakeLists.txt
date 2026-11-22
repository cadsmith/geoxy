add_executable(geoxy_suites
  runner.cpp
)

target_link_libraries(geoxy_suites
  PRIVATE
    geoxy_internal
    geoxy_support
    geoxy_doctest
)

target_compile_definitions(geoxy_suites
  PRIVATE
    GEOXY_SUITES_TREE="${CMAKE_CURRENT_SOURCE_DIR}"
)

add_test(NAME suites COMMAND geoxy_suites -s)
set_tests_properties(suites PROPERTIES LABELS suites COST 100)

find_package(Python3 COMPONENTS Interpreter QUIET)
if (Python3_Interpreter_FOUND)
  set(tools ${PROJECT_SOURCE_DIR}/tools/suites)
  set(python ${Python3_EXECUTABLE} -I -B)
  geoxy_stamp(oracle suites_oracle
    INPUTS ${tools}/*.py
    TOOLS ${Python3_EXECUTABLE}
    COMMAND ${python} ${tools}/check.py
  )
  geoxy_stamp(forged suites_forged
    INPUTS ${tools}/*.py ${CMAKE_CURRENT_SOURCE_DIR}/geoxy/cases/*
    TOOLS ${Python3_EXECUTABLE}
    COMMAND ${python} ${tools}/forge.py
      --check ${CMAKE_CURRENT_SOURCE_DIR}/geoxy/cases/intersect.txt
  )
  add_test(NAME suites_oracle COMMAND ${oracle})
  set_tests_properties(suites_oracle PROPERTIES LABELS suites COST 100)
  if (APPLE)
    add_test(NAME suites_forged COMMAND ${forged})
    set_tests_properties(suites_forged PROPERTIES LABELS suites COST 100)
  endif()
endif()

add_executable(geoxy_scoreboard
  score.cpp
)

target_link_libraries(geoxy_scoreboard
  PRIVATE
    geoxy_internal
    geoxy_support
)

target_compile_definitions(geoxy_scoreboard
  PRIVATE
    GEOXY_SUITES_TREE="${CMAKE_CURRENT_SOURCE_DIR}"
)

add_test(NAME suites_scoreboard COMMAND geoxy_scoreboard --check)
set_tests_properties(suites_scoreboard PROPERTIES LABELS suites)
