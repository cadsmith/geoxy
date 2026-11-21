if (CMAKE_SCRIPT_MODE_FILE)
  set(ENV{LLVM_PROFILE_FILE} ${CASES_OUTPUT}.profraw)
  execute_process(
    COMMAND ${CASES_PROGRAM} --list-test-cases --reporters=xml
    OUTPUT_VARIABLE listed
    COMMAND_ERROR_IS_FATAL ANY
  )
  set(text "")
  set(tests "\n")
  set(rest "${listed}")
  while (TRUE)
    string(FIND "${rest}" "<TestCase " start)
    if (start EQUAL -1)
      break()
    endif()
    string(SUBSTRING "${rest}" ${start} -1 rest)
    string(FIND "${rest}" ">" end)
    string(SUBSTRING "${rest}" 0 ${end} element)
    string(SUBSTRING "${rest}" ${end} -1 rest)
    string(REGEX MATCH " filename=\"([^\"]*)\"" found "${element}")
    file(TO_CMAKE_PATH "${CMAKE_MATCH_1}" source)
    file(RELATIVE_PATH file ${CASES_ROOT} "${source}")
    if (file MATCHES "^\\.\\.")
      continue()
    endif()
    string(REGEX REPLACE "\\.cpp$" "" file "${file}")
    string(REGEX MATCH "^<TestCase name=\"([^\"]*)\"" found "${element}")
    set(name "${CMAKE_MATCH_1}")
    string(REPLACE "&quot;" "\"" name "${name}")
    string(REPLACE "&apos;" "'" name "${name}")
    string(REPLACE "&lt;" "<" name "${name}")
    string(REPLACE "&gt;" ">" name "${name}")
    string(REPLACE "&amp;" "&" name "${name}")
    set(test "geoxy/${file}/${name}")
    string(FIND "${tests}" "\n${test}\n" again)
    if (NOT again EQUAL -1)
      message(FATAL_ERROR "Two test cases are named ${test}")
    endif()
    string(APPEND tests "${test}\n")
    string(REPLACE "\\" "\\\\" filter "${name}")
    string(REPLACE "," "\\," filter "${filter}")
    string(REGEX REPLACE "[^ -~]" "?" filter "${filter}")
    string(REPLACE "/" "?" path "${file}")
    string(APPEND text
      "add_test([==[${test}]==] [==[${CASES_PROGRAM}]==]"
      " [==[-tc=${filter}]==] [==[-sf=*tests?geoxy?${path}.cpp]==])\n"
      "set_tests_properties([==[${test}]==] PROPERTIES"
      " FAIL_REGULAR_EXPRESSION [==[test cases: +0 \\|]==]"
      " WORKING_DIRECTORY [==[${CASES_DIRECTORY}]==])\n"
    )
  endwhile()
  if (text STREQUAL "")
    message(FATAL_ERROR "${CASES_PROGRAM} lists no test case")
  endif()
  file(WRITE ${CASES_OUTPUT} "${text}")
  file(REMOVE ${CASES_OUTPUT}.profraw)
  return()
endif()

set(GEOXY_CASES_SCRIPT ${CMAKE_CURRENT_LIST_FILE})

function(geoxy_cases target root)
  set(listing ${CMAKE_CURRENT_BINARY_DIR}/${target}_cases.cmake)
  add_custom_command(
    OUTPUT ${listing}
    COMMAND ${CMAKE_COMMAND}
      -DCASES_PROGRAM=$<TARGET_FILE:${target}>
      -DCASES_ROOT=${root}
      -DCASES_DIRECTORY=${CMAKE_CURRENT_BINARY_DIR}
      -DCASES_OUTPUT=${listing}
      -P ${GEOXY_CASES_SCRIPT}
    DEPENDS ${target} ${GEOXY_CASES_SCRIPT}
    COMMENT "Listing the test cases of ${target}"
    VERBATIM
  )
  add_custom_target(${target}_cases ALL DEPENDS ${listing})
  set(include ${CMAKE_CURRENT_BINARY_DIR}/${target}_include.cmake)
  string(CONCAT text
    "if (EXISTS [==[${listing}]==])\n"
    "  include([==[${listing}]==])\n"
    "else()\n"
    "  add_test([==[${target}]==] [==[$<TARGET_FILE:${target}>]==])\n"
    "endif()\n"
  )
  file(GENERATE OUTPUT ${include} CONTENT "${text}")
  set_property(DIRECTORY APPEND PROPERTY TEST_INCLUDE_FILES ${include})
endfunction()
