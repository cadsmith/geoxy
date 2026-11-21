find_program(GEOXY_Z3 z3)
find_program(GEOXY_GXX NAMES g++-16 g++-15 g++-14)
set(GEOXY_PROVED_WITH ${CMAKE_CXX_COMPILER})
if (GEOXY_GXX)
  list(APPEND GEOXY_PROVED_WITH ${GEOXY_GXX})
endif()
find_package(Python3 COMPONENTS Interpreter QUIET)

file(GLOB GEOXY_PROOF_FILES CONFIGURE_DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/*.py
  ${CMAKE_CURRENT_SOURCE_DIR}/*.v
  ${PROJECT_SOURCE_DIR}/source/exact/*.cpp
  ${PROJECT_SOURCE_DIR}/source/exact/*.hpp
)
list(APPEND GEOXY_PROOF_FILES
  ${PROJECT_SOURCE_DIR}/source/abut.cpp
  ${PROJECT_SOURCE_DIR}/source/angle.cpp
  ${PROJECT_SOURCE_DIR}/source/degrees.hpp
  ${PROJECT_SOURCE_DIR}/source/grid.cpp
  ${PROJECT_SOURCE_DIR}/source/gridwork.hpp
  ${PROJECT_SOURCE_DIR}/source/expansion.cpp
  ${PROJECT_SOURCE_DIR}/source/extents.cpp
  ${PROJECT_SOURCE_DIR}/source/flatten.cpp
  ${PROJECT_SOURCE_DIR}/source/lean.hpp
  ${PROJECT_SOURCE_DIR}/source/chord.hpp
  ${PROJECT_SOURCE_DIR}/source/prepared.cpp
  ${PROJECT_SOURCE_DIR}/source/sure.hpp
  ${PROJECT_SOURCE_DIR}/source/swift/bernstein.cpp
  ${PROJECT_SOURCE_DIR}/source/swift/casteljau.hpp
  ${PROJECT_SOURCE_DIR}/source/swift/arcade.hpp
  ${PROJECT_SOURCE_DIR}/source/swift/braid.cpp
  ${PROJECT_SOURCE_DIR}/source/swift/circuit.cpp
  ${PROJECT_SOURCE_DIR}/source/swift/circuit.hpp
  ${PROJECT_SOURCE_DIR}/source/swift/gather.cpp
  ${PROJECT_SOURCE_DIR}/source/swift/piece.cpp
  ${PROJECT_SOURCE_DIR}/source/swift/piece.hpp
  ${PROJECT_SOURCE_DIR}/source/swift/rim.cpp
  ${PROJECT_SOURCE_DIR}/source/swift/rim.hpp
  ${PROJECT_SOURCE_DIR}/source/swift/ring.cpp
  ${PROJECT_SOURCE_DIR}/source/swift/terms.hpp
  ${PROJECT_SOURCE_DIR}/source/swift/vault.cpp
  ${PROJECT_SOURCE_DIR}/source/swift/section.cpp
  ${PROJECT_SOURCE_DIR}/source/swift/strand.cpp
  ${PROJECT_SOURCE_DIR}/source/swift/twin.cpp
  ${PROJECT_SOURCE_DIR}/source/twofold.hpp
  ${PROJECT_SOURCE_DIR}/source/legendre.hpp
  ${PROJECT_SOURCE_DIR}/source/rational.cpp
  ${PROJECT_SOURCE_DIR}/source/hodograph.hpp
  ${PROJECT_SOURCE_DIR}/source/kernel/piece.hpp
  ${PROJECT_SOURCE_DIR}/include/geoxy/spline.hpp
  ${PROJECT_SOURCE_DIR}/include/geoxy/angle.hpp
  ${PROJECT_SOURCE_DIR}/include/geoxy/api.hpp
  ${PROJECT_SOURCE_DIR}/include/geoxy/point.hpp
  ${PROJECT_SOURCE_DIR}/include/geoxy/vector.hpp
  ${CMAKE_CURRENT_LIST_FILE}
)
set(GEOXY_PROOF_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(GEOXY_PROOF_STAMPS)

function(geoxy_proof name)
  set(tools ${Python3_EXECUTABLE})
  foreach(word IN LISTS ARGN)
    if (EXISTS ${word} AND NOT IS_DIRECTORY ${word})
      list(APPEND tools ${word})
    endif()
  endforeach()
  geoxy_stamp(proving proof_${name}
    INPUTS ${GEOXY_PROOF_FILES}
    TOOLS ${tools}
    COMMAND ${Python3_EXECUTABLE} -I
      ${CMAKE_CURRENT_SOURCE_DIR}/${name}.py ${ARGN}
  )
  add_custom_command(
    OUTPUT ${GEOXY_PROOF_DIR}/${name}.stamp
    COMMAND ${proving}
    COMMAND ${CMAKE_COMMAND} -E touch ${GEOXY_PROOF_DIR}/${name}.stamp
    DEPENDS ${GEOXY_PROOF_FILES} ${tools}
    COMMENT "Proving ${name}"
    VERBATIM
  )
  list(APPEND GEOXY_PROOF_STAMPS ${GEOXY_PROOF_DIR}/${name}.stamp)
  set(GEOXY_PROOF_STAMPS ${GEOXY_PROOF_STAMPS} PARENT_SCOPE)
endfunction()

if (GEOXY_Z3 AND Python3_Interpreter_FOUND)
  geoxy_proof(circle ${GEOXY_Z3})
  geoxy_proof(first ${GEOXY_Z3})
  geoxy_proof(sure ${GEOXY_Z3})
  geoxy_proof(ring ${GEOXY_Z3})
  geoxy_proof(arch ${GEOXY_Z3})
  geoxy_proof(braid ${GEOXY_Z3})
  geoxy_proof(section ${GEOXY_Z3})
  geoxy_proof(twin ${GEOXY_Z3})
  geoxy_proof(abut ${GEOXY_Z3})
  geoxy_proof(pads ${GEOXY_Z3})
  geoxy_proof(lean ${GEOXY_Z3})
  geoxy_proof(bows ${GEOXY_Z3})
  geoxy_proof(prepared ${GEOXY_Z3})
  geoxy_proof(swift ${GEOXY_Z3})
  geoxy_proof(expansion ${GEOXY_Z3})
  geoxy_proof(flatten ${GEOXY_Z3})
  geoxy_proof(strand ${GEOXY_Z3})
endif()
if (APPLE AND Python3_Interpreter_FOUND
    AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  geoxy_proof(machine ${GEOXY_PROVED_WITH})
endif()
if (Python3_Interpreter_FOUND)
  geoxy_proof(convert)
  geoxy_proof(legendre)
  geoxy_proof(rational)
endif()
if (Python3_Interpreter_FOUND AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  geoxy_proof(translation ${CMAKE_CXX_COMPILER})
  geoxy_proof(integers ${GEOXY_PROVED_WITH})
endif()
if (GEOXY_ROCQ AND Python3_Interpreter_FOUND)
  FetchContent_GetProperties(flocq SOURCE_DIR GEOXY_FLOCQ)
  geoxy_proof(rocq ${GEOXY_ROCQ} ${GEOXY_FLOCQ} ${GEOXY_PROOF_DIR})
endif()

if (GEOXY_PROOF_STAMPS)
  add_custom_target(geoxy_proofs DEPENDS ${GEOXY_PROOF_STAMPS})
  add_test(NAME geoxy_proofs
    COMMAND ${CMAKE_COMMAND}
      --build ${CMAKE_BINARY_DIR} --target geoxy_proofs)
  set_tests_properties(geoxy_proofs PROPERTIES COST 1000)
endif()
