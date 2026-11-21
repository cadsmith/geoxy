set(API_ROOT ${PROJECT_SOURCE_DIR}/include)

file(GLOB_RECURSE headers RELATIVE ${API_ROOT}
  CONFIGURE_DEPENDS ${API_ROOT}/geoxy/*.hpp)

foreach(header IN LISTS headers)
  string(REGEX REPLACE "\\.hpp$" ".cpp" unit ${header})
  set(path ${CMAKE_CURRENT_BINARY_DIR}/alone/${unit})
  file(CONFIGURE OUTPUT ${path}
    CONTENT "#include \"${header}\"\n#include \"${header}\"\n")
  list(APPEND units ${path})
endforeach()

add_library(geoxy_headers OBJECT ${units})

target_link_libraries(geoxy_headers
  PRIVATE
    geoxy::geoxy
    geoxy_build
)

set(API_SCRIPT
  -DAPI_ROOT=${API_ROOT}
  -DAPI_SNAPSHOT=${CMAKE_CURRENT_SOURCE_DIR}/snapshot.txt
  -DAPI_CURRENT=${CMAKE_CURRENT_BINARY_DIR}/snapshot.txt
  -P ${CMAKE_CURRENT_SOURCE_DIR}/snapshot.cmake
)

add_test(NAME geoxy_api COMMAND ${CMAKE_COMMAND} ${API_SCRIPT})

add_custom_target(geoxy_api
  COMMAND ${CMAKE_COMMAND} -DAPI_UPDATE=ON ${API_SCRIPT}
  VERBATIM
)
