include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(GEOXY_INSTALL_CMAKEDIR ${CMAKE_INSTALL_LIBDIR}/cmake/geoxy)

install(TARGETS geoxy EXPORT geoxyTargets FILE_SET HEADERS)

install(EXPORT geoxyTargets
  NAMESPACE geoxy::
  DESTINATION ${GEOXY_INSTALL_CMAKEDIR}
)

configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/cmake/config.cmake.in
  ${PROJECT_BINARY_DIR}/geoxyConfig.cmake
  INSTALL_DESTINATION ${GEOXY_INSTALL_CMAKEDIR}
)

write_basic_package_version_file(
  ${PROJECT_BINARY_DIR}/geoxyConfigVersion.cmake
  COMPATIBILITY SameMinorVersion
)

install(
  FILES
    ${PROJECT_BINARY_DIR}/geoxyConfig.cmake
    ${PROJECT_BINARY_DIR}/geoxyConfigVersion.cmake
  DESTINATION ${GEOXY_INSTALL_CMAKEDIR}
)
