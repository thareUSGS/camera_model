message( "External project - Spice" )

IF(APPLE)
  SET(SPICE_URL "http://naif.jpl.nasa.gov/pub/naif/toolkit/C/MacIntel_OSX_AppleC_64bit/packages/cspice.tar.Z")
ELSEIF(UNIX)
  SET(SPICE_URL "http://naif.jpl.nasa.gov/pub/naif/toolkit//C/PC_Linux_GCC_64bit/packages/cspice.tar.Z")
ENDIF()

ExternalProject_Add( Spice
  URL "${SPICE_URL}"
  DOWNLOAD_NAME cspice.tar.gz
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND
    ${CMAKE_COMMAND} -E copy
      ${CMAKE_BINARY_DIR}/Spice-prefix/src/lib/cspice.a
      ${CMAKE_BINARY_DIR}/Spice-prefix/src/lib/libcspice.a &&
    ${CMAKE_COMMAND} -E copy_directory
      ${CMAKE_BINARY_DIR}/Spice-prefix/src/Spice/
      ${INSTALL_DEPENDENCIES_DIR}/cspice
)
