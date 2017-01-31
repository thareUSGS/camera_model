INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/include/mdis")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/include/csm")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/include/")

# Jenkins puts its builds under _build dir, probably a better way to do this
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/_build/INSTALL/cspice/include")
LINK_DIRECTORIES("${CMAKE_SOURCE_DIR}/_build/INSTALL/cspice/lib")

ADD_EXECUTABLE(set set.cpp)
TARGET_LINK_LIBRARIES(set dl ${CSMAPI_LIBRARY} IsdReader MdisNacSensorModel MdisPlugin gdal)
