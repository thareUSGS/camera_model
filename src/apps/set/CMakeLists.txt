INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/include/mdis")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/include/csm")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/include/")


LINK_DIRECTORIES("${CMAKE_SOURCE_DIR}/lib")
#LINK_DIRECTORIES("${CMAKE_SOURCE_DIR}/src/objs")

ADD_LIBRARY(libIsdReader.so SHARED "${CMAKE_SOURCE_DIR}/src/objs/IsdReader.cpp")
TARGET_LINK_LIBRARIES(libIsdReader.so libcsmapi.so)
ADD_LIBRARY(libMdisNacSensorModel.so SHARED "${CMAKE_SOURCE_DIR}/src/objs/MdisNacSensorModel.cpp")
TARGET_LINK_LIBRARIES(libMdisNacSensorModel.so libcsmapi.so)
ADD_LIBRARY(libMdisPlugin.so SHARED "${CMAKE_SOURCE_DIR}/src/objs/MdisPlugin.cpp")
TARGET_LINK_LIBRARIES(libMdisPlugin.so libcsmapi.so)


ADD_EXECUTABLE(set set.cpp)

TARGET_LINK_LIBRARIES(set libIsdReader.so)
TARGET_LINK_LIBRARIES(set libMdisNacSensorModel.so)
TARGET_LINK_LIBRARIES(set libMdisPlugin.so)

