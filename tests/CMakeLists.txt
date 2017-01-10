# Need line below for EXPECT_THROW to avoid:
# "warning: non-static data member initializers only 
# available with -std=c++11 or -std=gnu++11"
add_definitions(-std=c++11)

ADD_EXECUTABLE(mdisset mdisset.cpp)
TARGET_INCLUDE_DIRECTORIES(mdisset PUBLIC "${CMAKE_SOURCE_DIR}/include/csm")

ADD_EXECUTABLE(simpletest simpletest.cpp)
TARGET_LINK_LIBRARIES(simpletest gtest_main)

# We want to compile a test from the MdisNacSensorModelTest.cpp source.
ADD_EXECUTABLE(MdisNacSensorModelTest MdisNacSensorModelTest.cpp)

# We will need these include directories to find included headers.
TARGET_INCLUDE_DIRECTORIES(MdisNacSensorModelTest PUBLIC 
  "${CMAKE_SOURCE_DIR}/include/mdis"
  "${CMAKE_SOURCE_DIR}/include")
  
# Find libcsmapi.so
FIND_LIBRARY(CSMAPI_LIBRARY csmapi "${CMAKE_SOURCE_DIR}/lib")
# The test executable needs to link against the gtest library,
# csmapi library, and MdisNacSensorModel library
TARGET_LINK_LIBRARIES(MdisNacSensorModelTest
  gtest_main
  ${CSMAPI_LIBRARY}
  MdisNacSensorModel)

ADD_EXECUTABLE(MdisPluginTest MdisPluginTest.cpp)

TARGET_INCLUDE_DIRECTORIES(MdisPluginTest PUBLIC 
  "${CMAKE_SOURCE_DIR}/include/mdis"
  "${CMAKE_SOURCE_DIR}/include")
  
TARGET_LINK_LIBRARIES(MdisPluginTest
  gtest_main
  ${CSMAPI_LIBRARY}
  MdisPlugin
  MdisNacSensorModel)
  