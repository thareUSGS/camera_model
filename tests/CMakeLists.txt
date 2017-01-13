ADD_EXECUTABLE(mdisset mdisset.cpp)
TARGET_INCLUDE_DIRECTORIES(mdisset PUBLIC "${CMAKE_SOURCE_DIR}/include")


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
