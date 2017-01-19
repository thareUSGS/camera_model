
#ADD_EXECUTABLE(mdisset mdisset.cpp)
#TARGET_INCLUDE_DIRECTORIES(mdisset PUBLIC "${CMAKE_SOURCE_DIR}/include")


# Need line below for EXPECT_THROW to avoid:
# "warning: non-static data member initializers only
# available with -std=c++11 or -std=gnu++11"

# Find libcsmapi.so
FIND_LIBRARY(CSMAPI_LIBRARY csmapi "${CMAKE_SOURCE_DIR}/lib")

# Find the unit tests
FILE(GLOB SRC_FILES ${CMAKE_SOURCE_DIR}/tests/ *Test.cpp)

ADD_EXECUTABLE(runTests runTests.cpp ${SRC_FILES} ${HEADER_FILES})
TARGET_INCLUDE_DIRECTORIES(runTests PUBLIC
  ${CMAKE_SOURCE_DIR}/include/mdis
  ${CMAKE_SOURCE_DIR}/include/json
  ${CMAKE_SOURCE_DIR}/include/
  ${CMAKE_SOURCE_DIR}/include/transformations
  ${EIGEN3_INCLUDE_DIR})

TARGET_LINK_LIBRARIES(runTests gtest_main
                      MdisPlugin
                      MdisNacSensorModel
                      IsdReader
                      Transformations
                      ${CSMAPI_LIBRARY})
