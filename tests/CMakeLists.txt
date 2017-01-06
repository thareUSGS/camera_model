ADD_EXECUTABLE(mdisset mdisset.cpp)
TARGET_INCLUDE_DIRECTORIES(mdisset PUBLIC "${CMAKE_SOURCE_DIR}/include/csm")

ADD_EXECUTABLE(simpletest simpletest.cpp)
TARGET_LINK_LIBRARIES(simpletest gtest_main)

ADD_EXECUTABLE(MdisNacSensorModelTest MdisNacSensorModelTest.cpp)
TARGET_INCLUDE_DIRECTORIES(MdisNacSensorModelTest PUBLIC 
  "${CMAKE_SOURCE_DIR}/include/mdis"
  "${CMAKE_SOURCE_DIR}/include")
FIND_LIBRARY(CSMAPI_LIBRARY csmapi "${CMAKE_SOURCE_DIR}/lib")
TARGET_LINK_LIBRARIES(MdisNacSensorModelTest
  gtest_main
  ${CSMAPI_LIBRARY}
  MdisNacSensorModel)
