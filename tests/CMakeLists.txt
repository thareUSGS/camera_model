ADD_EXECUTABLE(mdisset mdisset.cpp)
TARGET_INCLUDE_DIRECTORIES(mdisset PUBLIC "${CMAKE_SOURCE_DIR}/include/csm")

ADD_EXECUTABLE(simpletest simpletest.cpp)
TARGET_LINK_LIBRARIES(simpletest gtest_main)
