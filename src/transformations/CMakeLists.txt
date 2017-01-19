INCLUDE_DIRECTORIES("${EIGEN3_INCLUDE_DIR}")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/include/transformations")

ADD_LIBRARY(Transformations SHARED transformations.cpp)
