INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/include/mdis")
INCLUDE_DIRECTORIES("../../../INSTALL/cspice/include")
LINK_DIRECTORIES("../../../INSTALL/cspice/lib")
LINK_DIRECTORIES("/usr/lib64")


ADD_LIBRARY(CSpiceIsd SHARED CSpiceIsd.cpp)
ADD_LIBRARY(SpiceController SHARED SpiceController.cpp)
TARGET_LINK_LIBRARIES(SpiceController libcspice.a)
TARGET_LINK_LIBRARIES(CSpiceIsd SpiceController)
ADD_EXECUTABLE(spice2isd spice2isd.cpp)
TARGET_LINK_LIBRARIES(spice2isd CSpiceIsd)
TARGET_LINK_LIBRARIES(spice2isd libgdal.so)    

