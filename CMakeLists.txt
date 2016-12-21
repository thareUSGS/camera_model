CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

# Define the version and project name
SET(CameraModel_Version_Major 0)
SET(CameraModel_Version_Minor 1)
PROJECT(CameraModel)

# whether not tests should be built
OPTION (ENABLE_TESTS "Build the tests?" OFF)

# Add the subdirs that are being built
ADD_SUBDIRECTORY(src)

# TODO: Add the variable, enable GTest or CTest and execute the tests
# To enable tests pass -DENABLE_TESTS=true to the cmake command
# TODO: Look at a canonical method for getting GTest added to CMake:  https://crascit.com/2015/07/25/cmake-gtest/
IF (ENABLE_TESTS)
	ADD_SUBDIRECTORY(tests)
	MESSAGE("Building tests")
ENDIF(ENABLE_TESTS)
