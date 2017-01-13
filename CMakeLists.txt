CMAKE_MINIMUM_REQUIRED (VERSION 3.0)
include(CheckCXXCompilerFlag)
# Define the version and project name
SET(CameraModel_Version_Major 0)
SET(CameraModel_Version_Minor 1)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

PROJECT(CameraModel)

# whether not tests should be built
OPTION (ENABLE_TESTS "Build the tests?" OFF)

# Add the subdirs that are being built
ADD_SUBDIRECTORY(src)

# To enable tests pass -DENABLE_TESTS=true to the cmake command
IF (ENABLE_TESTS)
	# Download and unpack googletest at configure time
	CONFIGURE_FILE(tests/CMakeGTEST.txt.in
	               googletest-download/CMakeLists.txt)
	EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
	  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
	EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} --build .
	  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )

	# Add googletest directly to our build. This adds
	# the following targets: gtest, gtest_main, gmock
	# and gmock_main
	ADD_SUBDIRECTORY(${CMAKE_BINARY_DIR}/googletest-src
	                 ${CMAKE_BINARY_DIR}/googletest-build)

	# Now simply link your own targets against gtest, gmock,
	# etc. as appropriate
	ADD_SUBDIRECTORY(tests)
ENDIF(ENABLE_TESTS)
