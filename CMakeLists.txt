CMAKE_MINIMUM_REQUIRED (VERSION 3.0)
include_directories(include)
include(CheckCXXCompilerFlag)
# Define the version and project name
SET(CameraModel_Version_Major 0)
SET(CameraModel_Version_Minor 1)

PROJECT(CameraModel)

# Set additional module search paths
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

#Setup to use CPP11
ADD_COMPILE_OPTIONS(-std=c++11)
# Setup for external dependencies (either use system or install)
if( NOT INSTALL_DEPENDENCIES_DIR )
    SET( INSTALL_DEPENDENCIES_DIR ${CMAKE_BINARY_DIR}/INSTALL CACHE STRING "Install directory for dependencies")
endif()

# Add external project capability
INCLUDE(ExternalProject)

# Provide the ability to download dependencies and build (default: download)
OPTION (USE_SYSTEM_EIGEN "Use system libraries for Eigen" OFF)
if (${USE_SYSTEM_EIGEN} MATCHES "OFF")
	include("${CMAKE_MODULE_PATH}/External-Eigen.cmake")
else()
	FIND_PACKAGE(Eigen3 REQUIRED 3.4)
endif()

# Find Spice
OPTION (USE_SYSTEM_SPICE "Use system libraries for the NAIFSpice Toolkit" OFF)
if (${USE_SYSTEM_SPICE} MATCHES "OFF")
  INCLUDE("${CMAKE_MODULE_PATH}/External-Spice.cmake")
endif()
FIND_PACKAGE(Spice)

# Find GDAL
#FIND_PACKAGE(GDAL REQUIRED)

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
