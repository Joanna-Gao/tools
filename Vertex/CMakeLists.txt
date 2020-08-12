cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(CreateVectors)

#----------------------------------------------------------------------------
# Tell CMake where to find the ROOT installation.
#
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS}/etc/cmake)

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
find_package(ROOT REQUIRED)

#include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Add the executable, and link it to the ROOT libraries
#
add_executable(CreateVectors CreateVectors.C)
target_link_libraries(CreateVectors ${ROOT_LIBRARIES})

install(TARGETS CreateVectors DESTINATION bin)
