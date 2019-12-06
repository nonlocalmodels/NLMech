# Copyright (c) 2019 Prashant K. Jha
# Copyright (c) 2019 Patrick Diehl
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
cmake_minimum_required(VERSION 3.10)

project(NLMech CXX)
set (VERSION_MAJOR 0)
set (VERSION_MINOR 1)
set (VERSION_UPDATE 0)

configure_file (
  "${PROJECT_SOURCE_DIR}/Config.h.in"
  "${PROJECT_SOURCE_DIR}/src/Config.h"
  )

# Set director for all exectuables
set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin")

# Set the common library directory
get_property(LIB64 GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS)

if ("${LIB64}" STREQUAL "TRUE")
    set(LIBSUFFIX 64)
else()
    set(LIBSUFFIX "")
endif()

set(LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib${LIBSUFFIX}")

# Add a postfix if the libraries are build in the debug mode
set(CMAKE_DEBUG_POSTFIX d)

#Enable ctest
enable_testing()
add_subdirectory(test)

# Add own cmake scripts
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")

# Set build type to release if non is specified
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# Find dependencies

# HPX
find_package(HPX REQUIRED)

if (NOT (HPX_BUILD_TYPE STREQUAL CMAKE_BUILD_TYPE))
    message("HPX_BUILD_TYPE: " ${HPX_BUILD_TYPE} " CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})
    message(FATAL_ERROR "Build types used for HPX and this application have to match")
endif ()

list(APPEND CMAKE_MODULE_PATH "${HPX_PREFIX}/../cmake/")

# YAML
find_package(YamlCpp REQUIRED)

# VTK
find_package(VTK REQUIRED)

if (VTK_VERSION VERSION_LESS "6")
    message(FATAL_ERROR "ERROR: Only VTK library 6+ supported")
endif ()

# BLAZE and BlazeIterative
find_package(blaze REQUIRED)
find_package(BlazeIterative REQUIRED)

# BLAS and LAPACK
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

# Add flags to compiler
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpthread")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif ()

# Set gnu standard to 14
if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "--std=gnu++14 ${CMAKE_CXX_FLAGS}")
    endif ()
else ()
    set(CMAKE_CXX_STANDARD 17)
endif ()

# include important directories
include_directories(src/external/)
include_directories(${VTK_INCLUDE_DIRS})
include_directories(${HPX_INCLUDE_DIRS})
include_directories(src)
link_directories(${HPX_LIBRARY_DIR})

# blaze include directory (via -Dblaze_INCLUDE_DIR)
include_directories(${blaze_INCLUDE_DIR})

# Enable HPX inside blaze
add_definitions(-DBLAZE_USE_HPX_THREADS)

# subdirectories
add_subdirectory(src/util)
add_subdirectory(src/rw)
add_subdirectory(src/inp)
add_subdirectory(src/material)
add_subdirectory(src/loading)
add_subdirectory(src/geometry)
add_subdirectory(src/fe)
add_subdirectory(src/model)
add_subdirectory(src/data)

# add test directory
add_subdirectory(src/test)

# add tool directory
add_subdirectory(tools/mesh)
add_subdirectory(tools/pp)
add_subdirectory(tools/dc)

# Build executable
add_hpx_executable(NLMech SOURCES src/main.cpp 
				   DEPENDENCIES Inp Model ${LAPACK_LIBRARIES})

# Build documentation
set(Enable_Documentation FALSE CACHE BOOL "Generates target for generating the documentation")
if(${Enable_Documentation})
add_subdirectory(docs)
endif()

