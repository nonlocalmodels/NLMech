# Copyright (c) 2019 Prashant K. Jha
# Copyright (c) 2019 Patrick Diehl
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#####################################################################
# Setup 
#####################################################################
project(NLMech
    DESCRIPTION "Framework for non-local mechanics"
    LANGUAGES CXX)

# Set the version 
set (VERSION_MAJOR 0)
set (VERSION_MINOR 1)
set (VERSION_UPDATE 0)

cmake_minimum_required(VERSION 3.10)



# Provide the configuration file to get the version
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

# Add own cmake scripts
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")

# Set build type to release if non is specified
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()


#####################################################################
# Setup the ctests
#####################################################################

enable_testing()
add_subdirectory(test)

#####################################################################
# Find dependencies
#####################################################################

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

# Gmsh
find_package(Gmsh REQUIRED)

# BLAZE and BlazeIterative
find_package(blaze REQUIRED NO_CMAKE_PACKAGE_REGISTRY)
include("${blaze_DIR}/blaze-config-version.cmake")
link_libraries(blaze::blaze)
find_package(BlazeIterative REQUIRED)

# BLAS and LAPACK
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)


set(Enable_PCL FALSE CACHE BOOL "Enables PCL for faster neighborsearch")
if(${Enable_PCL})
    find_package(PCL REQUIRED)
endif()

#####################################################################
# Set the compiler flags
#####################################################################

# Add flags to compiler
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpthread -fPIC")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread ")
endif ()

# Set gnu standard to 14
if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "--std=gnu++14 ${CMAKE_CXX_FLAGS} -fPIC")
    endif ()
else ()
    set(CMAKE_CXX_STANDARD 17)
endif ()

#####################################################################
# Set the include directories
#####################################################################

# include important directories
include_directories(src)

# Enable HPX inside blaze
add_definitions(-DBLAZE_USE_HPX_THREADS)


#####################################################################
# Set the subdirectories
#####################################################################

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
set(Enable_Expensive_Tests FALSE CACHE BOOL "Enables the computation intense tests")
add_subdirectory(src/test)

##############################################################################
# Tools
##############################################################################
set(Enable_Tools FALSE CACHE BOOL "Enables the tools to the build target")
if(${Enable_Tools})
add_subdirectory(tools)
endif()

##############################################################################
# Documentation
##############################################################################
set(Enable_Documentation FALSE CACHE BOOL "Generates target for generating the documentation")
if(${Enable_Documentation})
add_subdirectory(docs)
endif()

##############################################################################
# Build RPM
##############################################################################
set(Enable_RPM FALSE CACHE BOOL "Generates target for generating a rpm package")
if(${Enable_RPM})
add_subdirectory(cmake/packaging)
endif()


##############################################################################
# Build NLMech
##############################################################################
add_hpx_executable(NLMech SOURCES src/main.cpp 
                   DEPENDENCIES Inp Model ${LAPACK_LIBRARIES})



