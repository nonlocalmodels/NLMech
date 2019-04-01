# Copyright (c) 2019    Prashant K. Jha
#
# Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
# (See accompanying file LICENSE.txt)

cmake_minimum_required(VERSION 3.0)

project(NLMech CXX)
enable_testing()

# Add own cmake scripts
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")

#cmake_policy(SET CMP0023 OLD)

# Find dependencies

# HPX
find_package(HPX REQUIRED)

if(NOT (HPX_BUILD_TYPE STREQUAL CMAKE_BUILD_TYPE))
    message("HPX_BUILD_TYPE: " ${HPX_BUILD_TYPE} " CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})
    message(FATAL_ERROR "Build types used for HPX and this library have to match")
endif()

list(APPEND CMAKE_MODULE_PATH "${HPX_PREFIX}/../cmake/")

if(${HPX_FOUND})
    message("${HPX_INCLUDE_DIRS}")
endif()

# YAML
find_package(YamlCpp REQUIRED)

# VTK
find_package(VTK REQUIRED)

if (VTK_VERSION VERSION_LESS "6")
    message(FATAL_ERROR "ERROR: Only VTK library 6+ supported")
endif ()

# Blas and Lapack
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

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
    set(CMAKE_CXX_STANDARD 14)
endif ()

# include important directories
include_directories(src/external/blaze/include/)
include_directories(${VTK_INCLUDE_DIRS})
include_directories(${HPX_INCLUDE_DIRS})
include_directories(src)
link_directories(${HPX_LIBRARY_DIR})

# Add subdirectories
add_subdirectory(src/util)
add_subdirectory(src/inp)
add_subdirectory(src/material)
add_subdirectory(src/loading)
add_subdirectory(src/geometry)
add_subdirectory(src/model)
add_subdirectory(docs)

# Build executable
add_hpx_executable(NLMech SOURCES src/main.cpp DEPENDENCIES Util IO
        Material Loading Geometry Model ${YAML_CPP_LIB})