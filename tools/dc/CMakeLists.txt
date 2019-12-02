# Copyright (c) 2019 Prashant K. Jha
# Copyright (c) 2019 Patrick Diehl
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories(src)

AUX_SOURCE_DIRECTORY(./src SOURCES)

# build DC library
add_library(DC ${SOURCES})
target_link_libraries(DC PUBLIC FE RW Util -lpthread)

# Build dc executible
add_hpx_executable(dc SOURCES dc.cpp DEPENDENCIES DC)
