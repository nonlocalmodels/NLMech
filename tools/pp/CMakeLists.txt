# Copyright (c) 2019 Prashant K. Jha
# Copyright (c) 2019 Patrick Diehl
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

AUX_SOURCE_DIRECTORY(./src SOURCES)

# build Mesh library
add_library(PP ${SOURCES})
target_link_libraries(PP Util Inp Material RW Geometry FE Data)

# Build mesh executible
add_hpx_executable(pp SOURCES pp.cpp DEPENDENCIES PP)
