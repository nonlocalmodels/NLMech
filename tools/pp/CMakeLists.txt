# Copyright (c) 2019    Prashant K. Jha
#
# Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
# (See accompanying file LICENSE.txt)

AUX_SOURCE_DIRECTORY(./src SOURCES)

# build Mesh library
add_library(PP ${SOURCES})
target_link_libraries(PP Util Inp Material RW Geometry FE)

# Build mesh executible
add_hpx_executable(pp SOURCES pp.cpp DEPENDENCIES PP)
