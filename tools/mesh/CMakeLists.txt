# Copyright (c) 2019 Prashant K. Jha
# Copyright (c) 2019 Patrick Diehl
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

AUX_SOURCE_DIRECTORY(./src/1d SOURCES)
AUX_SOURCE_DIRECTORY(./src/2d SOURCES)

# build Mesh library
add_hpx_library(Mesh 
    DEPENDENCIES ${YAML_CPP_LIB} RW Util Data 
    SOURCES ${SOURCES})
#target_link_libraries(Mesh Util RW ${YAML_CPP_LIB} ${VTK_LIBRARIES} Data FE)

# Build mesh executible
add_hpx_executable(mesh SOURCES mesh.cpp DEPENDENCIES Mesh FE)
