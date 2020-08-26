# Copyright (c) 2019 Prashant K. Jha
# Copyright (c) 2019 Patrick Diehl
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

set(ProjectName "UScale")
set(ProjectExeName "uscale")

# sources
AUX_SOURCE_DIRECTORY(./src SOURCES)

# build library
add_library(${ProjectName} ${SOURCES})
target_include_directories(${ProjectName} PUBLIC src)
target_link_libraries(${ProjectName} PUBLIC FE RW Util -lpthread ${YAML_CPP_LIB} Data)

# Build executible
add_hpx_executable(${ProjectExeName} SOURCES main.cpp DEPENDENCIES ${ProjectName})
