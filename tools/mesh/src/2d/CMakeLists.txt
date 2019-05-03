# Copyright (c)     2017 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

AUX_SOURCE_DIRECTORY(. SOURCES) 

add_library(Mesh_TwoD ${SOURCES})
target_link_libraries(Mesh_TwoD Util IO)
