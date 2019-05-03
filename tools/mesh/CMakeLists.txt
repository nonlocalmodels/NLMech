# Copyright (c)    
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

AUX_SOURCE_DIRECTORY(./1d SOURCES)
AUX_SOURCE_DIRECTORY(./2d SOURCES)

# Build mesh
add_hpx_executable(mesh SOURCES Main.cpp DEPENDENCIES Util Inp)
