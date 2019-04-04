# Copyright (c) 2017
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# Build test
add_hpx_executable(TestFE SOURCES testFe.cpp testFeLib.cpp DEPENDENCIES Util FE)
