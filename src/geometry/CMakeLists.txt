# Copyright (c) 2019 Prashant K. Jha
# Copyright (c) 2019 Patrick Diehl
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

AUX_SOURCE_DIRECTORY(./ SOURCES)

if(${Enable_PCL})
    add_definitions(-DENABLE_PCL)
    include_directories(${PCL_INCLUDE_DIRS})
    link_directories(${PCL_LIBRARY_DIRS})
endif()

add_hpx_library(Geometry
    SOURCES ${SOURCES})

if(${Enable_PCL})
    target_link_libraries(Geometry PUBLIC ${PCL_LIBRARIES})
endif()
