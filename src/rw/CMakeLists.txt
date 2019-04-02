# Copyright (c) 2019    Prashant K. Jha
#
# Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
# (See accompanying file LICENSE.txt)

AUX_SOURCE_DIRECTORY(./ SOURCES)

add_library(RW ${SOURCES})
target_link_libraries(RW ${VTK_LIBRARIES})
