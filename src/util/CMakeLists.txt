# Copyright (c) 2019    Prashant K. Jha
#
# Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
# (See accompanying file LICENSE.txt)

AUX_SOURCE_DIRECTORY(./ SOURCES) 

add_library(Util ${SOURCES})
target_link_libraries(Util ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
