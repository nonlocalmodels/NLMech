# Copyright (c) 2019    Prashant K. Jha
#
# Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
# (See accompanying file LICENSE.txt)

AUX_SOURCE_DIRECTORY(./ SOURCES)
AUX_SOURCE_DIRECTORY(./decks SOURCES)

add_library(Inp ${SOURCES})
target_link_libraries(Inp Util ${YAML_CPP_LIB})
