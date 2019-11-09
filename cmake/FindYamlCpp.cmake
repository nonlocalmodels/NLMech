# Copyright (c) 2019    Prashant K. Jha
#               2019    Patrick Diehl
#
# Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
# (See accompanying file LICENSE.txt)
find_package(PkgConfig)

find_library(YAML_CPP_LIB
        NAMES libyaml-cpp.a libyaml-cpp.so libyaml-cpp.lib
        HINTS /usr/lib64 /usr/local/lib64 /usr/lib/ /usr/local/lib "${YAML_CPP_DIR}/lib/")
find_path(YAML_CPP_INCLUDE yaml-cpp/yaml.h HINTS /usr/include /usr/local/include "${YAML_CPP_DIR}/include/")

mark_as_advanced(YAML_CPP_LIBRARY_DIR)
mark_as_advanced(YAML_CPP_LIB)
mark_as_advanced(YAML_CPP_INCLUDE)

if (NOT YAML_CPP_LIB)
    message(FATAL_ERROR "YAML CPP Library not found: Specify the YAML_CPP_DIR where yaml cpp is located")
else ()
    include_directories(${YAML_CPP_INCLUDE})
endif ()
