# Copyright (c) 2020    Prashant K. Jha
#               2020    Patrick Diehl
#
# Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
# (See accompanying file LICENSE.txt)
find_package(PkgConfig)

find_library(GMSH_LIB
        NAMES libgmsh.a libgmsh.so
        HINTS /usr/lib64 /usr/local/lib64 /usr/lib/ /usr/local/lib "${GMSH_DIR}"
	    PATH_SUFFIXES lib lib64)
find_path(GMSH_INCLUDE gmsh.h HINTS /usr/include /usr/local/include "${GMSH_DIR}/include/")

mark_as_advanced(GMSH_DIR)
mark_as_advanced(GMSH_LIB)
mark_as_advanced(GMSH_INCLUDE)

if (NOT GMSH_LIB)
    message(FATAL_ERROR "Gmsh Library not found: Specify the GMSH_DIR where Gmsh is located")
else ()
    include_directories(${GMSH_INCLUDE})
endif ()
