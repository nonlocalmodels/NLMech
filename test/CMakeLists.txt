# Copyright (c) 2019 Prashant K. Jha
# Copyright (c) 2019 Patrick Diehl
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

##############################################################################
# Quasi-static tests
##############################################################################

# 1D
add_test(NAME quasistatic.1D.elastic.mesh 
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/mesh -i input_mesh.yaml -d 1
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/examples/qsModel/1D
)

add_test(NAME quasistatic.1D.elastic
COMMAND ${EXECUTABLE_OUTPUT_PATH}/NLMech -i input.yaml --hpx:threads=1
WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/examples/qsModel/1D
)

add_test(NAME quasistatic.1D.elastic.compare
COMMAND ${EXECUTABLE_OUTPUT_PATH}/dc -i input_compare.yaml -k fd_simple
WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/examples/qsModel/1D
)

# Depencies
set_tests_properties(quasistatic.1D.elastic PROPERTIES
  FIXTURES_SETUP quasistatic.1D.elastic.mesh)

set_tests_properties(quasistatic.1D.elastic.compare PROPERTIES
  FIXTURES_SETUP quasistatic.1D.elastic)

# The implicit test tooks too long for circle-ci and can only be run
# on desktops or clusters.
if(Enable_Expensive_Tests)

# 2D
add_test(NAME quasistatic.2D.elastic.mesh 
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/mesh -i input_mesh.yaml -d 2
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/examples/qsModel/2D
)

add_test(NAME quasistatic.2D.elastic
COMMAND ${EXECUTABLE_OUTPUT_PATH}/NLMech -i input.yaml --hpx:threads=5
WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/examples/qsModel/2D
)

add_test(NAME quasistatic.2D.elastic.compare
COMMAND ${EXECUTABLE_OUTPUT_PATH}/dc -i input_compare.yaml -k fd_simple
WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/examples/qsModel/2D
)

endif()

##############################################################################
# Tools
##############################################################################

# Compare tool
add_test(NAME tools.dc.fd_simple 
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/dc -i input.yaml -k fd_simple
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/examples/tools/dc/fd_simple
)

# Generate 2D mesh
add_test(NAME tools.mesh.2d
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/mesh -i input_mesh.yaml -d 2
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/examples/tools/mesh/uniform_triangle
)

# Read 2D traingle mesh generate with gmsh
add_test(NAME mesh.2d.read.triangle.gmsh
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/NLMech -i input_2d_gmsh_triangle.yaml --hpx:threads=1
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/examples/io
)

# Read 2D traingle mesh generate with gmsh
add_test(NAME mesh.2d.read.quad.gmsh
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/NLMech -i input_2d_gmsh_quad.yaml --hpx:threads=1
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/examples/io
)


##############################################################################
# IO
##############################################################################

# Test the legacy vtk writer 
add_test(NAME quasistatic.1D.elastic_legacy_vtk
COMMAND ${EXECUTABLE_OUTPUT_PATH}/NLMech -i input_legacy_vtk.yaml --hpx:threads=1
WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/examples/qsModel/1D
)

# Depencies
set_tests_properties(quasistatic.1D.elastic_legacy_vtk PROPERTIES
  FIXTURES_SETUP quasistatic.1D.elastic.mesh)