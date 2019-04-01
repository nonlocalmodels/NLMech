# Copyright (c) 2017
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

find_package(Doxygen)

add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile} conf.doxy
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/docs
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
        
add_custom_target(doc_develop
        COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile} conf_develop.doxy
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/docs
        COMMENT "Generating API documentation for development purpose with Doxygen"
        VERBATIM)
