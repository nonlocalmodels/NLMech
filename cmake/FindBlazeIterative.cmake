# Copyright (c) 2018 Patrick Diehl
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

find_package(PkgConfig QUIET)

find_path(BLAZEITERATIVE_INCLUDE BlazeIterative.hpp
	HINTS "${BLAZEITERATIVE_DIR}/include/" "${BLAZEITERATIVE_DIR}/include/BlazeIterative/" /usr/include/BlazeIterative/ /usr/local/include/BlazeIterative/ )

mark_as_advanced(BLAZEITERATIVE_DIR)
mark_as_advanced(BLAZEITERATIVE_INCLUDE)

if(NOT BLAZEITERATIVE_INCLUDE)
  message(FATAL_ERROR "Blazeiterative could not be found. Please specify BLAZEITERATIVE_DIR to assist locating it")
else()
  set(BLAZEITERATIVE_FOUND ON)
  include_directories(${BLAZEITERATIVE_INCLUDE})
#  message(Blaze iterative: ${BLAZEITERATIVE_INCLUDE})
endif()
