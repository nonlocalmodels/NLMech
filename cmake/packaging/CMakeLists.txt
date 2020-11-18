# Copyright (c) 2020 Patrick Diehl
#
# SPDX-License-Identifier: BSL-1.0
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

set(CPACK_PACKAGE_VERSION ${CMAKE_PROJECT_VERSION})
set(CPACK_GENERATOR "RPM")
set(CPACK_PACKAGE_NAME ${CMAKE_PROJECT_NAME})
set(CPACK_PACKAGE_RELEASE 1)
set(CPACK_PACKAGE_CONTACT "patrickdiehl@lsu.edu")
set(CPACK_PACKAGE_VENDOR "Nonlocalmodels")
set(CPACK_RPM_PACKAGE_LICENSE "BSL-1.0")
set(CPACK_RPM_PACKAGE_URL "https://github.com/nonlocalmodels/NLMech")
set(CPACK_RPM_PACKAGE_SUMMARY
${CMAKE_PROJECT_DESCRIPTION}
)
set(CPACK_RPM_PACKAGE_DESCRIPTION
${CMAKE_PROJECT_DESCRIPTION}
)
set(CPACK_RPM_PACKAGE_REQUIRES "hpx, boost, hwloc, gperftools, vtk, yaml-cpp, blaze, blaze_iterative")
set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
set(CPACK_PACKAGE_FILE_NAME
    "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_PACKAGE_RELEASE}.${CMAKE_SYSTEM_PROCESSOR}"
)
set(CPACK_RPM_EXCLUDE_FROM_AUTO_FILELIST_ADDITION
    /etc/ld.so.conf.d
    /usr/local
    /usr/local/bin
    /usr/local/.build-id
    /usr/local/include
    /usr/local/share
    /usr/local/lib
)
set(CPACK_RPM_CHANGELOG_FILE
    "${CMAKE_HOME_DIRECTORY}/cmake/packaging/Changelog.txt"
)
include(CPack)
