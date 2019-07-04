////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "testFeLib.h"
#include <hpx/hpx_main.hpp>

int main() {

  //
  // test quadrature method for triangle
  //
  {
    // test quad data for line element
    for (size_t i = 1; i < 6; i++)
      testLineElem(i);

    // test quad data for triangle element
    for (size_t i = 1; i < 6; i++)
      testTriElem(i);

    // test quad data for quadrangle element
    for (size_t i = 1; i < 6; i++)
      testQuadElem(i);

    // test additional time in computing quad points instead of storing it
    for (size_t i = 1; i < 6; i++) {
      testTriElemTime(i, 1000);

      testTriElemTime(i, 10000);

      testTriElemTime(i, 100000);

      testTriElemTime(i, 1000000);
    }
  }

  return EXIT_SUCCESS;
}
