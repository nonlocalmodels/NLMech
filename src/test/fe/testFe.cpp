// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "testFeLib.h"
#include <hpx/hpx_main.hpp>

int main() {

  //
  // test quadrature method for triangle
  //
  {
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