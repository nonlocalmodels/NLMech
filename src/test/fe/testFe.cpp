// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include <hpx/hpx_main.hpp>
#include "testFeLib.h"

int main() {

  //
  // test quadrature method for triangle
  //
  {
    // test quad data for triangle element
    for (size_t i=1; i<6; i++)
      testTriElem(i);

    // test quad data for quadrangle element
    for (size_t i=1; i<6; i++)
      testQuadElem(i);

  }

  return EXIT_SUCCESS;
}