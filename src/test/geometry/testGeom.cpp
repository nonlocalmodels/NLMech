// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include <hpx/hpx_main.hpp>
#include "testGeomLib.h"

int main() {

  //
  // test Fracture class and specifically the bitwise operation
  //
  testFracture();

  return EXIT_SUCCESS;
}