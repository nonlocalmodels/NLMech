// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "baseElem.h"
#include "../util/feElementDefs.h" // global definition of elements
#include <iostream>                // for std::cerr

fe::BaseElem::BaseElem(size_t order, size_t element_type)
    : d_quadOrder(order), d_elemType(element_type),
      d_numQuadPts(util::vtk_map_element_to_num_nodes[element_type]){};

void fe::BaseElem::init() {

  std::cerr << "Error: init() of BaseElem must be implemented in inheriting "
               "class.\n";
  exit(1);
}