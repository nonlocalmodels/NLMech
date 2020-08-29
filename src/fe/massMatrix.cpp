////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "massMatrix.h"

#include <iostream>  // for std::cerr

#include "inp/decks/massMatrixDeck.h"
#include "inp/policy.h"  // to check for constraints
#include "util/utilIO.h"

fe::MassMatrix::MassMatrix(inp::MassMatrixDeck *deck) {
  d_massMatrixDeck_p = deck;

  if (d_massMatrixDeck_p->d_MApproxType != "exact" and
      d_massMatrixDeck_p->d_MApproxType != "lumped") {
    std::cerr << "Error: Check M_Matrix_Approx tag in input data. We support "
                 " exact and lumped as two options for computation of mass "
                 "matrix.\n";
    exit(1);
  }

  // check memory control flag and it is 2 or higher enforce approximation of
  // mass matrix by lumping method
  if (inp::Policy::getInstance()->getMemoryControlFlag() >= 2)
    if (d_massMatrixDeck_p->d_MApproxType != "lumped")
      d_massMatrixDeck_p->d_MApproxType = "lumped";
}

std::string fe::MassMatrix::printStr(int nt, int lvl) const {
  auto tabS = util::io::getTabS(nt);
  std::ostringstream oss;
  oss << tabS << "------- MassMatrix --------" << std::endl << std::endl;
  oss << tabS << "Mass matrix deck address = " << d_massMatrixDeck_p << std::endl;
  oss << tabS << "Size of diagonal matrix = " << d_invMDiag.size() << std::endl;
  oss << tabS << std::endl;

  return oss.str();
}
