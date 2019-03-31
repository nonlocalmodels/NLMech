// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef PROBLEMDECK_H
#define PROBLEMDECK_H

#include <vector>
#include <iostream>       // error handling

namespace io {

/*! @brief Structure to read and store policy data */
struct ProblemDeck {

  /*! @brief Tag of spatial discretization. Possible values are: "",
   * "finite_difference", "weak_fe", "nodal_fe", "truss_fe"
   */
  std::string d_spatialDiscretization;
};
} // namespace io
#endif // PROBLEMDECK_H
