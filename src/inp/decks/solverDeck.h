// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef SOLVERDECK_H
#define SOLVERDECK_H

#include <vector>
#include <iostream>       // error handling

namespace inp {

/*! @brief Structure to read and store solver related input data */
struct SolverDeck {

  /**
   * \defgroup Data members
   */
  /**@{*/

  /*! @brief Solver type */
  std::string d_solverType;

  /*! @brief Maximum iterations */
  int d_maxIters;

  /*! @brief Tolerance */
  double d_tol;

  /** @}*/

  /*!
   * @brief Constructor
   */
  SolverDeck() : d_maxIters(0), d_tol(0.){};
};

} // namespace inp
#endif // SOLVERDECK_H
