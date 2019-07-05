// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef INP_SOLVERDECK_H
#define INP_SOLVERDECK_H

#include <iostream> // error handling
#include <vector>

namespace inp {

/**
 * \ingroup Input
 */
/**@{*/

/*! @brief Structure to read and store solver related input data */
struct SolverDeck {

  /*! @brief Solver type */
  std::string d_solverType;

  /*! @brief Maximum iterations */
  int d_maxIters;

  /*! @brief Tolerance */
  double d_tol;

  /*!
   * @brief Constructor
   */
  SolverDeck() : d_maxIters(0), d_tol(0.){};
};

/** @}*/

} // namespace inp

#endif // INP_SOLVERDECK_H
