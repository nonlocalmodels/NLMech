////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef INP_SOLVERDECK_H
#define INP_SOLVERDECK_H

#include <iostream> // error handling
#include <vector>
#include "util/utilIO.h"

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

  /*! @brief Perturbation for the finite difference approximation in the implicit time integration */
  double d_perturbation;

  /*!
   * @brief Constructor
   */
  SolverDeck() : d_maxIters(0), d_tol(0.){};

  /*!
   * @brief Returns the string containing information about the instance of
   * the object
   *
   * @param nt Number of tabs to append before each line of string
   * @param lvl Level of information sought (higher level means more
   * information)
   * @return string String containing information about this object
   * */
  std::string printStr(int nt = 0, int lvl = 0) const {
    auto tabS = util::io::getTabS(nt);
    std::ostringstream oss;
    oss << tabS << "------- SolverDeck --------" << std::endl << std::endl;
    oss << tabS << "Solver type = " << d_solverType << std::endl;
    oss << tabS << "Max iterations = " << d_maxIters << std::endl;
    oss << tabS << "Tolerance = " << d_tol << std::endl;
    oss << tabS << "Perturbation = " << d_perturbation << std::endl;
    oss << tabS << std::endl;

    return oss.str();
  };

  /*!
   * @brief Prints the information about the instance of the object
   *
   * @param nt Number of tabs to append before each line of string
   * @param lvl Level of information sought (higher level means more
   * information)
   * */
  void print(int nt = 0, int lvl = 0) const { std::cout << printStr(nt, lvl); };
};

/** @}*/

} // namespace inp

#endif // INP_SOLVERDECK_H
