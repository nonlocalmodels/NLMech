// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef LOADING_LOADING_H
#define LOADING_LOADING_H

#include <string>
#include <vector>

// forward declaration of loading deck
namespace inp {
struct LoadingDeck;
}

/*!
 * @brief Collection of methods and database related to loading
 *
 * This namespace provides methods and data members specific to application
 * of displacement and force boundary condition and also initial condition.
 *
 * @sa Loading, InitialCondition
 */
namespace loading {

/*!
 * @brief A class to apply displacement and force boundary condition
 *
 * In this class we process input data and apply complex boundary condition.
 * The boundary conditions can be specified in multiple sets, and in each set
 * one can specify the region where the boundary condition is to be applied,
 * and the type of function to be used.
 *
 * This class also provides method to set the fixity of nodes as fixed if the
 * displacement is specified on the dof of the node.
 */
class Loading {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   */
  Loading(inp::LoadingDeck *deck);

private:
  /**
   * @name Internal data
   */
  /**@{*/

  /*! @brief Vector of initial (reference) coordinates of nodes */
  //  std::vector<double> d_nd;

  /** @}*/
};

} // namespace loading

#endif // LOADING_LOADING_H
