// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef LOADING_INITIALCONDITION_H
#define LOADING_INITIALCONDITION_H

#include <string>
#include <vector>

// forward declaration of initial condition deck
namespace inp {
struct InitialConditionDeck;
}

namespace loading {

/*!
 * @brief A class to apply initial condition
 *
 * This class processes input data and provides method to apply initial
 * condition.
 */
class InitialCondition {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   */
  explicit InitialCondition(inp::InitialConditionDeck *deck);

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

#endif // LOADING_INITIALCONDITION_H
