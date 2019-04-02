// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef INITIALCONDITION_H
#define INITIALCONDITION_H

#include <string>
#include <vector>

// forward declaration of initial condition deck
namespace inp {
struct InitialConditionDeck;
}

namespace loading {

/*! @brief Methods and database associated to the mesh */
class InitialCondition {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   */
  InitialCondition(inp::InitialConditionDeck *deck);

private:
  /**
   * \defgroup Mesh related data
   */
  /**@{*/

  /*! @brief Vector of initial (reference) coordinates of nodes */
  //  std::vector<double> d_nd;

  /** @}*/
};

} // namespace loading

#endif // NEIGHBOR_H
