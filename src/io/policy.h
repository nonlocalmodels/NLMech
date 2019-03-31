// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef POLICY_H
#define POLICY_H

#include <vector>
#include <string>

// forward declaration of policy deck
namespace io {
struct PolicyDeck;
}

namespace io {

/*! @brief Methods and database associated to the mesh */
class Policy {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   */
  Policy(io::PolicyDeck *deck);

private:
  /**
   * \defgroup Mesh related data
   */
  /**@{*/

  /*! @brief Vector of initial (reference) coordinates of nodes */
//  std::vector<double> d_nd;

  /** @}*/

};

} // namespace io

#endif // NEIGHBOR_H
