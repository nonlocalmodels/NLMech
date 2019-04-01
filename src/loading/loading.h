// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef LOADING_H
#define LOADING_H

#include <vector>
#include <string>

// forward declaration of loading deck
namespace inp {
struct LoadingDeck;
}

namespace loading {

/*! @brief Methods and database associated to the mesh */
class Loading {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   */
  Loading(inp::LoadingDeck *deck);

private:
  /**
   * \defgroup Mesh related data
   */
  /**@{*/

  /*! @brief Vector of initial (reference) coordinates of nodes */
//  std::vector<double> d_nd;

  /** @}*/

};

} // namespace inp

#endif // NEIGHBOR_H
