// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef FRACTURE_H
#define FRACTURE_H

#include <string>
#include <vector>

// forward declaration of fracture deck
namespace inp {
struct FractureDeck;
}

//! Collection of methods and database purely related to geometry
namespace geometry {

/*! @brief Methods and database associated to the mesh */
class Fracture {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   */
  Fracture(inp::FractureDeck *deck);

private:
  /**
   * \defgroup Mesh related data
   */
  /**@{*/

  /*! @brief Vector of initial (reference) coordinates of nodes */
  //  std::vector<double> d_nd;

  /** @}*/
};

} // namespace geometry

#endif // FRACTURE_H
