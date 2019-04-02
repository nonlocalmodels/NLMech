// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef MATERIAL_H
#define MATERIAL_H

#include <string>
#include <vector>

// forward declaration of material deck
namespace inp {
struct MaterialDeck;
}

namespace material {

/*! @brief Methods and database associated to the mesh */
class Material {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   */
  Material(inp::MaterialDeck *deck);

  bool isStateActive();

private:
  /**
   * \defgroup Mesh related data
   */
  /**@{*/

  /*! @brief Vector of initial (reference) coordinates of nodes */
  //  std::vector<double> d_nd;

  /** @}*/
};

} // namespace material

#endif // NEIGHBOR_H
