// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef MATERIAL_MATERIAL_H
#define MATERIAL_MATERIAL_H

// forward declaration of material deck
namespace inp {
struct MaterialDeck;
}

/*!
 * @brief Collection of methods and database related to material
 *
 * This namespace provides collection of different materials within the
 * Peridynamics theory.
 *
 * At present we have implemented state-based model which can also simulate
 * bond-based model. We consider \b RNP regularized potential as well \b PMB
 * material.
 *
 * @sa Material
 */
namespace material {

/*! @brief A class to model a peridynamic material */
class Material {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   */
  explicit Material(inp::MaterialDeck *deck);

  /*!
   * @brief Returns true if state-based potential is active
   * @return True/false
   */
  bool isStateActive();

private:
  /**
   * @name Internal data
   */
  /**@{*/

  /*! @brief Vector of initial (reference) coordinates of nodes */
  //  std::vector<double> d_nd;

  /** @}*/
};

} // namespace material

#endif // MATERIAL_MATERIAL_H
