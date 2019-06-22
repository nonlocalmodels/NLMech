// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef GEOM_FRACTURE_H
#define GEOM_FRACTURE_H

#include "util/point.h" // definition of Point3
#include <inp/decks/fractureDeck.h>
#include <stdint.h> // uint8_t type
#include <string.h> // size_t type
#include <vector>

// forward declaration of fracture deck
namespace inp {
struct EdgeCrack;
struct FractureDeck;
} // namespace inp

/*!
 * @brief Collection of methods and database related to geometry
 *
 * This namespace provides methods and data members specific to geometry. It
 * consists of class Fracture, InteriorFlags, and Neighbor.
 *
 * @sa Fracture, InteriorFlags, Neighbor
 */
namespace geometry {

/*! @brief A class for fracture state of bonds
 *
 * In this class fracture state of each bonds (i.e. whether the bond is
 * broken or not) is stored. It also comes with the access to the state
 * of bond.
 */
class Fracture {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   * @param nodes Pointer to nodal coordinates
   * @param neighbor_list Pointer to neighbor list
   */
  Fracture(inp::FractureDeck *deck, const std::vector<util::Point3> *nodes,
           const std::vector<std::vector<size_t>> *neighbor_list);

  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   */
  explicit Fracture(inp::FractureDeck *deck);

  /*!
   * @brief Sets the bond state
   *
   * @param i Nodal id
   * @param j Local id of bond in neighbor list of i
   * @param state State which is applied to the bond
   */
  void setBondState(const size_t &i, const size_t &j, const bool &state);

  /*!
   * @brief Sets the bond state
   *
   * @param i Nodal id
   * @param j Local id of bond in neighbor list of i
   * @return state True if bond is fractured otherwise false
   */
  bool getBondState(const size_t &i, const size_t &j);

  /*!
   * @brief Returns the bonds of given node i
   *
   * @param i Nodal id
   * @return List Bonds of node i
   */
  const std::vector<uint8_t> getBonds(const size_t &i);

private:
  /*!
   * @brief Sets flag of bonds of i as fractured which intersect the
   * pre-crack
   *
   * @param i Nodal id
   * @param crack Pointer to the pre-crack
   * @param nodes Pointer to nodal coordinates
   * @param neighbors Pointer to neighbors of node i
   */
  void computeFracturedBondFd(const size_t &i, inp::EdgeCrack *crack,
                              const std::vector<util::Point3> *nodes,
                              const std::vector<size_t> *neighbors);

  /*! @brief Interior flags deck */
  inp::FractureDeck *d_fractureDeck_p;

  /*! @brief Vector which stores the state of bonds
   *
   * Given node i, vector d_fracture[i] is the list of state of bonds of node
   * i.
   *
   * This is the most memory efficient data where 1 bit represents the state
   * of bond.
   */
  std::vector<std::vector<uint8_t>> d_fracture;
};

} // namespace geometry

#endif // GEOM_FRACTURE_H
