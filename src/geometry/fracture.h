// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef GEOM_FRACTURE_H
#define GEOM_FRACTURE_H

#include "inp/decks/fractureDeck.h"
#include "util/feElementDefs.h"
#include "util/point.h"
#include <stdint.h>
#include <string.h>
#include <vector>

// forward declaration of fracture deck
// namespace inp {
// struct FractureDeck;
//}

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
   */
  Fracture(inp::FractureDeck *deck, const std::vector<util::Point3> *nodes,
           const std::vector<std::vector<size_t>> *neighbor_list);

  /*!
   * @brief Marks the bond as broken
   *
   * @param i Nodal id
   * @param j Local id of bond in neighborlist of i
   */
  void markBondBroken(const size_t &i, const size_t &j);

  /*!
   * @brief Sets the bond state
   *
   * @param i Nodal id
   * @param j Local id of bond in neighborlist of i
   * @param state State which is applied to the bond
   */
  void setBondState(const size_t &i, const size_t &j, bool state);

  /*!
   * @brief Sets the bond state
   *
   * @param i Nodal id
   * @param j Local id of bond in neighborlist of i
   * @return state Current state of the bond
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
   */
  void computeFracturedBondFdAlongY(const size_t &i,
                                    inp::EdgeCrack *crack,
                                    const std::vector<util::Point3> *nodes,
                                    const std::vector<size_t> *neighbors);

  /*! @brief Interior flags deck */
  inp::FractureDeck *d_fractureDeck_p;

  /*! @brief Vector which stores the state of bonds
   *
   * Given node i, vector d_fracture[i] is the list of state of bonds of node
   * i.
   *
   * @note 1. Since we only need single bit, as fracture state is either 0 (not
   * broken) or 1 (broken), we can further optimize the memory. Currently, we
   * store 8 bit per bond. In future we can try to reduce it to just 2 or 1
   * bit per bond.
   *
   * @note 2. According to <a
   * href="https://codereview.stackexchange
   * .com/questions/117880/comparing-stdvectorbool-to-stdvectorchar
   * ">Comparing vector<bool> and vector<char></a>
   * the vector<char> is faster for small size of vector. Since list of bonds
   * for each node is smaller number, we use vector<char>.
   *
   * @note 3. Both speed optimization and memory optimization is required to
   * store the fracture state efficiently.
   */
  std::vector<std::vector<uint8_t>> d_fracture;

  /** @}*/
};

} // namespace geometry

#endif // GEOM_FRACTURE_H
