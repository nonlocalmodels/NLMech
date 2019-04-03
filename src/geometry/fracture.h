// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef GEOM_FRACTURE_H
#define GEOM_FRACTURE_H

#include <vector>

// forward declaration of fracture deck
namespace inp {
struct FractureDeck;
}

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
  explicit Fracture(inp::FractureDeck *deck);

private:
  /*! @brief Interior flags deck */
  inp::FractureDeck *d_fractureDeck_p;

  /*! @brief Vector which stores the state of bonds
   *
   * Given node i, vector d_fracture[i] is the list of state of bonds of node
   * i.
   *
   * @note 1. Since we only need single bit, as fracture state is either 0 (not
   * broken) or 1 (broken), we can further optimize the memory. Currently, we
   * store 4 bit per bond. In future we can try to reduce it to just 2 or 1
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
  std::vector<std::vector<char>> d_fracture;

  /** @}*/
};

} // namespace geometry

#endif // GEOM_FRACTURE_H
