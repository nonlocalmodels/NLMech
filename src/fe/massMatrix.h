// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef FE_MASSMATRIX_H
#define FE_MASSMATRIX_H

#include "util/matrixBlaze.h" // definition of SymMatrixFij

// forward declaration
namespace inp {
struct MassMatrixDeck;
}

namespace fe {

/*! @brief A class for mass matrix
 *
 * In this class we compute and store inverse of a mass matrix. We can either
 * compute exact mass matrix and its inverse (using Blaze library) or we can
 * use row-sum approximation (lumping) to approximate the matrix as diagonal
 * matrix.
 *
 * User can specify the order of quadrature approximation in computation of
 * elements of mass matrix.
 *
 * @note This class is needed only when the discretization is of
 * weak_finite_element type. It depends on the Mesh::d_gMap and Mesh::d_gInvMap
 * of Mesh class.
 *
 */

class MassMatrix {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   */
  explicit MassMatrix(inp::MassMatrixDeck *deck);

private:
  /*! @brief Mass matrix deck
   *
   * Deck contains information such as order of quadrature approximation and
   * type of approximation to be used for mass matrix.
   *
   * @sa inp::MassMatrixDeck
   * */
  inp::MassMatrixDeck *d_massMatrixDeck_p;

  /**
   * @name Mass matrix
   */
  /**@{*/

  /*! @brief Inverse of a diagonal mass matrix stored in a vector */
  std::vector<double> d_invMDiag;

  /*! @brief Inverse of exact mass matrix stored
   *
   * We use Blaze matrix data type.
   *
   * Currently we use float for each element of matrix to reduce the memory
   * load.
   *
   * @sa util::SymMatrixFij
   */
  util::SymMatrixFij d_invM;

  /** @}*/
};

} // namespace fe

#endif // FE_MASSMATRIX_H
