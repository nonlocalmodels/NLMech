////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

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
 * In this class we compute and store the inverse of a mass matrix. We can
 * either compute the exact mass matrix and its inverse (using Blaze library) or
 * we can use row-sum approximation (lumping) to approximate the matrix as diagonal
 * matrix.
 *
 * Order of quadrature approximation in computation of elements of mass
 * matrix can be specified in the input file, see inp::MassMatrixDeck.
 *
 * @note This class is needed only when the discretization is of
 * \b weak_finite_element type. This class requires data fe::Mesh::d_gMap and
 * fe::Mesh::d_gInvMap.
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
   */
  inp::MassMatrixDeck *d_massMatrixDeck_p;

  /*! @brief Inverse of a diagonal mass matrix stored as a vector */
  std::vector<double> d_invMDiag;

  /*! @brief Inverse of the exact mass matrix
   *
   * We use Blaze matrix library with float data type.
   */
  util::SymMatrixFij d_invM;
};

} // namespace fe

#endif // FE_MASSMATRIX_H
