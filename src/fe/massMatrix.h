// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

//
// Created by prashant on 4/1/19.
//

#ifndef FE_MASSMATRIX_H
#define FE_MASSMATRIX_H

#include "util/matrixBlaze.h" // definition of SymMatrixFij

// forward declaration
namespace inp {
struct MassMatrixDeck;
}

namespace fe {
class MassMatrix {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   */
  explicit MassMatrix(inp::MassMatrixDeck *deck);

private:
  /**
   * \defgroup Mass matrix
   */
  /**@{*/

  /*! @brief Mass matrix deck */
  inp::MassMatrixDeck *d_massMatrixDeck_p;

  /*! @brief Inverse of diagonal mass matrix stored in a vector */
  std::vector<double> d_invMDiag;

  /*! @brief Inverse of exact mass matrix stored in a Blaze matrix data type.
   * We store data in float to reduce memory load. However if results appear
   * to be not so accurate, this should be changed to double
   */
  util::SymMatrixFij d_invM;

  /** @}*/
};

} // namespace fe

#endif // FE_MASSMATRIX_H
