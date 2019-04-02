// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef INP_MASSMATRIXDECK_H
#define INP_MASSMATRIXDECK_H

#include <string>

namespace inp {

/*! @brief Structure to read and store mass matrix related input data */
struct MassMatrixDeck {

  /**
   * \defgroup Data members
   */
  /**@{*/

  /*! @brief Mass matrix approximation type.
   * List of allowed values are: "exact" (no approximation), "lumped"
   * (lumping of mass matrix)
   */
  std::string d_MApproxType;

  /** @}*/

  /*!
   * @brief Constructor
   */
  MassMatrixDeck() = default;
};

} // namespace inp
#endif // INP_MASSMATRIXDECK_H
