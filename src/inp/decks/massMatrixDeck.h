// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef INP_MASSMATRIXDECK_H
#define INP_MASSMATRIXDECK_H

#include <string>

namespace inp {

/**
 * \ingroup Input
 */
/**@{*/

/*! @brief Structure to read and store mass matrix related input data */
struct MassMatrixDeck {

  /**
   * @name Data members
   */
  /**@{*/

  /*!
   * @brief Mass matrix approximation type
   *
   * List of allowed values are:
   * - \a exact -- no approximation
   * - \a lumped -- lumping of mass matrix
   */
  std::string d_MApproxType;

  /** @}*/

  /*!
   * @brief Constructor
   */
  MassMatrixDeck() = default;
};

/** @}*/

} // namespace inp

#endif // INP_MASSMATRIXDECK_H
