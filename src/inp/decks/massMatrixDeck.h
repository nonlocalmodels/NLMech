////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

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
  /*!
   * @brief Mass matrix approximation type
   *
   * List of allowed values are:
   * - \a exact -- no approximation
   * - \a lumped -- lumping of mass matrix
   */
  std::string d_MApproxType;

  /*!
   * @brief Constructor
   */
  MassMatrixDeck() = default;
};

/** @}*/

} // namespace inp

#endif // INP_MASSMATRIXDECK_H
