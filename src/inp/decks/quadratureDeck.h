////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef INP_QUADRATUREDECK_H
#define INP_QUADRATUREDECK_H

namespace inp {

/**
 * \ingroup Input
 */
/**@{*/

/*! @brief Structure to read and store quadrature point related input data */
struct QuadratureDeck {

  /**
   * @name Data members
   */
  /**@{*/

  /*! @brief Order of quadrature point integration approximation */
  size_t d_quadOrder;

  /*!
   * @brief Order of quadrature point integration approximation for
   * mass matrix
   */
  size_t d_quadOrderM;

  /** @}*/

  /*!
   * @brief Constructor
   */
  QuadratureDeck() : d_quadOrder(0), d_quadOrderM(0){};
};

/** @}*/

} // namespace inp

#endif // INP_QUADRATUREDECK_H
