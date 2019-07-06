// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef INP_QUADRATUREDECK_H
#define INP_QUADRATUREDECK_H

namespace inp {

/**
 * \ingroup Input
 */
/**@{*/

/*! @brief Structure to read and store quadrature point related input data */
struct QuadratureDeck {

  /*! @brief Order of quadrature point integration approximation */
  size_t d_quadOrder;

  /*!
   * @brief Order of quadrature point integration approximation for
   * mass matrix
   */
  size_t d_quadOrderM;

  /*!
   * @brief Constructor
   */
  QuadratureDeck() : d_quadOrder(0), d_quadOrderM(0){};
};

/** @}*/

} // namespace inp

#endif // INP_QUADRATUREDECK_H
