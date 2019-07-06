////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef INP_MESHDECK_H
#define INP_MESHDECK_H

#include <string>

namespace inp {

/**
 * \ingroup Input
 */
/**@{*/

/*! @brief Structure to read and store mesh related input data */
struct MeshDeck {

  /**
   * @name Data members
   */
  /**@{*/

  /*! @brief Dimension */
  size_t d_dim;

  /*!
   * @brief Tag for spatial discretization
   *
   * List of allowed values are:
   * - \a finite_difference
   * - \a weak_finite_element
   * - \a nodal_finite_element
   * - \a truss_finite_element
   */
  std::string d_spatialDiscretization;

  /*! @brief Filename to read mesh data */
  std::string d_filename;

  /*! @brief Flag which indicates if mesh size is to be computed */
  bool d_computeMeshSize;

  /*! @brief Mesh size */
  double d_h;

  /** @}*/

  /*!
   * @brief Constructor
   */
  MeshDeck() : d_dim(0), d_computeMeshSize(false), d_h(0.){};
};

/** @}*/

} // namespace inp

#endif // INP_MESHDECK_H
