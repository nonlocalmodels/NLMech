// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef INP_MESHDECK_H
#define INP_MESHDECK_H

#include <string>

namespace inp {

/*! @brief Structure to read and store geometry related input data */
struct MeshDeck {

  /**
   * \defgroup Data members
   */
  /**@{*/

  /*! @brief Dimension */
  size_t d_dim;

  /*! @brief Tag for spatial discretization.
   * List of allowed values are: "", "finite_difference",
   * "weak_finite_element", "nodal_finite_element", "truss_finite_element"
   */
  std::string d_spatialDiscretization;

  /*! @brief Filename to read mesh data */
  std::string d_filename;

  /** @}*/

  /*!
   * @brief Constructor
   */
  MeshDeck() : d_dim(0){};
};

} // namespace inp
#endif // INP_MESHDECK_H
