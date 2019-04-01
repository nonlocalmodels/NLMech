// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef GEOMETRYDECK_H
#define GEOMETRYDECK_H

#include <string>

namespace inp {

/*! @brief Structure to read and store geometry related input data */
struct GeometryDeck {

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

  /*! @brief Order of quadrature point integration approximation */
  size_t d_quadOrder;

  /*! @brief Order of quadrature point integration approximation for
   * computation of mass matrix
   */
  size_t d_quadOrderM;

  /*! @brief Mass matrix approximation type.
   * List of allowed values are: "exact" (no approximation), "lumped"
   * (lumping of mass matrix)
   */
  std::string d_MApproxType;

  /** @}*/

  /*!
   * @brief Constructor
   */
  GeometryDeck() : d_dim(0), d_quadOrder(0), d_quadOrderM(0){};
};

} // namespace inp
#endif // GEOMETRYDECK_H
