////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef INP_MODELDECK_H
#define INP_MODELDECK_H

#include <string>

namespace inp {

/**
 * \ingroup Input
 */
/**@{*/

/*! @brief Structure to read and store model related input data */
struct ModelDeck {

  /*!
   * @brief Simulation type
   *
   * List of allowed values are:
   * - \a explicit
   * - \a implicit
   */
  std::string d_simType;

  /*! @brief Flag indicating if this restart problem */
  bool d_isRestartActive;

  /*!
   * @brief Tag for spatial discretization
   *
   * List of allowed values are:
   * - \a **finite_difference**
   * - \a **weak_finite_element**
   * - \a **nodal_finite_element**
   * - \a **truss_finite_element**
   */
  std::string d_spatialDiscretization;

  /*!
   * @brief Tag for time discretization
   *
   * List of allowed values are:
   * - *empty string*
   * - \a **central_difference**
   * - \a **velocity_verlet**
   */
  std::string d_timeDiscretization;

  /*! @brief Dimension */
  size_t d_dim;

  /*! @brief Final simulation time */
  double d_tFinal;

  /*! @brief Size of time steps */
  double d_dt;

  /*! @brief Number of time steps */
  size_t d_Nt;

  /*! @brief Horizon */
  double d_horizon;

  /*!
   * @brief Ratio of Horizon to mesh size
   *
   * E.g. ratio = 4 means mesh size is 1/4th of the horizon.
   */
  int d_rh;

  /*! @brief Mesh size */
  double d_h;

  /*!
   * @brief Constructor
   */
  ModelDeck()
      : d_dim(0), d_isRestartActive(false), d_tFinal(0.), d_dt(0.), d_Nt(0),
        d_horizon(0.), d_rh(0), d_h(0.){};
};

/** @}*/

} // namespace inp

#endif // INP_MODELDECK_H
