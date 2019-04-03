// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef INP_LOADINGDECK_H
#define INP_LOADINGDECK_H

#include <iostream> // error handling
#include <vector>

namespace inp {

/*! @brief Structure for displacement/force boundary condition data */
struct BCData {

  /**
   * @name Data members
   */
  /**@{*/

  /*!
   * @brief Type of region over which the bc will be applied
   *
   * List of allowed values are:
   *
   * - \a rectangle
   * - \a angled_rectangle
   */
  std::string d_regionType;

  /*!
   * @brief Name of the formula with respect to time
   *
   * List of allowed values are:
   * - "" (none)
   * - \a constant
   * - \a linear
   * - \a linear_step
   * - \a linear_slow_fast
   */
  std::string d_timeFnType;

  /*!
   * @brief Name of the formula of with respect to spatial coordinate
   *
   * List of allowed values are:
   * - "" (none)
   * - \a constant
   * - \a hat_x
   * - \a hat_y
   * - \a sin
   */
  std::string d_spatialFnType;

  /*! @brief X coordinate of left-bottom point of rectangle/angled rectangle */
  double d_x1;

  /*! @brief Y coordinate of left-bottom point of rectangle/angled rectangle */
  double d_y1;

  /*! @brief X coordinate of right-top point of rectangle/angled rectangle */
  double d_x2;

  /*! @brief Y coordinate of right-top point of rectangle/angled rectangle */
  double d_y2;

  /*! @brief Angle of the rectangle (for angled rectangle) */
  double d_theta;

  /*!
   * @brief List of dofs on which this bc will be applied
   *
   * E.g. if bc is only applied on x-component, d_direction will be 1. If
   * bc is applied on x- and y-component, d_direction will be vector with
   * elements 1 and 2.
   */
  std::vector<int> d_direction;

  /*! @brief List of parameters for function wrt time */
  std::vector<double> d_timeFnParams;

  /*! @brief List of parameters for function wrt spatial coordinate */
  std::vector<double> d_spatialFnParams;

  /** @}*/

  /*!
   * @brief Constructor
   */
  BCData() : d_x1(0.), d_y1(0.), d_x2(0.), d_y2(0.), d_theta(0.){};
};

/**
 * \ingroup Input
 */
/**@{*/

/*! @brief Structure to read and store policy data */
struct LoadingDeck {

  /**
   * @name Data members
   */
  /**@{*/

  /*! @brief List of displacement bcs */
  std::vector<inp::BCData> d_uBCData;

  /*! @brief List of force bcs */
  std::vector<inp::BCData> d_fBCData;

  /** @}*/

  /*!
   * @brief Constructor
   */
  LoadingDeck() = default;
};

/** @}*/

} // namespace inp

#endif // INP_LOADINGDECK_H
