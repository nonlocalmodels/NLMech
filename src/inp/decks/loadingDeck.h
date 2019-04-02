// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef LOADINGDECK_H
#define LOADINGDECK_H

#include <iostream> // error handling
#include <vector>

namespace inp {

/*! @brief Structure for displacement/force boundary condition data */
struct BCData {

  /**
   * \defgroup Data members
   */
  /**@{*/

  /*! @brief Tag specifying the type of region over which the bc will be
   * applied.
   * List of allowed values are: "rectangle" and "angled_rectangle"
   */
  std::string d_regionType;

  /*! @brief Tag which specifies the formula of boundary condition with
   * respect to time.
   * List of allowed values are: "", "constant", "linear", "linear_step",
   * "linear_slow_fast"
   */
  std::string d_timeFnType;

  /*! @brief Tag which specifies the formula of boundary condition with
   * respect to spatial coordinate.
   * List of allowed values are: "", "constant", "hat_x", "hat_y", "sin"
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

  /*! @brief Angle of the rectangle in case the region is of angled rectangle
   * type
   */
  double d_theta;

  /*! @brief List of dofs on which this bc will be applied.
   * E.g. if bc is only applied on x-component, d_direction will be 1. If
   * bc is applied on x- and y-component, d_direction will be vector with
   * elements 1 and 2
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

/*! @brief Structure to read and store policy data */
struct LoadingDeck {

  /**
   * \defgroup Data members
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

} // namespace inp
#endif // LOADINGDECK_H
