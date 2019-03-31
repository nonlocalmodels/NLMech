// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef FRACTUREDECK_H
#define FRACTUREDECK_H

#include <vector>
#include <iostream>       // error handling

namespace io {

/*! @brief A structure to edge crack of any orientation */
struct EdgeCrack {

  /**
   * \defgroup Data members
   */
  /**@{*/

  /**< @brief Value 1 for orientation along horizontal axis (x-axis),
   * value -1 for orientation along vertical axis (y-axis),
   * value 0 for any crack which makes nonzero angle with the horizontal axis
   */
  int d_o;

  /**< @brief Angle the edge crack makes with the horizontal axis. I.e. 0 if
   * orientation = 1, PI/2 if orientation = -1, theta for orientation = 0
   */
  double d_theta;

  /**< @brief Total current length of crack */
  double d_l;

  /**< @brief Current length of crack on top side (right side) */
  double d_lt;

  /**< @brief Current length of crack on bottom side (left side) */
  double d_lb;

  /**< @brief Velocity of top (right) crack tip */
  std::vector<double> d_vt;

  /**< @brief Velocity of bottom (left) crack tip */
  std::vector<double> d_vb;

  /**< @brief Closest id of node to top (right) crack tip */
  int d_it;

  /**< @brief Closest id of node to bottom (left) crack tip */
  int d_ib;

  /**< @brief Closest id of node to old top (right) crack tip */
  int d_iot;

  /**< @brief Closest id of node to old bottom (left) crack tip */
  int d_iob;

  /**< @brief Top (right) crack tip location */
  std::vector<double> d_pt;

  /**< @brief Bottom (left) crack tip location */
  std::vector<double> d_pb;

  /**< @brief Old top (right) crack tip location */
  std::vector<double> d_pot;

  /**< @brief Old bottom (left) crack tip location */
  std::vector<double> d_pob;

  /** @}*/

  /*!
   * @brief Constructor
   */
  EdgeCrack() {
    d_o = 1;
    d_theta = 0.0;
    d_l = 0.0;
    d_lt = 0.0;
    d_lb = 0.0;
    d_it = -1;
    d_ib = -1;
    d_iot = -1;
    d_iob = -1;

    std::vector<double> zpt(3,0.0);
    d_vt = zpt;
    d_vb = zpt;
    d_pt = zpt;
    d_pb = zpt;
    d_pot = zpt;
    d_pob = zpt;
  }
};

/*! @brief Structure to read and store fracture data */
struct FractureDeck {

  /*! @brief Vector of pre-crack data*/
  std::vector<io::EdgeCrack> d_cracks;

  /*! @brief Flag which indicates if pre-crack is present */
  bool d_crackActive;

  /*! @brief Output interval. I.e. code perform crack data output every N
   * number of steps, when N is given by this data */
  size_t d_crackOutput;
};
} // namespace io
#endif // FRACTUREDECK_H
