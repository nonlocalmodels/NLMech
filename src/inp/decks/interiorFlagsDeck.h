// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef INTERIORFLAGSDECK_H
#define INTERIORFLAGSDECK_H

namespace inp {

/*! @brief Structure to read and store interior flags (no-fail region) related
 * input data
 */
struct InteriorFlagsDeck {

  /**
   * \defgroup Data members
   */
  /**@{*/

  /*! @brief Flag which indicates if no-fail region is active */
  bool d_noFailActive;

  /*! @brief Tolerance to compare which nodes fall in interior and which
   * fall in exterior
   */
  double d_noFailTol;

  /** @}*/

  /*!
   * @brief Constructor
   */
  InteriorFlagsDeck() : d_noFailActive(false), d_noFailTol(0.){};
};

} // namespace inp
#endif // INTERIORFLAGSDECK_H
