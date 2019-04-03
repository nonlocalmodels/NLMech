// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef INP_INTERIORFLAGSDECK_H
#define INP_INTERIORFLAGSDECK_H

namespace inp {

/**
 * \ingroup Input
 */
/**@{*/

/*! @brief Structure to read and store interior flags (no-fail region) */
struct InteriorFlagsDeck {

  /**
   * @name Data members
   */
  /**@{*/

  /*! @brief Flag which indicates if no-fail region is active */
  bool d_noFailActive;

  /*! @brief Tolerance to decide if the point is in interior/exterior */
  double d_noFailTol;

  /** @}*/

  /*!
   * @brief Constructor
   */
  InteriorFlagsDeck() : d_noFailActive(false), d_noFailTol(0.){};
};

/** @}*/

} // namespace inp

#endif // INP_INTERIORFLAGSDECK_H
