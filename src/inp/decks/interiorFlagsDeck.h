////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef INP_INTERIORFLAGSDECK_H
#define INP_INTERIORFLAGSDECK_H

namespace inp {

/**
 * \ingroup Input
 */
/**@{*/

/*!
 * @brief Structure to read and store inout data for interior flags (**no-fail
 * region**)
 */
struct InteriorFlagsDeck {

  /*! @brief Flag which indicates if no-fail region is active */
  bool d_noFailActive;

  /*! @brief Flag which indicates if we compute the flag inline instead of
   * storing it. This is effective only if d_noFailActive is set to true. */
  bool d_computeAndNotStoreFlag;

  /*! @brief Tolerance to decide if the point is in interior/exterior */
  double d_noFailTol;

  /*!
   * @brief Constructor
   */
  InteriorFlagsDeck()
      : d_noFailActive(false), d_computeAndNotStoreFlag(false),
        d_noFailTol(0.){};
};

/** @}*/

} // namespace inp

#endif // INP_INTERIORFLAGSDECK_H
