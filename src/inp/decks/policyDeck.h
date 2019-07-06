////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef INP_POLICYDECK_H
#define INP_POLICYDECK_H

#include <vector>

namespace inp {

/**
 * \ingroup Input
 */
/**@{*/

/*! @brief Structure to read and store policy related input data */
struct PolicyDeck {

  /**
   * @name Data members
   */
  /**@{*/

  /*!
   * @brief Flag which indicates level of memory control to be enforced
   *
   * Default is 0 which means no control. Max at present is 2 which means as
   * much control as possible.
   */
  int d_memControlFlag;

  /*!
   * @brief Enable post-processing calculation
   *
   * Default is true.
   */
  bool d_enablePostProcessing;

  /** @}*/

  /*!
   * @brief Constructor
   */
  PolicyDeck() : d_memControlFlag(0), d_enablePostProcessing(true){};
};

/** @}*/

} // namespace inp

#endif // INP_POLICYDECK_H
