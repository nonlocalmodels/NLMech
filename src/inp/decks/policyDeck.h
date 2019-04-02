// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef POLICYDECK_H
#define POLICYDECK_H

#include <vector>

namespace inp {

/*! @brief Structure to read and store policy related input data */
struct PolicyDeck {

  /**
   * \defgroup Data members
   */
  /**@{*/

  /*! @brief Flag which indicates level of memory control to be enforced */
  int d_memControlFlag;

  /** @}*/

  /*!
   * @brief Constructor
   */
  PolicyDeck() : d_memControlFlag(0){};
};

} // namespace inp
#endif // POLICYDECK_H
