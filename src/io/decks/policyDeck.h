// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef POLICYDECK_H
#define POLICYDECK_H

#include <vector>
#include <iostream>       // error handling

namespace io {

/*! @brief Structure to read and store policy data */
struct PolicyDeck {

  /*! @brief Flag which indicates level of memory control to be enforced */
  int d_memControlFlag;
};
} // namespace io
#endif // POLICYDECK_H
