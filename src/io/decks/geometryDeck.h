// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef GEOMETRYDECK_H
#define GEOMETRYDECK_H

#include <vector>
#include <iostream>       // error handling

namespace io {

/*! @brief Structure to read and store policy data */
struct GeometryDeck {

  /*! @brief Flag which indicates level of memory control to be enforced */
  int d_memControlFlag;
};
} // namespace io
#endif // GEOMETRYDECK_H
