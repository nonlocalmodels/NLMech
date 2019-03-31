// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef INTERIORFLAGS_H
#define INTERIORFLAGS_H

#include <vector>
#include <string>

// forward declaration of interior flags deck
namespace io {
struct InteriorFlagsDeck;
}

//!Collection of methods and database purely related to geometry
namespace geometry {

/*! @brief Methods and database associated to the mesh */
class InteriorFlags {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   */
  InteriorFlags(io::InteriorFlagsDeck *deck);

private:
  /**
   * \defgroup Mesh related data
   */
  /**@{*/

  /*! @brief Vector of initial (reference) coordinates of nodes */
//  std::vector<double> d_nd;

  /** @}*/

};

} // namespace geometry

#endif // INTERIORFLAGS_H
