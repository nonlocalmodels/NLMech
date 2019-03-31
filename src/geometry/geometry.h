// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <string>

// forward declaration of geometry deck
namespace io {
struct GeometryDeck;
}

//!Collection of methods and database purely related to geometry
namespace geometry {

/*! @brief Methods and database associated to the mesh */
class Geometry {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   */
  Geometry(io::GeometryDeck *deck);

private:
  /**
   * \defgroup Mesh related data
   */
  /**@{*/

  /*! @brief Vector of initial (reference) coordinates of nodes */
  std::vector<double> d_nd;

  /*! @brief Element-node connectivity data. I.e. list of ids of nodes which
   * are also vertex of the element */
  std::vector<size_t> d_enc;

  /*! @brief Node-element connectivity data. I.e. list of ids of elements
   * which have particular node as its vertex */
  std::vector<size_t> d_nec;

  /*! @brief Vector of fixity mask of each node */
  std::vector<bool> d_fix;

  /*! @brief Vector of volume of each node. For uniform square mesh, the
   * volume is simply h^2 (h^3 in 3-d), however, for general mesh, the volume
   * is computed using linear interpolation
   * */
  std::vector<double> d_vol;

  /** @}*/

};

} // namespace geometry

#endif // GEOMETRY_H
