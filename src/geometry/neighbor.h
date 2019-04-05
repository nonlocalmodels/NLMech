// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef GEOM_NEIGHBOR_H
#define GEOM_NEIGHBOR_H

#include "util/point.h"         // definition of Point3
#include <string>
#include <vector>

// forward declaration of neighbor deck
namespace inp {
struct NeighborDeck;
}

namespace geometry {

/*! @brief A class to store neighbor list and provide access to the list
 *
 * Currently, nested vector is used for list. However, this is not memory
 * efficient, as the vectors have small memory overhead and the total
 * overhead then is N times the overhead of vector, where N is the number of
 * nodes. When N is large, the total overhead becomes very large.
 *
 * @note Require further memory optimization.
 */
class Neighbor {

public:
  /*!
   * @brief Constructor
   * @param horizon Horizon
   * @param deck Input deck which contains user-specified information
   * @param nodes Pointer to nodal positions
   */
  Neighbor(const double &horizon, inp::NeighborDeck *deck,
      const std::vector<util::Point3> *nodes);

  /*!
   * @brief Get neighbor list of node i (element i incase of
   * weak_finite_element)
   *
   * @param i Id of node
   * @return Vec vector of neighboring nodes
   */
  const std::vector<size_t> &getNeighbors(const size_t &i);

  /*!
   * @brief Get pointer to full neighbor list
   * @return Vec Neighborlist
   */
  const std::vector<std::vector<size_t>> *getNeighborsP();

private:
  /*! @brief Interior flags deck */
  inp::NeighborDeck *d_neighborDeck_p;

  /*! @brief Vector of list of neighbors for each node */
  std::vector<std::vector<size_t>> d_neighbors;
};

} // namespace geometry

#endif // GEOM_NEIGHBOR_H
