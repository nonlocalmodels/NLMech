////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

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
 * overhead is N times the overhead of vector, where N is the number of
 * nodes. When N is large, the total overhead becomes very large.
 *
 * @todo Require further memory optimization.
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
   * @brief Get neighbor list of node i (element i in case of
   * **weak_finite_element**)
   *
   * @param i Id of node
   * @return vector Vector of neighboring nodes
   */
  const std::vector<size_t> &getNeighbors(const size_t &i);

  /*!
   * @brief Get the pointer to full neighbor list
   * @return pointer Pointer of neighbor list data
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
