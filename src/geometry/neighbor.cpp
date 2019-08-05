////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "neighbor.h"
#include "inp/decks/neighborDeck.h"
#include "util/compare.h"
#include <hpx/include/parallel_algorithm.hpp>

geometry::Neighbor::Neighbor(const double &horizon, inp::NeighborDeck *deck,
                             const std::vector<util::Point3> *nodes)
    : d_neighborDeck_p(deck) {

  d_neighbors.resize(nodes->size());

  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      nodes->size(), [this, horizon, nodes](boost::uint64_t i) {
        util::Point3 xi = (*nodes)[i];

        // loop over all the nodes and check which nodes are
        // within the horizon ball of i_node
        std::vector<size_t> neighs;
        for (size_t j = 0; j < nodes->size(); j++) {

          if (j == i)
            continue;

          if (util::compare::definitelyLessThan(xi.dist((*nodes)[j]),
                                                horizon + 1.0E-10))
            neighs.push_back(j);
        } // loop over nodes j

        this->d_neighbors[i] = neighs;
      }); // end of parallel for loop

  f.get();
}

const std::vector<size_t> &geometry::Neighbor::getNeighbors(const size_t &i) {
  return d_neighbors[i];
}
const std::vector<std::vector<size_t>> *geometry::Neighbor::getNeighborsP() {
  return &d_neighbors;
}
