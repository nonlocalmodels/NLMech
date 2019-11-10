////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "volumeCorrection.h"

#include <hpx/include/parallel_algorithm.hpp>

#include "util/compare.h"

geometry::VolumeCorrection::VolumeCorrection(data::DataManager *dataManager) {
  correctVolume(dataManager->getModelDeckP()->d_horizon,
                dataManager->getMeshP()->getMeshSize(),
                dataManager->getNeighborP(),
                dataManager->getMeshP()->getNodesP());
  weightedVolume(dataManager->getNeighborP(),
                 dataManager->getMeshP()->getNodesP(), dataManager->getMeshP());
}

void geometry::VolumeCorrection::correctVolume(
    const double &horizon, const double &dx, geometry::Neighbor *neighbors,
    const std::vector<util::Point3> *nodes) {
  d_volumeCorrection_p = new std::vector<std::vector<double>>(nodes->size());

  hpx::parallel::for_loop(
      hpx::parallel::execution::par, 0, nodes->size(), [&](boost::uint64_t i) {
        (*d_volumeCorrection_p)[i] =
            std::vector<double>(neighbors->getNeighbors(i).size(), 1.);

        size_t k = 0;
        for (auto j : neighbors->getNeighbors(i)) {
          util::Point3 X = (*nodes)[j] - (*nodes)[i];

          double r = dx * 0.5;
          if (util::compare::definitelyGreaterThan(X.length(), horizon - r))
            (*d_volumeCorrection_p)[i][k] = (horizon + r - X.length()) / dx;

          k++;
        }
      });
}

void geometry::VolumeCorrection::weightedVolume(
    geometry::Neighbor *neighbors, const std::vector<util::Point3> *nodes,
    fe::Mesh *p_mesh) {
  d_weightedVolume_p = new std::vector<double>(nodes->size(), 1);

  hpx::parallel::for_loop(
      hpx::parallel::execution::par, 0, nodes->size(), [&](boost::uint64_t i) {
        double tmp = 0;
        size_t k = 0;

        for (auto j : neighbors->getNeighbors(i)) {
          util::Point3 X = (*nodes)[j] - (*nodes)[i];
          tmp += (X.length() * X.length() * (*d_volumeCorrection_p)[i][k] *
                  p_mesh->getNodalVolume(j));
          k++;
        }

        (*d_weightedVolume_p)[i] = tmp;
      });
}
