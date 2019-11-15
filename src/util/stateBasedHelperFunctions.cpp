////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "stateBasedHelperFunctions.h"

#include <hpx/include/parallel_algorithm.hpp>

#include "util/compare.h"

util::StateBasedHelperFunctions::StateBasedHelperFunctions(
    data::DataManager *dataManager, double factor) {
      
      dataManager->setExtensionP(new std::vector<std::vector<double>>(
      dataManager->getMeshP()->getNodesP()->size()));

  dilatation(dataManager,dataManager->getModelDeckP()->d_dim, factor);

}

void util::StateBasedHelperFunctions::dilatation(
    data::DataManager *dataManager, size_t dim, double factor) {

   dataManager->setDilatationP(new std::vector<double>(dataManager->getMeshP()->getNodesP()->size(), 0.));

  

   hpx::parallel::for_loop(hpx::parallel::execution::par, 0, dataManager->getMeshP()->getNodesP()->size(),
  	[&](boost::uint64_t i) {

    size_t k = 0;
    double w = 1;
    for (auto j : dataManager->getNeighborP()->getNeighbors(i)) {
      util::Point3 Y = ((*dataManager->getMeshP()->getNodesP())[j] + (*dataManager->getDisplacementP())[j]) -
                       ((*dataManager->getMeshP()->getNodesP())[i] + (*dataManager->getDisplacementP())[i]);
      util::Point3 X = (*dataManager->getMeshP()->getNodesP())[j] - (*dataManager->getMeshP()->getNodesP())[i];

      (*dataManager->getExtensionP())[i].push_back(Y.length() - X.length());  

      switch (dim) {
        case 1:
          (*dataManager->getDilatationP())[i] += (1. / (*dataManager->getVolumeCorrectionP()->d_weightedVolume_p)[i]) * w * X.length() *
                                  (*dataManager->getExtensionP())[i][k] *
                                  (*dataManager->getVolumeCorrectionP()->d_volumeCorrection_p)[i][k] * (*dataManager->getMeshP()->getNodalVolumeP())[j];
          break;
        case 2:
          (*dataManager->getDilatationP())[i] += (2. / (*dataManager->getVolumeCorrectionP()->d_weightedVolume_p)[i]) * factor * w *
                                  X.length() *
                                  (*dataManager->getExtensionP())[i][k] *
                                  (*dataManager->getVolumeCorrectionP()->d_volumeCorrection_p)[i][k] * (*dataManager->getMeshP()->getNodalVolumeP())[j];
          break;
        case 3:
          (*dataManager->getDilatationP())[i] += (3. / (*dataManager->getVolumeCorrectionP()->d_weightedVolume_p)[i]) * w * X.length() *
                                  (*dataManager->getExtensionP())[i][k] *
                                  (*dataManager->getVolumeCorrectionP()->d_volumeCorrection_p)[i][k] * (*dataManager->getMeshP()->getNodalVolumeP())[j];
          break;
      }

      k++;
    }

    });   
}
