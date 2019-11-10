////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "stateBasedHelperFunctions.h"
#include "util/compare.h"
#include <hpx/include/parallel_algorithm.hpp>

util::StateBasedHelperFunctions::StateBasedHelperFunctions(
		data::DataManager *dataManager, double factor) {

	dataManager->setExtensionP(
			new std::vector<std::vector<double>>(
					dataManager->getMeshP()->getNodesP()->size()));

	dilatation(dataManager,dataManager->getNeighborP(),
			dataManager->getMeshP()->getNodesP(),
			dataManager->getDisplacementP(),
			dataManager->getVolumeCorrectionP()->d_weightedVolume_p,
			dataManager->getVolumeCorrectionP()->d_volumeCorrection_p,
			dataManager->getMeshP()->getNodalVolumeP(),
			dataManager->getModelDeckP()->d_dim, factor);

}

void util::StateBasedHelperFunctions::dilatation(data::DataManager *dataManager,geometry::Neighbor *neighbors,
		const std::vector<util::Point3> *nodes,
		const std::vector<util::Point3> *displacement,
		const std::vector<double> *weightedVolume,
		const std::vector<std::vector<double>> *volumeCorrection,
		const std::vector<double> *volumes, size_t dim, double factor) {

	d_dilatation_p = new std::vector<double>(nodes->size(), 0.);
	//d_extension_p = new std::vector<std::vector<double>>(nodes->size());
	//d_extension_p = std::make_shared<std::vector<std::vector<double>>>();

	//hpx::parallel::for_loop(hpx::parallel::execution::par, 0, nodes->size(),
	//	[&](boost::uint64_t i) {

	for (size_t i = 0; i < nodes->size(); i++) {
		size_t k = 0;
		double w = 1;
		for (auto j : neighbors->getNeighbors(i)) {

			util::Point3 Y = ((*nodes)[j] + (*displacement)[j])
					- ((*nodes)[i] + (*displacement)[i]);
			util::Point3 X = (*nodes)[j] - (*nodes)[i];

			(*dataManager->getExtensionP())[i].push_back(Y.length() - X.length()); //

			switch (dim) {
			case 1:
				(*d_dilatation_p)[i] += (1. / (*weightedVolume)[i]) * w
						* X.length() * (*dataManager->getExtensionP())[i][k]
						* (*volumeCorrection)[i][k] * (*volumes)[j];
				break;
			case 2:
				(*d_dilatation_p)[i] += (2. / (*weightedVolume)[i]) * factor * w
						* X.length() * (*dataManager->getExtensionP())[i][k]
						* (*volumeCorrection)[i][k] * (*volumes)[j];
				break;
			case 3:
				(*d_dilatation_p)[i] += (3. / (*weightedVolume)[i]) * w
						* X.length() * (*dataManager->getExtensionP())[i][k]
						* (*volumeCorrection)[i][k] * (*volumes)[j];
				break;
			}

			k++;

		}

		//});

	}

}

