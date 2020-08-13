////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef GEOM_VOLUMECORRECTION_H
#define GEOM_VOLUMECORRECTION_H

#include "util/point.h"         // definition of Point3
#include <string>
#include <vector>
#include "neighbor.h"
#include "fe/mesh.h"
#include "data/DataManager.h"
#include "inp/decks/modelDeck.h"

// forward declaration of class
namespace data {
class DataManager;
} // namespace data

namespace inp {
struct ModelDeck;
}

// forward declaration of class
namespace fe {
class Mesh;
} // namespace fe

namespace geometry {

/*! @brief A class to to compute the volume correction for state-based models
 *
 * This class implements the so-called volume correction for state-based PD models
 * using the approach described by [Silling](https://doi.org/10.1007/s10659-007-9125-1)
 * and implementation details are described by [Littlewood](https://www.sandia.gov/~djlittl/docs/PeridynamicSoftwareRoadmap.pdf).
 *
 * @see https://doi.org/10.1007/s10659-007-9125-1
 */
class VolumeCorrection {

public:

	/*! @brief Weighted volume of nodes */
	std::vector<double>* d_weightedVolume_p;
	/*! @brief Volume correction for the neighborhood of each node */
	std::vector<std::vector<double>>* d_volumeCorrection_p;

	/*!
	 * @brief Constructor
	 * @param dataManager Pointer to the data manager object
	 */
	VolumeCorrection(data::DataManager* dataManager);

private:

	/*! @brief Computes the volume correction for each node inside the neighborhood. Since
	 * nodes close the boundary of the neighborhood have less volume inside the neighborhood
	 * and should have less influence while computing the forces.
	 * @param horizon Horizon
	 * @param dx Nodal spacing
	 * @param neighbors Class containing the neighborhoods
	 * @param nodes Pointer to nodal positions
	 * */
	//std:string d_type;
	void correctVolume(const double &horizon, const double &dx,
			geometry::Neighbor *neighbors,
			const std::vector<util::Point3> *nodes);

	/*! @brief Computes the weighted volume \f$m_i = \sum\limits_{j\in B_\delta(x_i)} w \vert x_j - x_i \vert V_j\f$
	 * with \f$w\f$ as the influence function for each discrete node.
	 * @param neighbors Class containing the neighborhoods
	 * @param nodes Pointer to nodal positions
	 * @param p_mesh Pointer to the mesh object
	 * */
	void weightedVolume(geometry::Neighbor *neighbors,
			const std::vector<util::Point3> *nodes, fe::Mesh* p_mesh);

};

} // namespace geometry

#endif //GEOM_VOLUMECORRECTION_H
