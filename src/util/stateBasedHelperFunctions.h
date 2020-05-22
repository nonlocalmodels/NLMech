////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTIL_STATEBASEHELPER_H
#define UTIL_STATEBASEHELPER_H

#include "util/point.h"         // definition of Point3
#include <string>
#include <vector>
#include "fe/mesh.h"
#include "geometry/neighbor.h"
#include "data/DataManager.h"
#include "inp/decks/modelDeck.h"
#include "geometry/volumeCorrection.h"


// forward declaration of class
namespace fe {
class Mesh;
} // namespace fe

// forward declaration of class
namespace data {
class DataManager;
} // namespace data

namespace inp {
struct ModelDeck;
}

namespace geometry {
class VolumeCorrection;
} // namespace geometry

namespace util {

/*! @brief A class for global methods of state-based material models
 *
 * This class implements the computation of the extension and the dilatation
 * of a state-based material models as described by [Silling](https://doi.org/10.1007/s10659-007-9125-1)
 * and implementation details are described by [Littlewood](https://www.sandia.gov/~djlittl/docs/PeridynamicSoftwareRoadmap.pdf).
 *
 *@see https://doi.org/10.1007/s10659-007-9125-1
 */
class StateBasedHelperFunctions {

public:

	/*!
	 * @brief Constructor
	 * @param DataManager Pointer to the data manager object
	 */
	StateBasedHelperFunctions(data::DataManager* dataManager,
			double factor);
private:
	/*!
	 * @brief Computes the extension \f$ e_i = \vert \eta + \xi \vert - \vert \xi \vert \f$ and
	 * dilatation \f$ \sum\limits_{B_\delta(x_i)} \frac{3}{m_i} \vert \xi \vert e_i V_j \f$ for each node where $m_i$ is the weighted volume.
	 * @param dataMamager Class holding all the global simulaiton data
	 * @param dim Dimension of the problem
	 * @param factor Dimensional depended material property
	 */
	void dilatation(data::DataManager *dataManager, size_t dim, double factor);

};

} // namespace util

#endif //UTIL_STATEBASEHELPER_H
