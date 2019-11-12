// Copyright (c) 2019    Patrick Diehl
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef MATERIAL_PD_ELASTIC_STATEMATERIAL_H
#define MATERIAL_PD_ELASTIC_STATEMATERIAL_H

#include "baseMaterial.h"
#include "util/point.h"
#include "util/matrixBlaze.h"
#include "geometry/neighbor.h"
#include "fe/mesh.h"
#include "geometry/volumeCorrection.h"

#include <vector>

#include <hpx/include/parallel_algorithm.hpp>

// forward declaration
namespace inp {
struct MaterialDeck;
struct OutputDeck;
}

// forward declaration
namespace data {
class DataManager;
}

// forward declaration of class
namespace fe {
class Mesh;
} // namespace fe

namespace geometry {
class VolumeCorrection;
} // namespace geometry

namespace util {
class StateBasedHelperFunctions;
} // namespace material

namespace material {

namespace pd {

/*! @brief A class implementing the state-based elastic material model
 *
 * This class implements the the state-based elastic material model
 * as described by [Silling](https://doi.org/10.1007/s10659-007-9125-1)
 * and implementation details are described by [Littlewood](https://www.sandia.gov/~djlittl/docs/PeridynamicSoftwareRoadmap.pdf).
 *
 * @see https://doi.org/10.1007/s10659-007-9125-1
 */
class ElasticState: public BaseMaterial {

public:
	/*!
	 * @brief Constructor
	 * @param deck Pointer to the input deck
	 * @param DataManager Pointer to the data manager object
	 */
	ElasticState(inp::MaterialDeck *deck,data::DataManager* dataManager);

	/*!
	 * @brief Returns energy and force state between node i and node j
	 * @param i Id of node i
	 * @param j node of j
	 * @return Value Pair of energy and force
	 */

	std::pair<util::Point3, double> getBondEF(size_t i, size_t j);

	/*!
	 * @brief Returns critical bond strain between node i and node j
	 *
	 * @param i Id of node i
	 * @param j Id of node j
	 * @return Value Critical strain
	 */
	double getSc(size_t i, size_t j);

	/*!
	 * @brief Computes the stress tensor for one node
	 * @param i Id of node i
	 * @return The stress tensor for node i
	 */
	util::Matrix33 getStress(size_t i);

	/*!
	 * @brief Computes the strain vector for one node
	 * @param i ID of node i
	 * @return The strain tensor for node I
	 */
	util::Matrix33 getStrain(size_t i);

	/*!
	 * @brief return the factor for two-dimensional problems
	 */
	double getFactor2D();

private:
	/*!
	 * @brief Computes elastic state-based material parameters from elastic constants
	 *
	 * Either Young's modulus E or bulk modulus K, and Poisson ratio \f$ \nu \f$,
	 * are needed.
	 *
	 * @param deck Input material deck
	 * @param M Moment of influence function
	 */
	void computeParameters(inp::MaterialDeck *deck, size_t dim);

	/**
	 * @name Helper functions to compute the stress and strain
	 */
	/**@{*/

	/*!
	 * @brief Computes the deformation gradient for node i
	 * @param i Id of node
	 * @return The deformation gradient for node i
	 *
	 */
	util::Matrix33 deformation_gradient(size_t i);

	/*!
	 * @brief Computes the x vector state between node i and node j
	 * @param i Id of node
	 * @param j Id of node
	 * @return The x vector state
	 *
	 */
	util::Point3 X_vector_state(size_t i, size_t j);

	/*!
	 * @brief Computes the y vector state between node i and node j
	 * @param i Id of node
	 * @param j Id of node
	 * @return The y vector state
	 *
	 */
	util::Point3 Y_vector_state(size_t i, size_t j);

	/*!
	 * @brief Computes the K shape tensor for node i
	 * @param i Id of node
	 * @return The  K shape tensor for node i
	 *
	 */
	util::Matrix33 K_shape_tensor(size_t i);

	/*!
	 * @brief Computes the K modulus tensor
	 * @param i Node i
	 * @param j Neighbor j
	 * @param k Neighbor k
	 * @param m Index for the volume correction
	 * @return The K modulus tensor
	 *
	 */
	util::Matrix33 K_modulus_tensor(size_t i, size_t j, size_t k, size_t m);

	/*!
	 * @brief Provide the image of x under the Dirac Delta Function
	 * @param deck The input deck
	 * @param problem The corresponding problem
	 * @param x Node x
	 * @param i Node i
	 * @param j Neighbor j
	 * @param k Neighbor k
	 * @param m Index for the volume correction
	 * @return 1 if x is a null-vector, otherwise 0
	 */
	double dirac_delta(util::Point3 x, size_t i, size_t j, size_t m);

	/*@}*/



	/*! @brief Correction factor for using plain stress */
	double d_factor2D;

	/*! @brief Dimension of the problem */
	size_t dim;

	double horizon;

	/*! @brief Compute strain energy */
	bool strainEnergy;


	/**
	 * @name Pointers for the function parameters
	 */
	/**@{*/

	/*! @brief Pointer to the material deck */
	const inp::MaterialDeck *d_deck;

	/*! @brief Pointer to the nodes */
	const std::vector<util::Point3> *d_nodes;

	/*! @brief Pointer to the weighted volumes of the nodes */
	const std::vector<double> *d_weightedVolume;

	/*! @brief Pointer to the volume correction of the nodes */
	const std::vector<std::vector<double>> *d_volumeCorrection;

	/*! @brief Pointer to the volumes of the nodes */
	const std::vector<double> * d_volumes;

	/*! @brief Pointer to the neighbors of the nodes */
	geometry::Neighbor *d_neighbors;

	data::DataManager* d_dataManager_p;

	/*@}*/


};

} // namespace pd

} // namespace material

#endif // MATERIAL_PD_ELASTIC_STATEMATERIAL_H
