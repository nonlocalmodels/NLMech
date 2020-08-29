////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef MODEL_QSMODEL_H
#define MODEL_QSMODEL_H

#include <hpx/config.hpp>
#include "model/model.h"
#include "inp/decks/modelDeck.h"

// forward declaration of class
namespace fe {
class MassMatrix;
class Mesh;
// class Quadrature;
}// namespace fe

namespace geometry {
class Fracture;
class InteriorFlags;
class Neighbor;
class VolumeCorrection;
} // namespace geometry

namespace inp {
struct ModelDeck;
struct RestartDeck;
struct OutputDeck;
struct SolverDeck;
class Input;
class Policy;
} // namespace inp

namespace loading {
class InitialCondition;
class ULoading;
class FLoading;
struct BCData;
} // namespace loading

namespace material {
namespace pd {
class BaseMaterial;
}
} // namespace material

namespace util {
class StateBasedHelperFunctions;
} // namespace util

namespace data {
class DataManager;
} // namespace data

namespace rw {
namespace writer {
class Writer;
}
}

namespace model {

/**
 * \defgroup Explicit Explicit
 */
/**@{*/

/*! @brief A class for *finite difference approximation* of **Peridynamics**
 *
 * We consider *explicit* schemes such as *central difference* and
 * *velocity verlet* for time integration.
 *
 * This class acts as a holder of lower rank classes, such as Mesh, Loading,
 * InitialCondition, Fracture, etc, and uses the methods and data of the
 * lower rank classes to perform calculation.
 *
 * @note 1. We can run finite difference on any finite element mesh as long as
 * the mesh consists of only one type of elements. We can mesh the domain
 * using **Gmsh** and use its **.msh** file to run the finite difference
 * approximation.
 *
 * @note 2. Currently only dimension 2 is supported.
 *
 * @note 3. Either triangle or quadrangle elements are supported.
 */
template<class T>
class QuasiStaticModel: public Model {

public:

	/*! @brief Constructor
	 *  @param deck Pointer to the input deck
	 */
	QuasiStaticModel(inp::Input *deck);

	/*! @brief Destructor */
	~QuasiStaticModel();

private:

	/*!
	 * @brief Initialize all data members
	 */
	void initHObjects();

	/*!
	 * @brief Computes the forces of all nodes
	 * @param full If true the Strain and Stress tensors are computed
	 */
	void computeForces(bool full=false);

	/*!
	 * @brief Computes the forces of all nodes using the pertubated displacement
	 * @param thread The thread which is doing the acutal computation
	 */
	void computePertubatedForces(size_t thread); 

	/*! @brief Assembles the Jacobian matrix
	 */
	void assembly_jacobian_matrix();

	/*! @brief Assembles the Jacobian matrix
	 * @param begin First node of the chunk
	 * @param end Last node of the chunk
	 * @param thread Id of the thread handling this chunk
	 */
	void assembly_jacobian_matrix_part(size_t begin, size_t end, size_t thread);

	/*! @brief Computes the new displacement of Newton step
	 * @param res Residual vector
	 * @return The updated displacement
	 */
	util::VectorXi newton_step(util::VectorXi &res);

	/*!
	 * @brief Starts the simulation and controls the solver
	 */
	void solver();

	/*!
	 * @brief Computes the residual for the Newton step
	 * @return The residual vector
	 */
	util::VectorXi computeResidual();

	/**
	 * @name Functions to manipulate the tangent stiffness matrix
	 *
	 */
	/**@{*/

	/*!
	 * @brief Removes the i-th row of a matrix
	 * @param matrix The matrix
	 * @param rowToRemove Id of the row to remove
	 */
	void removeRow(util::Matrixij &matrix, size_t rowToRemove);

	/*!
	 * @brief Removes the i-th column of a matrix
	 * @param matrix The matrix
	 * @param colToRemove Id of the column to remove
	 */
	void removeCol(util::Matrixij &matrix, size_t colToRemove);

	/*!
	 * @brief Removes the i-th row of a vector
	 * @param vector The vector
	 * @param rowToRemove Id of the row to remove
	 */
	void removeRow(util::VectorXi &vector, size_t rowToRemove);

	/** @}*/

	/*! @brief Number of nodes */
	size_t d_nnodes;

	/*! @brief Current simulation time */
	double d_time;

	/*! @brief Name of the model */
	std::string d_name;

	/*! @brief Number of available os threads */
	size_t d_osThreads;

	/*! Jacobian matrix */
	util::Matrixij jacobian;

	/*! @brief Data manager objects for the assembly of the stiffness matrix */
	std::vector<data::DataManager*> d_dataManagers;

	/*! @brief Model deck */
	inp::ModelDeck *d_modelDeck_p;

	/*! @brief Output deck */
	inp::OutputDeck *d_outputDeck_p;

	/*! @brief Pointer to Input object */
	inp::Input *d_input_p;

	/*! @brief Pointer to Material object */
	material::pd::BaseMaterial *d_material_p;

	/*! @brief Data Manager */
	data::DataManager *d_dataManager_p;

};

}

#endif
