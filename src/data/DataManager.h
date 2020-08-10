////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
///////////////////////////////////////////////////////////////////////////////

#ifndef SRC_DATA_DATAMANAGER_H_
#define SRC_DATA_DATAMANAGER_H_

#include "util/point.h"

#include <vector>

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
class Input;
class Policy;
} // namespace inp

namespace loading {
class InitialCondition;
class ULoading;
class FLoading;
} // namespace loading

namespace material {
namespace pd {
class BaseMaterial;
}
} // namespace material

namespace util {
class StateBasedHelperFunctions;
} // namespace material

namespace data {

static size_t instances = 0;

class DataManager {

public:

	DataManager();

	/**
	 * @name Access to the deck objects
	 *
	 */
	/**@{*/

	void setModelDeckP(inp::ModelDeck *pointer);

	inp::ModelDeck* getModelDeckP();

	void setOutputDeckP(inp::OutputDeck *pointer);

	inp::OutputDeck* getOutputDeckP();

	/** @}*/

	/**
	 * @name Major simulation data
	 *
	 */
	/**@{*/

	void setBodyForceP(std::vector<util::Point3> *pointer);

	std::vector<util::Point3>* getBodyForceP();

	void setForceP(std::vector<util::Point3> *pointer);

	std::vector<util::Point3>* getForceP();


	void setReactionForceP(std::vector<util::Point3> *pointer);

	std::vector<util::Point3>* getReactionForceP();


	void setVelocityP(std::vector<util::Point3> *pointer);

	std::vector<util::Point3>* getVelocityP();

	void setDisplacementP(std::vector<util::Point3> *pointer);

	std::vector<util::Point3>* getDisplacementP();

	/** @}*/

	void setMeshP(fe::Mesh *d_mesh_p);

	fe::Mesh* getMeshP();

	void setNeighborP(geometry::Neighbor *d_neighbor_p);

	geometry::Neighbor* getNeighborP();

	void setVolumeCorrectionP(geometry::VolumeCorrection* pointer);

	geometry::VolumeCorrection* getVolumeCorrectionP();


	void setStateBasedHelperFunctionsP(util::StateBasedHelperFunctions* pointer);

	util::StateBasedHelperFunctions* getStateBasedHelperFunctionsP();

	void setDisplacementLoadingP(loading::ULoading * pointer);

	loading::ULoading* getDisplacementLoadingP();


	void setForceLoadingP(loading::FLoading * pointer);

	loading::FLoading* getForceLoadingP();


	void setExtensionP(std::vector<std::vector<double>>* pointer);

	std::vector<std::vector<double>>* getExtensionP();


	void setStressTensorP(std::vector<util::Matrix33>* pointer);

	std::vector<util::Matrix33>* getStressTensorP();

	void setStrainTensorP(std::vector<util::Matrix33>* pointer);

	std::vector<util::Matrix33>* getStrainTensorP();

	void setDilatationP(std::vector<double> * pointer);

	std::vector<double> * getDilatationP();

	void setTotalReactionForceP(std::vector<double>* pointer);

	std::vector<double>* getTotalReactionForceP();

	void setStrainEnergyP(std::vector<float>* pointer);

	std::vector<float>* getStrainEnergyP();

	void setWorkDoneP(std::vector<float>* pointer);

	std::vector<float>* getWorkDoneP();

	void setPhiP(std::vector<float>* pointer);

	std::vector<float>* getPhiP();

	void setDamageFunctionP(std::vector<float>* pointer);

	std::vector<float>* getDamageFunctionP();

private:

	/**
	 * @name Pointers to the deck objects
	 *
	 */
	/**@{*/

	/*! @brief Model deck */
	inp::ModelDeck *d_modelDeck_p = nullptr;

	/*! @brief Output deck */
	inp::OutputDeck *d_outputDeck_p = nullptr;

	/** @}*/

	/**
	 * @name Pointers to the major simulation data
	 *
	 */
	/**@{*/

	/*! @brief Pointer to the body force vector */
	std::vector<util::Point3> *d_b_p = nullptr;

	/*! @brief Pointer to the displacement vector */
	std::vector<util::Point3> *d_u_p = nullptr;

	/*! @brief Pointer to the velocity vector */
	std::vector<util::Point3> *d_v_p = nullptr;

	/*! @brief Pointer to the force vector */
	std::vector<util::Point3> *d_f_p = nullptr;

	/*! @brief Pointer to the reaction force */
	std::vector<util::Point3> *d_reaction_force_p = nullptr;

	/*! @brief Extension for the neighborhood of each node */
	std::vector<std::vector<double>>* d_extension_p;

	/*! @brief Dilatation of nodes */
	std::vector<double> *d_dilatation_p = nullptr;

	/*! @brief Pointer to the strain energy vector */
	std::vector<float>* d_e_p = nullptr;

	/*! @brief Pointer to the total reaction force vector */
	std::vector<double>* d_total_reaction_force_p = nullptr;

	/*! @brief Pointer to the Work done on each of the nodes */
  	std::vector<float>* d_w_p = nullptr;

	 /*! @brief Damage function \f$ \phi \f$ at the nodes */
  	std::vector<float>* d_phi_p = nullptr;

	  /*! @brief Damage function \f$ Z \f$ at the nodes */
  	std::vector<float>* d_Z_p = nullptr;

	/*! @brief Pointer to the strain tensor vector */
	std::vector<util::Matrix33>* d_strain_p = nullptr;

	/*! @brief Pointer to the stress tensor vector */
    std::vector<util::Matrix33>* d_stress_p = nullptr;


	


	/**@{*/

	/*! @brief Pointer to Mesh object */
	fe::Mesh *d_mesh_p = nullptr;

	/*! @brief Pointer to Neighbor object */
	geometry::Neighbor *d_neighbor_p = nullptr;

	/*! @brief Pointer to Volume Correction object */
	geometry::VolumeCorrection *d_volumeCorrection_p = nullptr;

	/*! @brief Pointer to StatebasedHelper function object */
	util::StateBasedHelperFunctions* d_sbhelper_p = nullptr;


	/*! @brief Pointer to displacement Loading object */
	loading::ULoading *d_uLoading_p = nullptr;

	/*! @brief Pointer to force Loading object */
	loading::FLoading *d_fLoading_p = nullptr;

};

}

#endif /* SRC_DATA_DATAMANAGER_H_ */
