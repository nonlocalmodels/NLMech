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

/*! @brief Data mamanger to share the global simulation data between the classes */
namespace data {

/*! @brief Data manager to collect the global simulation data */
class DataManager {

public:

	/*! @brief Constructor */
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
	 * @name Access methods
	 *
	 */
	/**@{*/

	/*! Sets the body force pointer 
	* @param pointer Pointer
	*/
	void setBodyForceP(std::vector<util::Point3> *pointer);

	/*! Get the pointer to the body force
	 * @return pointer
	 */
	std::vector<util::Point3>* getBodyForceP();

	/*! Sets the force pointer 
	* @param pointer Pointer
	*/
	void setForceP(std::vector<util::Point3> *pointer);

	/*! Get the pointer to force
	 * @return pointer
	 */
	std::vector<util::Point3>* getForceP();

	/*! Sets the body reaction force pointer 
	* @param pointer Pointer
	*/
	void setReactionForceP(std::vector<util::Point3> *pointer);

	/*! Get the pointer to reaction force
	 * @return pointer
	 */
	std::vector<util::Point3>* getReactionForceP();

	/*! Sets the  velocity pointer 
	* @param pointer Pointer
	*/
	void setVelocityP(std::vector<util::Point3> *pointer);

    /*! Get the pointer to velocity 
    *@return pointer
    */
	std::vector<util::Point3>* getVelocityP();

	/*! Sets the displacement pointer 
	 * @param pointer Pointer
	 */
	void setDisplacementP(std::vector<util::Point3> *pointer);


	/*! Get the pointer to displacement
	 * @return pointer
	 */
	std::vector<util::Point3>* getDisplacementP();

	/*! Sets the mesh pointer 
	 * @param pointer Pointer
	 */
	void setMeshP(fe::Mesh *d_mesh_p);


	/*! Get the pointer to the mesh object
	 * @return pointer
	 */
	fe::Mesh* getMeshP();

	/*! Sets the neighborhood pointer
	* @param pointer Pointer
	*/
	void setNeighborP(geometry::Neighbor *d_neighbor_p);

	/*! Get the pointer to neighborhood object
	 * @return pointer
	 */
	geometry::Neighbor* getNeighborP();

	/*! Sets the volume correction pointer 
	* @param pointer Pointer
	*/
	void setVolumeCorrectionP(geometry::VolumeCorrection* pointer);


	/*! Get the pointer to volume correction
	 * @return pointer
	 */
	geometry::VolumeCorrection* getVolumeCorrectionP();

	/*! Sets the pointer to the state-based helper
	 * @param pointer Pointer
	 */
	void setStateBasedHelperFunctionsP(util::StateBasedHelperFunctions* pointer);

	/*! Get the pointer to state-nased helper functions
	 * @return pointer
	 */
	util::StateBasedHelperFunctions* getStateBasedHelperFunctionsP();

	/*! Sets the pointer to the displacement loading object
	 * @param pointer Pointer
	 */
	void setDisplacementLoadingP(loading::ULoading * pointer);


	/*! Get the pointer to displacement loading object
	 * @return pointer
	 */
	loading::ULoading* getDisplacementLoadingP();

	/*! Sets the pointer to the force loading object
	 * @param pointer Pointer
	 */
	void setForceLoadingP(loading::FLoading * pointer);

	/*! Get the pointer to force loading object
	 * @return pointer
	 */
	loading::FLoading* getForceLoadingP();

	/*! Sets the pointer to the extension state
	 * @param pointer Pointer
	 */
	void setExtensionP(std::vector<std::vector<double>>* pointer);


	/*! Get the pointer to extension state
	 * @return pointer
	 */
	std::vector<std::vector<double>>* getExtensionP();

	/*! Sets the pointer to the stress tensor
	 * @param pointer Pointer
	 */
	void setStressTensorP(std::vector<util::Matrix33>* pointer);


	/*! Get the pointer to stress tensor
	 * @return pointer
	 */
	std::vector<util::Matrix33>* getStressTensorP();

	/*! Sets the pointer to the strain tensor
	 * @param pointer Pointer
	 */
	void setStrainTensorP(std::vector<util::Matrix33>* pointer);

	/*! Get the pointer to strain tensor
	 * @return pointer
	 */
	std::vector<util::Matrix33>* getStrainTensorP();

	/*! Sets the pointer to the dilatation state
	 * @param pointer Pointer
	 */
	void setDilatationP(std::vector<double> * pointer);

	/*! Get the pointer to dilatation
	 * @return pointer
	 */
	std::vector<double> * getDilatationP();

	/*! Sets the pointer to the total reaction force 
	 * @param pointer Pointer
	 */
	void setTotalReactionForceP(std::vector<double>* pointer);


	/*! Get the pointer to total reaction force
	 * @return pointer
	 */ 
	std::vector<double>* getTotalReactionForceP();

	/*! Sets the pointer to the strain energy
	 * @param pointer Pointer
	 */
	void setStrainEnergyP(std::vector<float>* pointer);

	/*! Get the pointer to strain energy
	 * @return pointer
	 */
	std::vector<float>* getStrainEnergyP();

	/*! Sets the pointer to the work done vector
	 * @param pointer Pointer
	 */
	void setWorkDoneP(std::vector<float>* pointer);

	/*! Get the pointer to work done
	 * @return pointer
	 */
	std::vector<float>* getWorkDoneP();

	/*! Sets the pointer to the phi vector
	 * @param pointer Pointer
	 */
	void setPhiP(std::vector<float>* pointer);



	/*! Get the pointer to phi
 	 * @return pointer
	 */
	std::vector<float>* getPhiP();

	/*! Sets the pointer to damage function
	 * @param pointer Pointer
	 */
	void setDamageFunctionP(std::vector<float>* pointer);


	/*! Get the pointer to damage function
	 * @return pointer
	 */
	std::vector<float>* getDamageFunctionP();

	/*! Sets the pointer to the bond-based fracture energy
	 * @param pointer Pointer
	 */
	void setBBFractureEnergyP(std::vector<float>*pointer);


	/*! Get the pointer to bond-based fracture energy
	 * @return pointer
	 */
    std::vector<float>* getBBFractureEnergyP();

	/*! Sets the pointer to the fracture energy
	* @param pointer Pointer
	*/
	void setFractureEnergyP(std::vector<float>*pointer);

	/*! Get the pointer to fracture energy
	 * @return pointer
	 */
    std::vector<float>* getFractureEnergyP();

	/*! Sets the pointer to the kinetic energy
	* @param pointer Pointer
	*/
	void setKineticEnergyP(std::vector<float>* pointer);

	/*! Get the pointer to kinetic energy
	 * @return pointer
	 */
	std::vector<float>* getKineticEnergyP();

	/** @}*/

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
	std::vector<std::vector<double>>* d_extension_p = nullptr;

	/*! @brief Dilatation of nodes */
	std::vector<double> *d_dilatation_p = nullptr;

	/*! @brief Pointer to the strain energy vector */
	std::vector<float>* d_e_p = nullptr;

	/*! @brief Pointer to the total reaction force vector */
	std::vector<double>* d_total_reaction_force_p = nullptr;

	/*! @brief Pointer to the Work done on each of the nodes */
  	std::vector<float>* d_w_p = nullptr;

	/*! @brief Pointer to damage function \f$ \phi \f$ at the nodes */
  	std::vector<float>* d_phi_p = nullptr;

	/*! @brief Pointer to damage function \f$ Z \f$ at the nodes */
  	std::vector<float>* d_Z_p = nullptr;

	/*! @brief Pointer to bond-based fracture energy of the nodes */
  	std::vector<float>* d_eFB_p = nullptr;

	/*! @brief Pointer to fracture energy of the nodes */
  	std::vector<float>* d_eF_p = nullptr;

	/*! @brief Pointer to the kinetic energy of the nodes */
	std::vector<float>* d_ke_p = nullptr;

	/*! @brief Pointer to the strain tensor vector */
	std::vector<util::Matrix33>* d_strain_p = nullptr;

	/*! @brief Pointer to the stress tensor vector */
    std::vector<util::Matrix33>* d_stress_p = nullptr;

	/** @}*/

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

	/** @}*/

};

}

#endif /* SRC_DATA_DATAMANAGER_H_ */
