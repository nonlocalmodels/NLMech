// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef FDMODEL_H
#define FDMODEL_H

#include <hpx/config.hpp>
#include <vector>

// include model abstraction class
#ifndef MODEL_ABSTRACTIONS_H
#include "../model.h"
#endif

// forward declaration of class
namespace fe {
class MassMatrix;
class Mesh;
class Quadrature;
} // namespace fe

namespace geometry {
class Fracture;
class InteriorFlags;
class Neighbor;
} // namespace geometry

namespace inp {
class Input;
class Policy;
} // namespace inp

namespace loading {
class InitialCondition;
class Loading;
} // namespace loading

namespace material {
class Material;
}

namespace model {

/*! @brief Implements finite difference discretization of Peridynamics */
class FDModel : public Model {

public:
  /*!
   * @brief Constructor
   * @param deck The input deck
   */
  explicit FDModel(inp::Input *deck);

  /**
   * \defgroup Method providing view of the state of the Model
   */
  /**@{*/

  /*!
   * @brief Return the current time step
   * @return Time step
   */
  size_t currentStep() override;

  /*!
   * @brief Return the total energy
   * @return Total energy
   */
  float getEnergy() override;

  /** @}*/

private:
  /**
   * \defgroup Methods to initialize the data
   */
  /**@{*/

  /*!
   * @brief Initialize high level data members
   */
  void initHObjects();

  /*!
   * @brief Initialize remaining data members
   */
  void init();

  /** @}*/

  /**
   * \defgroup Methods to implement explicit time integration
   */
  /**@{*/

  /*!
   * @brief Perform time integration
   */
  void integrate();

  /*!
   * @brief Perform time integration using central-difference scheme
   */
  void integrateCD();

  /*!
   * @brief Perform time integration using velocity-verlet scheme
   */
  void integrateVerlet();

  /** @}*/

  /**
   * \defgroup Methods to apply boundary condition and initial condition
   */
  /**@{*/

  /*! @brief Performs setup of boundary condition, such as filling node ids
   * in loading object and setting of fixity of nodes
   */
  void setupBoundaryCondition();

  /*!
   * @brief Apply displacement boundary condition to current position
   */
  void applyDisplacementBC();

  /*!
   * @brief Apply external loading to the nodes
   */
  void applyForceBC();

  /*! @brief Writes the initial condition to the current position and current
   * velocity
   */
  void applyInitialCondition();

  /** @}*/

  /**
   * \defgroup Methods to handle output and debug
   */
  /**@{*/

  /*!
   * @brief Output the snapshot of data at current time step
   */
  void output();

  /*!
   * @brief Performs debug operations and outputs message to the screen
   * @param Energy at previous time step
   */
  void debug(float e_old);

  /** @}*/

private:
  /**
   * \defgroup High level objects as member data
   */
  /**@{*/

  /*! @brief Pointer to Mass matrix object containing mass matrix (if any) */
  fe::MassMatrix *d_massMatrix_p;

  /*! @brief Pointer to Mesh object containing list of node,
   * element-node connectivity and other fem related information
   */
  fe::Mesh *d_mesh_p;

  /*! @brief Pointer to Quadrature object providing methods for quadrature
   * point approximation of integration */
  fe::Quadrature *d_quadrature_p;

  /*! @brief Pointer to Fracture object containing fracture state
   * information of each interacting bond
   */
  geometry::Fracture *d_fracture_p;

  /*! @brief Pointer to Neighbor object containing list of neighboring nodes
   * for each node
   */
  geometry::Neighbor *d_neighbor_p;

  /*! @brief Pointer to InteriorFlags object containing flags indicating
   * whether the nodes are in the interior or near the boundary
   */
  geometry::InteriorFlags *d_interiorFlags_p;

  /*! @brief Pointer to Input object containing all data specified in input
   * files
   */
  inp::Input *d_input_p;

  /*! @brief Pointer to Policy object which implements the policies in the
   * code. E.g. whether to allocate space for minor simulation data, whether
   * to allocate space for interior flags, etc
   */
  inp::Policy *d_policy_p;

  /*! @brief Pointer to InitialCondition object which handles application of
   * initial condition to velocity and displacement
   */
  loading::InitialCondition *d_initialCondition_p;

  /*! @brief Pointer to Loading object which handles application of
   * displacement and force boundary conditions
   */
  loading::Loading *d_loading_p;

  /*! @brief Pointer to Material object which provides the method to compute
   * forces between interacting bonds (bond-based peridynamics), and forces due
   * to the hydrostatic strain (state-based peridynamic)
   */
  material::Material *d_material_p;

  /** @}*/
};

} // namespace model

#endif // FDMODEL_H
