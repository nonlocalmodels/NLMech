// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef MODEL_FDMODEL_H
#define MODEL_FDMODEL_H

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
//class Quadrature;
} // namespace fe

namespace geometry {
class Fracture;
class InteriorFlags;
class Neighbor;
} // namespace geometry

namespace inp {
struct ModelDeck;
class Input;
class Policy;
} // namespace inp

namespace loading {
class InitialCondition;
class Loading;
} // namespace loading

namespace material {
class Material;
} // namespace material

namespace model {

/**
 * \defgroup Explicit Explicit
 */
/**@{*/

/*! @brief A class for \a finite \a difference \a approximation of
 * \b Peridynamics
 *
 * In this class we implement the \a finite \a difference \a approximation of
 * \b peridynamics.
 *
 * We consider \a explicit \a scheme such as \a central \a difference and
 * \a velocity \a verlet for time integration.
 *
 * This class acts as a link to lower rank classes, such as Mesh, Loading,
 * InitialCondition, Fracture, etc, and uses the methods and data of the
 * lower rank classes to run simulations.
 *
 * @note 1. We can run finite difference on any finite element mesh as long as
 * mesh consists of only one type of elements. Therefore, we are restricted
 * to run finite difference simulation only on uniform grids. User can
 * prepare a mesh using \b Gmsh and use its .msh file to run the finite
 * difference approximation.
 *
 * @note 2. Currently only dimension 2 is supported.
 *
 * @note 3. Either triangle or quadrangle element mesh are supported.
 */
class FDModel : public Model {

public:
  /*!
   * @brief Constructor
   * @param deck The input deck
   */
  explicit FDModel(inp::Input *deck);

  /**
   * @name Common methods
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
   * @name Methods to initialize the data
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
   * @name Methods to implement explicit time integration
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
   * @name Methods to apply boundary condition and initial condition
   */
  /**@{*/

  /*!
   * @brief Performs setup of boundary condition, such as filling node ids
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
   * @name Methods to handle output and debug
   */
  /**@{*/

  /*!
   * @brief Output the snapshot of data at current time step
   */
  void output();

  /*!
   * @brief Performs debug operations and outputs message to the screen
   * @param e_old at previous time step
   */
  void debug(float e_old);

  /** @}*/

private:

  /*! @brief Model deck */
  inp::ModelDeck *d_modelDeck_p;

  /**
   * @name Data: High level objects
   */
  /**@{*/

  /*!
   * @brief Pointer to Mass matrix object containing mass matrix (if any)
   *
   * @sa MassMatrix
   */
  fe::MassMatrix *d_massMatrix_p;

  /*!
   * @brief Pointer to Mesh object
   *
   * @sa Mesh
   */
  fe::Mesh *d_mesh_p;

  /*!
   * @brief Pointer to Quadrature object
   *
   * @sa Quadrature
   */
//  fe::Quadrature *d_quadrature_p;

  /*!
   * @brief Pointer to Fracture object
   *
   * @sa Fracture
   */
  geometry::Fracture *d_fracture_p;

  /*!
   * @brief Pointer to Neighbor object
   *
   * @sa Neighbor
   */
  geometry::Neighbor *d_neighbor_p;

  /*! @brief Pointer to InteriorFlags object
   *
   * @sa InteriorFlags
   */
  geometry::InteriorFlags *d_interiorFlags_p;

  /*! @brief Pointer to Input object
   *
   * @sa Input
   */
  inp::Input *d_input_p;

  /*! @brief Pointer to Policy object
   *
   * @sa Policy
   */
  inp::Policy *d_policy_p;

  /*! @brief Pointer to InitialCondition object
   *
   * @sa InitialCondition
   */
  loading::InitialCondition *d_initialCondition_p;

  /*! @brief Pointer to Loading object
   *
   * @sa Loading
   */
  loading::Loading *d_loading_p;

  /*! @brief Pointer to Material object
   *
   * @sa Material
   */
  material::Material *d_material_p;

  /** @}*/
};

/** @}*/

} // namespace model

#endif // MODEL_FDMODEL_H
