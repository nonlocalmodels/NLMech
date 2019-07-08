// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef MODEL_FDMODEL_H
#define MODEL_FDMODEL_H

#include <hpx/config.hpp>
#include "../model.h"
#include <vector>

// forward declaration of class
namespace fe {
class MassMatrix;
class Mesh;
// class Quadrature;
} // namespace fe

namespace geometry {
class Fracture;
class InteriorFlags;
class Neighbor;
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
class Material;
}
} // namespace material

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
  /*!
   * @brief Main driver to simulate
   * @param deck Input deck
   */
  void run(inp::Input *deck);

  /*!
   * @brief Restarts the simulation from previous state
   * @param deck Input deck
   */
  void restart(inp::Input *deck);

  /*!
   * @brief Computes peridynamic forces
   */
  void computeForces();

  /*!
   * @brief Computes hydrostatic strains for force calculation
   */
  void computeHydrostaticStrains();

  /*!
   * @brief Computes postprocessing quantities
   */
  void computePostProcFields();

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
   *
   * Central difference scheme
   * \f[ u_{new} = \Delta t^2 (f_{int} + f_{ext}) / \rho  +
   * \Delta t v_{old} + u_{old} \f]
   * \f[ v_{new} = \frac{u_{new} - u_{old}}{\Delta t}. \f]
   */
  void integrateCD();

  /*!
   * @brief Perform time integration using velocity-verlet scheme
   *
   * Velocity verlet scheme
   * 1. \f$ v_{mid} = v_{old} + \frac{\Delta t}{2} (f_{int,old} + f_{ext,
   * old}) / \rho
   * \f$
   *
   * 2. \f$ u_{new} = u_{old} + \Delta t * v_{mid} \f$
   *
   * 3. \f$ v_{new} = v_{mid} +  \frac{\Delta t}{2} (f_{int,new} + f_{ext,
   * new}) / \rho \f$
   */
  void integrateVerlet();

  /** @}*/

  /**
   * @name Methods to handle output and debug
   */
  /**@{*/

  /*!
   * @brief Output the snapshot of data at current time step
   *
   * Name of data in the output simulation file and the name of data they
   * refer to
   *
   * Major simulation data
   * | Name in output file          | Data name             |
   * | ---------------------------- | --------------------- |
   * | Displacement                 | model::Model::d_u     |
   * | Velocity                     | model::Model::d_v     |
   * | Force                        | model::Model::d_f * nodal volume |
   * | time                         | model::Model::d_time  |
   *
   * Minor simulation data (if these are populated and computed). Calculation
   * and population of post-processing data can be controlled by specifying
   * memory control level in inp::PolicyDeck::d_memControlFlag or by
   * enabling/disabling post-processing calculation in
   * inp::PolicyDeck::d_enablePostProcessing.
   * | Name in output file          | Data name             |
   * | ---------------------------- | --------------------- |
   * | Force_Density                | model::Model::d_f     |
   * | Strain_Energy                | model::Model::d_e    |
   * | Work_Done                    | model::Model::d_w     |
   * | Fixity                       | fe::Mesh::getFixityP() |
   * | Node_Volume                  | fe::Mesh::getNodalVolumeP() |
   * | Damage_Phi                   | model::Model::d_phi    |
   * | Damage_Z                     | model::Model::d_Z      |
   * | Fracture_Perienergy_Bond     | model::Model::d_eFB    |
   * | Fracture_Perienergy_Total    | model::Model::d_eF     |
   * | Total_Energy                 | model::Model::d_Z      |
   * | Total_Fracture_Perienergy_Bond   | model::Model::d_teFB   |
   * | Total_Fracture_Perienergy_Total  | model::Model::d_teF    |
   *
   */
  void output();

  /*! @brief Checks if output frequency is to be modified
   *
   * 1. If valid criteria is specified, this method modifes the current output
   * interval to new output interval.
   *
   * 2. This is useful when we require snapshot of simulation at small intervals
   * only when results are interesting. For example, for crack propagation
   * problem, we need more snapshots of simulation during time when the crack
   * is propagating. We do not need many snapshots of time when crack is not
   * growing or moving.
   *
   * 3. Currently we have implemented `max_Z` criteria which checks the maximum
   * value of damage at nodes. If maximum value exceeds the user specified
   * value, we change the output interval to small value as specified by user
   * . Until maximum of damage touches value given by user, the code produces
   * output at larger interval.
   *
   * 4. We next discuss the naming of output files. We will take the
   * smaller output interval as a basis for naming the output files. Let us
   * look at simple example:
   *
   *    - Suppose total number of time steps `N=1500`,
   *    - `dt_out_large = 150` when the criteria `max_Z` is not met
   *    - `dt_out_small = 50` when criteria is met
   *    - Suppose `max_Z` criteria is not met till `1200`th step, then upto
   *    `1200`th step, output files will be:
   *        * `output_0.vtu, output_3.vtu, output_6.vtu, ... output_21.vtu,
   *        output_24.vtu`
   *        * Files above correspond to time steps
   *    `0, 3*dt_out_small = 150, 6*dt_out_small = 300, ..., 1050, 1200`
   *    - From `1200` step onwards output files will be:
   *        * `output_25.vtu, output_26.vtu, output_27.vtu, ..., output_30.vtu`
   *        * Files above correspond to time steps
   *    `25*dt_out_small = 1250, 26*dt_out_small = 1300, ..., 1500`
   *
   *    - So in brief, we consider `dt_out_small` to get the tag for output file
   * . Output intervals however depend on the criteria, and can be at
   * `dt_out_large`, or it can be at the same interval as `dt_out_small`.
   *
   */
  void checkOutputCriteria();

  /** @}*/

private:
  /*! @brief Model deck */
  inp::ModelDeck *d_modelDeck_p;

  /*! @brief Restart deck */
  inp::RestartDeck *d_restartDeck_p;

  /*! @brief Output deck */
  inp::OutputDeck *d_outputDeck_p;

  /**
   * @name Data: High level objects
   */
  /**@{*/

  /*! @brief Pointer to Mass matrix object containing mass matrix (if any) */
  fe::MassMatrix *d_massMatrix_p;

  /*! @brief Pointer to Mesh object */
  fe::Mesh *d_mesh_p;

  /*! @brief Pointer to Fracture object */
  geometry::Fracture *d_fracture_p;

  /*! @brief Pointer to Neighbor object */
  geometry::Neighbor *d_neighbor_p;

  /*! @brief Pointer to InteriorFlags object */
  geometry::InteriorFlags *d_interiorFlags_p;

  /*! @brief Pointer to Input object */
  inp::Input *d_input_p;

  /*! @brief Pointer to Policy object */
  inp::Policy *d_policy_p;

  /*! @brief Pointer to InitialCondition object */
  loading::InitialCondition *d_initialCondition_p;

  /*! @brief Pointer to displacement Loading object */
  loading::ULoading *d_uLoading_p;

  /*! @brief Pointer to force Loading object */
  loading::FLoading *d_fLoading_p;

  /*! @brief Pointer to Material object */
  material::pd::Material *d_material_p;

  /** @}*/
};

/** @}*/

} // namespace model

#endif // MODEL_FDMODEL_H
