// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef MODEL_MODEL_H
#define MODEL_MODEL_H

#include <cstdlib>
#include <hpx/config.hpp>
#include <vector>

#include "../util/matrix.h" // definition of SymMatrix3
#include "../util/point.h"  // definition of Point3

/*!
 * @brief Collection of Peridynamic models
 *
 * This namespace provides collection of Peridynamic models. Depending on the
 * spatial discretization, e.g. finite difference, weak finite element, nodal
 * finite element, and truss finite element, we get different implementation
 * of model.
 *
 * E.g., in FDModel we implement finite difference discretization with
 * explicit time integration of Peridynamic equation.
 *
 * @sa FDModel
 */
namespace model {

/*!
 * @brief A base class for different models
 *
 * This class provides a base for specific models, such as FDModel.
 */
class Model {

public:
  /**
   * @name Common methods
   */
  /**@{*/

  /*!
   * @brief Return the current time step
   * @return Time step
   */
  virtual size_t currentStep();

  /*!
   * @brief Return the total energy
   * @return Total energy
   */
  virtual float getEnergy();

  /** @}*/

protected:
  /**
   * @name Major simulation data
   *
   * Major simulation data are those which play direct role in the simulation
   * and without declaring these simulation will not work. For Peridynamics,
   * displacement at current time, velocity at current time, and total force
   * are the major simulation data.
   *
   * For state-based peridynamics, we also need hydrostatic strain at current
   * time. If we choose we can remove this data from major simulation list,
   * however, this will increase the computational load.
   *
   * In addition to major simulation data listed here, there are data in
   * other classes, for example nodal data fe::Mesh::d_nodes in fe::Mesh,
   * which are also necessary for the simulation.
   *
   */
  /**@{*/

  /*! @brief Displacement of the nodes */
  std::vector<util::Point3> d_u;

  /*! @brief Velocity of the nodes */
  std::vector<util::Point3> d_v;

  /*! @brief Total force on the nodes */
  std::vector<util::Point3> d_f;

  /*! @brief Hydrostatic strains at the nodes */
  std::vector<double> d_hS;

  /*! @brief Current time step */
  size_t d_n;

  /** @}*/

  /**
   * @name Minor simulation data
   *
   * These datas are postprocessing data and have no role in simulation.
   * Since they do not play direct role in simulation, we can compromise in
   * their accuracy and use 'float' instead of 'double'.
   */
  /**@{*/

  /*! @brief Energy of the nodes */
  std::vector<float> d_e;

  /*! @brief Work done on each of the nodes */
  std::vector<float> d_w;

  /*! @brief Damage function \f$ \phi \f$ at the nodes */
  std::vector<float> d_phi;

  /*! @brief Damage function \f$ Z \f$ at the nodes */
  std::vector<float> d_Z;

  /*! @brief Fracture energy of the nodes */
  std::vector<float> d_eF;

  /*! @brief Bond-based fracture energy of the nodes */
  std::vector<float> d_eFB;

  /*! @brief Strains of the nodes */
  std::vector<util::SymMatrix3> d_strain;

  /*! @brief Stress of the nodes */
  std::vector<util::SymMatrix3> d_stress;

  /*! @brief Total internal energy */
  float d_te;

  /*! @brief Total work done */
  float d_tw;

  /*! @brief Total kinetic energy */
  float d_tk;

  /*! @brief Total fracture energy */
  float d_teF;

  /*! @brief Total bond-based fracture energy */
  float d_teFB;

  /** @}*/
};

} // namespace model

#endif // MODEL_MODEL_H