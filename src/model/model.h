// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef MODEL_ABSTRACTIONS_H
#define MODEL_ABSTRACTIONS_H

#include <cstdlib>
#include <hpx/config.hpp>
#include <vector>

#include "../util/matrix.h"
#include "../util/point.h"

// namespace util {
//
// forward declaration
// struct Point3;
// struct SymMatrix3;
//}

//! Collection of models
namespace model {

/*! @brief Abstraction of a model */
class Model {

public:
  //  /*!
  //   * @brief Constructor
  //   * @param deck The input deck
  //   */
  //  Model();

  /**
   * \defgroup Common methods
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
   * \defgroup Major simulation data
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
   * \defgroup Minor simulation data
   * These data are postprocessing data and have no role in simulation.
   * Use "float" to store values.
   */
  /**@{*/

  /*! @brief Energy of the nodes */
  std::vector<float> d_e;

  /*! @brief Work done on each of the nodes */
  std::vector<float> d_w;

  /*! @brief Damage function \phi at the nodes */
  std::vector<float> d_phi;

  /*! @brief Damage function Z at the nodes */
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

#endif