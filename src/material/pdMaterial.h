// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef MATERIAL_PD_MATERIAL_H
#define MATERIAL_PD_MATERIAL_H

#include "util/point.h"         // definition of Point3
#include <string>
#include <vector>

// forward declaration
namespace inp {
struct MaterialDeck;
}

namespace material {
namespace pd {
class BaseMaterial;
class BaseInfluenceFn;
} // namespace pd
} // namespace material

/*!
 * @brief Collection of methods and database related to material
 *
 * This namespace provides collection of different materials, such as
 * peridynamic material, etc.
 *
 * At present we have implemented both bond-based and state-based model which
 * can also. We consider \b RNP regularized potential as well \b PMB material.
 */
namespace material {

/*!
 * @brief Collection of methods and database related to peridynamic material
 *
 * At present we have implemented both bond-based and state-based model which
 * can also. We consider \b RNP regularized potential as well \b PMB material.
 *
 * @sa Material
 */

namespace pd {

/*!
 * @brief A class providing methods to compute energy density and force of
 * peridynamic material
 */
class Material {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   * @param dim Dimension
   * @param horizon Horizon
   */
  Material(inp::MaterialDeck *deck, const size_t &dim, const double &horizon);

  /*!
   * @brief Returns true if state-based potential is active
   * @return True/false
   */
  bool isStateActive();

  /*!
   * @brief Returns energy and force between bond
   *
   * @param r Reference (initial) bond length
   * @param s Bond strain
   * @param fs Bond fracture state
   * @param break_bonds Flag specifying if bond's fracture state is to be
   * modified (needed to implement \a no-fail region
   * @return Value Pair of energy and force
   */
  std::pair<double, double> getBondEF(const double &r, const double &s,
                                      bool &fs, const bool &break_bonds);

  /*!
   * @brief Returns hydrostatic energy and force
   *
   * @param theta Hydrostatic strain
   * @return Value Pair of energy and force
   */
  std::pair<double, double> getStateEF(const double &theta);

  /*!
   * @brief Returns hydrostatic force density
   *
   * @param g_prime Derivative of hydrostatic potential function
   * @param r Reference (initial) bond length
   * @return Value Force density
   */
  double getStateForce(const double &g_prime, const double &r);

  /*!
   * @brief Returns contribution of bond to hydrostatic strain
   *
   * @param S Bond strain
   * @param r Reference (initial) bond length
   * @return Value Contribution to hydrostatic strain
   */
  double getBondContribToHydroStrain(const double &S, const double &r);

  /*!
   * @brief Returns the value of influence function
   *
   * @param r Reference (initial) bond length
   * @return Value Influence function at r
   */
  double getInfFn(const double &r);

  /*!
   * @brief Returns the moment of influence function
   *
   * If \f$ J(r) \f$ is the influence function for \f$ r\in [0,1)\f$ then \f$
   * i^{th}\f$ moment is given by \f[ M_i = \int_0^1 J(r) r^i dr. \f]
   *
   * @param i ith moment
   * @return Value Moment
   */
  double getMoment(const size_t &i);

  /*!
   * @brief Returns the bond strain
   * @param dx Reference bond vector
   * @param du Difference of displacement
   * @return Value Bond strain \f$ S = \frac{du \cdot dx}{|dx|^2} \f$
   */
  double getS(const util::Point3 &dx, const util::Point3 &du);

  /*!
   * @brief Returns the density of the material
   * @return density Density of the material
   */
  double getDensity();

  /*!
   * @brief Returns true if bond contributes to hydrostatic force
   * @param S Bond strain
   * @param r Reference bond length
   * @return True/false
   */
  bool addBondContribToState(const double &S, const double &r);

  /*!
   * @brief Returns critical bond strain
   *
   * @param r Reference length of bond
   * @return Value Critical strain
   */
  double getSc(const double &r);

  /*!
   * @brief Returns the material deck
   * @return deck Material deck
   */
  inp::MaterialDeck * getMaterialDeck();

private:
  /**
   * @name Internal data
   */
  /**@{*/

  /*! @brief Flag indicating if peridynamic state potential is active */
  bool d_stateActive;

  /*! @brief Horizon */
  double d_horizon;

  /*! @brief Dimension */
  size_t d_dimension;

  /*! @brief Density */
  double d_density;

  /** @}*/

  /**
   * @name Base objects implementing particular material model
   */
  /**@{*/

  /*! @brief Base object for pd force and energy */
  material::pd::BaseMaterial *d_baseMaterial_p;

  /*! @brief Base object for influence function */
  material::pd::BaseInfluenceFn *d_baseInfluenceFn_p;

  /** @}*/

  /*! @brief Input material deck */
  inp::MaterialDeck *d_deck_p;
};

} // namespace pd

} // namespace material

#endif // MATERIAL_PD_MATERIAL_H
