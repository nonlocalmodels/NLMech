////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef MATERIAL_PD_MATERIAL_H
#define MATERIAL_PD_MATERIAL_H

#include "util/point.h" // definition of Point3
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
 * This namespace provides collection of different materials. Currently we
 * only have peridynamic material models.
 */
namespace material {

/*!
 * @brief Collection of methods and database related to peridynamic material
 *
 * At present we have implemented both bond-based and state-based model. We
 * consider \b RNP regularized potential proposed and studied in [Lipton
 * 2016](https://link.springer.com/article/10.1007/s10659-015-9564-z),
 * [Jha and Lipton 2018](https://doi.org/10.1137/17M1112236), [Diehl et al
 * 2018](https://arxiv.org/abs/1806.06917),
 * [Jha and Lipton 2019](https://doi.org/10.1016/j.cma.2019.03.024). We have
 * also implemented PMB material model (Prototypical micro-elastic brittle
 * material), see [Silling 2000](https://www.sciencedirect
 * .com/science/article/pii/S0022509699000290).
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
   * @return bool True/false
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
   * @return pair Pair of energy and force
   */
  //std::pair<double, double> getBondEF(const double &r, const double &s,
   //                                   bool &fs, const bool &break_bonds);

  std::pair<util::Point3, double> getBondEF(size_t i , size_t j);

  util::Matrix33 getStrain(size_t i);

  util::Matrix33 getStress(size_t i);

  /*!
   * @brief Returns hydrostatic energy density
   *
   * @param theta Hydrostatic strain
   * @return energy Energy density
   */
  double getStateEnergy(const double &theta);

  /*!
   * @brief Returns hydrostatic force density
   *
   * @param theta Hydrostatic strain
   * @param r Reference (initial) bond length
   * @return force Force density
   */
  double getStateForce(const double &theta, const double &r);

  /*!
   * @brief Returns true if bond contributes to hydrostatic force
   * @param S Bond strain
   * @param r Reference bond length
   * @return bool True/false
   */
  bool doesBondContribToState(const double &S, const double &r);

  /*!
   * @brief Returns contribution of bond to hydrostatic strain
   *
   * @param S Bond strain
   * @param r Reference (initial) bond length
   * @return strain Contribution to hydrostatic strain
   */
  double getBondContribToHydroStrain(const double &S, const double &r);

  /*!
   * @brief Returns the bond strain
   * @param dx Reference bond vector
   * @param du Difference of displacement
   * @return strain Bond strain \f$ S = \frac{du \cdot dx}{|dx|^2} \f$
   */
  double getS(const util::Point3 &dx, const util::Point3 &du);

  /*!
   * @brief Returns critical bond strain
   *
   * @param r Reference length of bond
   * @return strain Critical strain
   */
  double getSc(const double &r);

  /*!
   * @brief Returns the density of the material
   * @return density Density of the material
   */
  double getDensity();

  /*!
   * @brief Returns the value of influence function
   *
   * @param r Reference (initial) bond length
   * @return value Influence function at r
   */
  double getInfFn(const double &r);

  /*!
   * @brief Returns the moment of influence function
   *
   * If \f$ J(r) \f$ is the influence function for \f$ r\in [0,1)\f$ then \f$
   * i^{th}\f$ moment is given by \f[ M_i = \int_0^1 J(r) r^i dr. \f]
   *
   * @param i ith moment
   * @return value Moment
   */
  double getMoment(const size_t &i);

  /*!
   * @brief Returns the material deck
   * @return deck Material deck
   */
  inp::MaterialDeck *getMaterialDeck();

private:

  /*! @brief Flag indicating if peridynamic state potential is active */
  bool d_stateActive;

  /*! @brief Horizon */
  double d_horizon;

  /*! @brief Dimension */
  size_t d_dimension;

  /*! @brief Density */
  double d_density;

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
