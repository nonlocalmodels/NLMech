////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef MATERIAL_PD_BASEMATERIAL_H
#define MATERIAL_PD_BASEMATERIAL_H

#include <cstring>
#include <vector>

namespace material {

namespace pd {

/*! @brief A base class providing methods to compute energy density and force
 *
 * This is a base class which provides method to compute pairwise energy and
 * force, as well as hydrostatic force and energy.
 */
class BaseMaterial {

public:
  /*!
   * @brief Constructor
   * @param dim Dimension
   * @param horizon Horizon
   */
  BaseMaterial(const size_t &dim, const double &horizon)
      : d_horizon(horizon), d_dimension(dim){};

  /*!
   * @brief Returns energy and force between bond
   *
   * @param r Reference (initial) bond length
   * @param s Bond strain
   * @param J Influence function at r
   * @param fs Bond fracture state
   * @return pair Pair of energy and force
   */
  virtual std::pair<double, double> getBondEF(const double &r, const double &s,
                                              const double &J,
                                              bool &fs) {
    return std::make_pair(0., 0.);
  };

  /*!
   * @brief Returns energy and force between bond for \a no-fail region bonds
   *
   * @param r Reference (initial) bond length
   * @param s Bond strain
   * @param J Influence function at r
   * @return pair Pair of energy and force
   */
  virtual std::pair<double, double>
  getBondEFNoFail(const double &r, const double &s, const double &J) {
    return std::make_pair(0., 0.);
  };

  /*!
   * @brief Returns hydrostatic energy density
   *
   * @param theta Hydrostatic strain
   * @return energy Energy density
   */
  virtual double getStateEnergy(const double &theta) {
    return 0.;
  };

  /*!
   * @brief Returns hydrostatic force density
   *
   * @param theta Hydrostatic strain
   * @param J Influence function at r
   * @return force Force density
   */
  virtual double getStateForce(const double &theta, const double &J) {
    return 0.;
  };

  /*!
   * @brief Returns true if bond contributes to hydrostatic force
   * @param S Bond strain
   * @param r Reference bond length
   * @return bool True/false
   */
  virtual bool doesBondContribToState(const double &S, const double &r) {
    return false;
  };

  /*!
   * @brief Returns contribution of bond to hydrostatic strain
   *
   * @param S Bond strain
   * @param r Reference (initial) bond length
   * @param J Influence function at r
   * @return strain Contribution to hydrostatic strain
   */
  virtual double getBondContribToHydroStrain(const double &S, const double &r, const
  double &J) {
    return 0.;
  };

  /*!
   * @brief Returns critical bond strain
   *
   * @param r Reference length of bond
   * @return strain Critical strain
   */
  virtual double getSc(const double &r) { return 0.; };

protected:
  /*! @brief Horizon */
  double d_horizon;

  /*! @brief Dimension */
  size_t d_dimension;
};

} // namespace pd

} // namespace material

#endif // MATERIAL_PD_BASEMATERIAL_H
