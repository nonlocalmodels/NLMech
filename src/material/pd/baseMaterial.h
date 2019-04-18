// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef MATERIAL_PD_BASEMATERIAL_H
#define MATERIAL_PD_BASEMATERIAL_H

#include <cstring>
#include <vector>

namespace material {

namespace pd {

/*! @brief A base class providing methods to compute energy density and force
 *
 * This is a base class which provides method to compute pairwise energy and
 * force, as well as hydrostatic force and energy. Later on we can introduce
 * a damage model as another abstraction.
 */
class BaseMaterial {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   */
  BaseMaterial(const size_t &dim, const double &horizon)
      : d_horizon(horizon), d_dimension(dim){};

  /*!
   * @brief Returns energy and force between bond
   *
   * @param r Reference (initial) bond length
   * @param s Bond strain
   * @param influence Value of influence function at r
   * @param fs Bond fracture state
   * @return Value Pair of energy and force
   */
  virtual std::pair<double, double> getBondEF(const double &r, const double &s,
                                              const double &influence,
                                              bool &fs) {
    return std::make_pair(0., 0.);
  };

  /*!
   * @brief Returns energy and force between bond for \a no-fail region bonds
   *
   * @param r Reference (initial) bond length
   * @param s Bond strain
   * @param influence Value of influence function at r
   * @return Value Pair of energy and force
   */
  virtual std::pair<double, double>
  getBondEFNoFail(const double &r, const double &s, const double &influence) {
    return std::make_pair(0., 0.);
  };

  /*!
   * @brief Returns hydrostatic energy and force
   *
   * @param theta Hydrostatic strain
   * @return Value Pair of energy and force
   */
  virtual std::pair<double, double> getStateEF(const double &theta) {
    return std::make_pair(0., 0.);
  };

  /*!
   * @brief Returns hydrostatic force density
   *
   * @param g_prime Derivative of hydrostatic potential function
   * @param r Reference (initial) bond length
   * @return Value Force density
   */
  virtual double getStateForce(const double &g_prime, const double &r) {
    return 0.;
  };

  /*!
   * @brief Returns contribution of bond to hydrostatic strain
   *
   * @param S Bond strain
   * @param r Reference (initial) bond length
   * @return Value Contribution to hydrostatic strain
   */
  virtual double getBondContribToHydroStrain(const double &S, const double &r) {
    return 0.;
  };

  /*!
   * @brief Returns critical bond strain
   *
   * @param r Reference length of bond
   * @return Value Critical strain
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
