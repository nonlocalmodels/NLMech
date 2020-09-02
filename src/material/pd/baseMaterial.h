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

#include "util/point.h"

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
      : d_horizon(horizon), d_dimension(dim), d_stateActive(false), d_name(""){};

  /*!
   * @brief Returns true if state-based potential is active
   * @return bool True/false
   */
  bool isStateActive() const { return d_stateActive; };

  /*!
   * @brief Returns name of the material
   * @return string Name
   */
  std::string name() const {return d_name;};

  /*!
   * @brief Returns energy and force between bond
   *
   * @param i Id of node 1
   * @param j Local id in the neighborlist of node i
   * @return pair Pair of energy and force
   */
  virtual std::pair<util::Point3,double> getBondEF(size_t i , size_t j){

	  return std::make_pair(util::Point3(), 0.);

  };

  /*!
   * @brief Returns the bond strain
   * @param dx Reference bond vector
   * @param du Difference of displacement
   * @return strain Bond strain
   */
  virtual double getS(const util::Point3 &dx, const util::Point3 &du) {
    return
        0.;
  };

  /*!
   * @brief Returns the bond strain
   * @param i Id of node 1
   * @param j Id of node 2
   * @return strain Bond strain
   */
  virtual double getS(size_t i, size_t j) { return 0.; };

  /*! @brief Returns critical bond strain
   *
   * @param r Reference length of bond
   * @return strain Critical strain
   */
  virtual double getSc(const double &r) { return 0.; };

  /*!
   * @brief Returns critical bond strain
   *
   * @param i Id of node 1
   * @param j Id of node 2
   * @return strain Critical strain
   */
  virtual double getSc(size_t i , size_t j) { return 0.; };

  /*!
   * @brief Get direction of bond force
   * @param dx Relative bond vector (reference configuration)
   * @param du Relative bond displacement vector
   * @return vector Unit vector along the bond force
   */
  virtual util::Point3 getBondForceDirection(const util::Point3 &dx,
                                     const util::Point3 &du) const {return
                                     util::Point3(); }

  /*!
   * @brief Returns strain tensor
   *
   * @param i Id of node
   * @return strain tensor
   */
  virtual util::Matrix33 getStrain(size_t i) {

    return util::Matrix33();
  }

  /*!
   * @brief Returns stress tensor
   *
   * @param i Id of node
   * @return stress tensor
   */
  virtual util::Matrix33 getStress(size_t i) {

    return util::Matrix33();
  }

 /*!
   * @brief Let the material class in the quasi-static case know that there is a new 
   * loading step
   */
  virtual void update(){}

  /*!
   * @brief Returns the value of influence function
   *
   * @param r Reference (initial) bond length
   * @return value Influence function at r
   */
  virtual double getInfFn(const double &r) const {return 0.; };

  /*!
   * @brief Returns horizon
   *
   * @return horizon Horizon
   */
  double getHorizon() const {
    return d_horizon;
  }

  /*!
   * @brief Returns the density of the material
   * @return density Density of the material
   */
  double getDensity() const {return d_density; };

  /*!
* @brief Returns the string containing information about the instance of
* the object
*
* @param nt Number of tabs to append before each line of string
* @param lvl Level of information sought (higher level means more
* information)
* @return string String containing information about this object
* */
  virtual std::string printStr(int nt = 0, int lvl = 0) const {
    // TODO implement
    return "";
  };

  /*!
   * @brief Prints the information about the instance of the object
   *
   * @param nt Number of tabs to append before each line of string
   * @param lvl Level of information sought (higher level means more
   * information)
   * */
  virtual void print(int nt = 0, int lvl = 0) const { std::cout << printStr(nt,
      lvl); };

protected:

  /*! @brief Horizon */
  double d_horizon;

  /*! @brief Dimension */
  size_t d_dimension;

  /*! @brief Flag indicating if peridynamic state potential is active */
  bool d_stateActive;

  /*! @brief Density */
  double d_density;

  /*! @brief Name of material */
  std::string d_name;
};

} // namespace pd

} // namespace material

#endif // MATERIAL_PD_BASEMATERIAL_H
