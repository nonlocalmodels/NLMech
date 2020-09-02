////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef MATERIAL_PD_BONDMATERIAL_NEW_H
#define MATERIAL_PD_BONDMATERIAL_NEW_H

#include "baseMaterial.h"
#include "util/point.h"
#include "util/matrixBlaze.h"
#include "geometry/neighbor.h"
#include "fe/mesh.h"

#include <vector>

#include <hpx/include/parallel_algorithm.hpp>

// forward declaration
namespace inp {
struct MaterialDeck;
struct OutputDeck;
}

// forward declaration
namespace data {
class DataManager;
}

// forward declaration of class
namespace fe {
class Mesh;
} // namespace fe

namespace material {
namespace pd {
class BaseInfluenceFn;
} // namespace pd
} // namespace material

namespace material {

namespace pd {

/*! @brief A Class implementing regularized nonlinear peridynamic model
 *
 * Provides method to compute energy and force using nonlinear bond-based
 * model introduced and studied in
 * [Lipton 2016](https://doi.org/10.1007/s10659-015-9564-z),
 * [Jha and Lipton 2018](https://epubs.siam.org/doi/10.1137/17M1112236), and
 * [Lipton and Jha 2019](https://doi.org/10.1007/s42102-019-00010-0).
 *
 * 1. The pairwise energy is
 * \f[ E_{bond} = \int_D \frac{1}{|B_\epsilon(x)|} \int_{B_\epsilon(x)} |y-x|
 * W_{bond}(S(y,x; u)) dy dx \f]
 * where \f$ B_\epsilon(x) \f$ is a ball (circle in 2-d) centered at \f$ x\f$
 * of  radius \f$ \epsilon\f$, \f$ |B_\epsilon(x)| \f$ is the volume (area in
 * 2-d) of ball, \f$ S(y,x;u) = \frac{u(y) - u(x)}{|y-x|} \cdot
 * \frac{y-x}{|y-x|} \f$
 * is the linearized bond strain (assuming small deformation).
 *
 * 2. Bond energy density \f$ W_{bond} \f$ is given by
 * \f[ W_{bond} (S(y,x;u)) = J^\epsilon(|y-x|) \frac{1}{\epsilon |y-x|} \psi
 * (|y-x| S(y,x;u)^2) \f]
 * where \f$ \psi(r) : \textsf{R}^+ \to \textsf{R} \f$ is positive, smooth,
 * concave
 * function with following properties
 * \f[ \lim_{r\to 0^+} \frac{\psi(r)}{r} = \psi'(0), \quad \lim_{r\to \infty}
 * \psi(r) = \psi_\infty < \infty. \f]
 * \f$ J^\epsilon(|y-x|)\f$ is the influence function.
 *
 * 3. Force at material point \f$ x \f$ is given by
 * \f[ f_{bond}(x) = \frac{4}{|B_\epsilon(x)|} \int_{B_\epsilon(x)}
 * \frac{J^\epsilon(|y-x|)}{\epsilon} \psi'(|y-x| S(y,x;u)^2) S(y,x;u)
 * \frac{y-x}{|y-x|} dy.\f]
 *
 * 4. In this class, we will assume
 *  \f[ \psi(r) = C ( 1-\exp[-\beta r] \f]
 * where \f$ C, \beta \f$ are the peridynamic material parameter
 * determined from the elastic and fracture properties of the material.
 */
class RNPBond : public BaseMaterial {

public:
  /*!
	 * @brief Constructor
	 * @param deck Pointer to the input deck
	 * @param DataManager Pointer to the data manager object
	 */
  RNPBond(inp::MaterialDeck *deck,data::DataManager* dataManager);

  /*!
	 * @brief Returns energy and force state between node i and node j
	 * @param i Id of node i
	 * @param j node of j
	 * @return Value Pair of energy and force
	 */
  std::pair<util::Point3, double> getBondEF(size_t i, size_t j);

  /*!
   * @brief Returns the bond strain
   * @param dx Reference bond vector
   * @param du Difference of displacement
   * @return strain Bond strain
   */
  double getS(const util::Point3 &dx, const util::Point3 &du);

  /*!
   * @brief Returns the bond strain
   * @param i Id of node 1
   * @param j Id of node 2
   * @return strain Bond strain
   */
  double getS(size_t i, size_t j);

  /*! @brief Returns critical bond strain
   *
   * @param r Reference length of bond
   * @return strain Critical strain
   */
  double getSc(const double &r);

  /*!
   * @brief Returns critical bond strain
   *
   * @param i Id of node 1
   * @param j Id of node 2
   * @return strain Critical strain
   */
  double getSc(size_t i , size_t j);

  /*!
   * @brief Get direction of bond force
   * @param dx Relative bond vector (reference configuration)
   * @param du Relative bond displacement vector
   * @return vector Unit vector along the bond force
   */
  util::Point3 getBondForceDirection(const util::Point3 &dx,
                                     const util::Point3 &du) const;

  /*!
   * @brief Returns the value of influence function
   *
   * @param r Reference (initial) bond length
   * @return value Influence function at r
   */
  double getInfFn(const double &r) const ;

private:
  /*!
   * @brief Computes rnp material parameters from elastic constants
   *
   * Either Young's modulus E or bulk modulus K, and either critical energy
   * release rate Gc or critical stress intensity factor KIc are needed.
   * Assuming Poisson's ratio \f$ \nu = \frac{1}{4} \f$ we compute lame
   * parameters \f$ \lambda, \mu \f$ where \f$ \mu = \lambda \f$ for bond-based.
   *
   * With lame parameters, we use following formula, see Equation (5.7) & (5
   * .8) of [Lipton 2016](https://doi.org/10.1007/s10659-015-9564-z)
   * - if \f$ d=2 \f$
   * \f[ \lambda = \mu = \frac{\psi'(0)}{4} M_2, \qquad Gc = \frac{4
   * \psi_{\infty}}{\pi} M_2. \f]
   * - if \f$ d=3\f$
   * \f[ \lambda = \mu = \frac{\psi'(0)}{5} M_3, \qquad Gc = \frac{3
   * \psi_{\infty}}{2} M_3. \f]
   *
   * Where \f$M_2, M_3\f$ are defined by
   * \f[ M_2 =\int_0^1 r^2 J(r) dr, \qquad M_3 = \int_0^1 r^3 J(r) dr. \f]
   *
   * For potential function \f$ \psi(r) = c ( 1-\exp[-\beta r])\f$, we have \f$
   * \psi'(0) = c\beta, \psi_{\infty} = c \f$. Thus, the values of \f$ c, \beta \f$
   * are given by
   * - if \f$ d=2 \f$
   * \f[ c = \frac{\pi G_c}{4} \frac{1}{M_2}, \qquad
   *     \beta = \frac{4 \lambda}{c} \frac{1}{M_2} .\f]
   * - if \f$ d=3 \f$
   * \f[ c = \frac{2 G_c}{3} \frac{1}{M_3}, \qquad
   *     \beta = \frac{5 \lambda}{c} \frac{1}{M_3} .\f]
   * @param deck Input material deck
   * @param M Moment of influence function
   */
  void computeParameters(inp::MaterialDeck *deck, const double &M);

  /*!
   * @brief Computes elastic and fracture properties from the rnp material
   * parameters
   *
   * This function does opposite of pd::Material::RNPBond::computeParameters.
   * From peridynamic material properties, it uses the relation between lame
   * parameters and peridynamic parameters to compute the lame parameters, and
   * from lame parameters it computes the elastic constants.
   *
   * @param deck Input material deck
   * @param M Moment of influence function
   */
  void computeMaterialProperties(inp::MaterialDeck *deck, const double &M);

  /**
   * @name Material parameters
   */
  /**@{*/

  /*! @brief Parameter C */
  double d_C;

  /*! @brief Parameter \f$ \beta \f$ */
  double d_beta;

  /** @}*/

  /*! @brief Inflection point of nonlinear function = \f$ 1/\sqrt{2\beta}\f$ */
  double d_rbar;

  /*! @brief Inverse of factor = \f$ \epsilon |B_\epsilon(0)|\f$ */
  double d_invFactor;

  /*! @brief Factor to multiply to critical strain to check if bond is
   * fractured
   *
   * For nonlinear model, we consider bond is broken when it
   * exceeds 10 times of critical strain. Typical value of factor is 10.
   */
  double d_factorSc;

  /*! @brief Flag which indicates if the breaking of bond is irreversible */
  bool d_irrevBondBreak;

  /*! @brief Kn coefficient for normal contact force between broken bonds */
  double d_contact_Kn = 0.;

  /*! @brief Contact radius for normal contact force between broken bonds */
  double d_contact_Rc = 0.;

  /*! @brief Pointer to the material deck */
  const inp::MaterialDeck *d_deck;

  /*! @brief Store pointer to datamanager */
  data::DataManager* d_dataManager_p;

  /*! @brief Base object for influence function */
  material::pd::BaseInfluenceFn *d_baseInfluenceFn_p;
};

} // namespace pd

} // namespace material

#endif // MATERIAL_PD_BONDMATERIAL_NEW_H
