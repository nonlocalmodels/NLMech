////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef MATERIAL_PD_STATEMATERIAL_H
#define MATERIAL_PD_STATEMATERIAL_H

#include "baseMaterial.h"
#include <vector>

// forward declaration
namespace inp {
struct MaterialDeck;
}

namespace material {

namespace pd {

/*! @brief A Class implementing regularized nonlinear peridynamic model
 *
 * Provides method to compute energy and force using nonlinear state-based
 * model introduced and studied in [Lipton 2018](https://doi.org/10
 * .1007/978-3-319-22977-5_33-1), [Jha and Lipton 2019](https://www
 * .sciencedirect.com/science/article/pii/S0045782519301537).
 *
 * 1. The pairwise energy and hydrostatic energy are
 * \f[ E_{bond} = \int_D \frac{1}{|B_\epsilon(x)|} \int_{B_\epsilon(x)} |y-x|
 * W_{bond}(S(y,x; u)) dy dx \f]
 * \f[ E_{state} = \int_D W_{state}(\theta(x;u)) dx \f]
 * where \f$ B_\epsilon(x) \f$ is a ball (circle in 2-d) centered at \f$ x\f$
 * of  radius \f$ \epsilon\f$, \f$ |B_\epsilon(x)| \f$ is the volume (area in
 * 2-d) of ball, \f$ S(y,x;u) = \frac{u(y) - u(x)}{|y-x|} \cdot
 * \frac{y-x}{|y-x|} \f$
 * is the linearized bond strain (assuming small deformation) and \f$ \theta
 * (x;u) \f$ is the volumetric (hydrostatic) strain at material point \f$ x\f$
 * defined as
 * \f[ \theta(x;u) = \frac{1}{|H_\epsilon(x)|} \int_{H_\epsilon(x)}
 * J^\epsilon(|y-x|) S(y,x;u) |y-x| dy. \f]
 * \f$ J^\epsilon(|y-x|)\f$ is the influence function.
 *
 * 2. Bond energy density \f$ W_{bond} \f$ is given by
 * \f[ W_{bond} (S(y,x;u)) = J^\epsilon(|y-x|) \frac{1}{\epsilon |y-x|} \psi
 * (|y-x| S(y,x;u)^2) \f]
 * where \f$ \psi(r) : \textsf{R}^+ \to \textsf{R} \f$ is positive, smooth,
 * concave
 * function with following properties
 * \f[ \lim_{r\to 0^+} \frac{\psi(r)}{r} = \psi'(0), \quad \lim_{r\to \infty}
 * \psi(r) = \psi_\infty < \infty. \f]
 *
 * 3. State energy density \f$ W_{state}\f$ is given by
 * \f[ W_{state}(\theta(x;u)) = \frac{g(\theta(x;u))}{\epsilon^2} \f]
 * where \f$ g : \mathsf{R} \to \mathsf{R} \f$ is the hydrostatic potential
 * function. It can be either quadratic function or can be convex for small
 * strain and concave for large strain.
 *
 * 4. Force at material point \f$ x \f$ due to bond-based potential is given by
 * \f[ f_{bond}(x) = \frac{4}{|B_\epsilon(x)|} \int_{B_\epsilon(x)}
 * \frac{J^\epsilon(|y-x|)}{\epsilon} \psi'(|y-x| S(y,x;u)^2) S(y,x;u)
 * \frac{y-x}{|y-x|} dy.\f]
 *
 * 5. Force at material point \f$ x\f$ due to state-based potential is given by
 * \f[ f_{state}(x) = \frac{1}{|B_\epsilon(x)|} \int_{B_\epsilon(x)}
 * \frac{J^\epsilon(|y-x|)}{\epsilon^2} [g'(\theta(x;u)) + g'(\theta(y;u))]
 * \frac{y-x}{|y-x|} dy.\f]
 *
 * 6. In this class, we will assume
 *  - \f[ \psi(r) = C ( 1-\exp[-\beta r] \f]
 *  - \f[ g(r) = \frac{\bar{C}}{2} r^2 \f]
 * where \f$ C, \beta, \bar{C} \f$ are the peridynamic material parameter
 * determined from the elastic and fracture properties of the material.
 */
class RNPState : public BaseMaterial {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   * @param dim Dimension
   * @param horizon Horizon
   * @param M Moment of influence function
   */
  RNPState(inp::MaterialDeck *deck, const size_t &dim, const double &horizon,
           const double &M);

  /*!
   * @brief Returns energy and force between bond due to pairwise interaction
   *
   * Peridynamic energy at point \f$ x \f$ is
   * \f[ e(x) = \frac{1}{|B_\epsilon(0)|} \int_{B_\epsilon(x)}
   * \frac{J^\epsilon(|y-x|)}{\epsilon} \psi(|y-x|S^2) dy \f]
   * and force at point x is
   * \f[ f(x) = \frac{4}{|B_\epsilon(0)|} \int_{B_\epsilon(x)}
   * \frac{J^\epsilon(|y-x|)}{\epsilon} \psi'(|y-x|S^2) S \frac{y-x}{|y-x|}
   * dy. \f]
   * where \f$\psi(r) = C(1-\exp(-\beta r))\f$.
   *
   * For given initial bond length \f$ r \f$ and bond strain \f$ s\f$, this
   * function returns pair of
   * \f[ \hat{e} =  \frac{J^\epsilon(r)}{\epsilon |B_\epsilon(0)|} \psi(r s^2)
   * \f]
   * and
   * \f[ \hat{f} = \frac{4 J^\epsilon(r) s}{\epsilon |B_\epsilon(0)|}
   * \psi'(r s^2). \f]
   *
   * @param r Reference (initial) bond length
   * @param s Bond strain
   * @param J Influence function at r, i.e. \f$ J^\epsilon(r) \f$
   * @param fs Bond fracture state
   * @return Value Pair of energy and force
   */
  std::pair<double, double> getBondEF(const double &r, const double &s,
                                      const double &J, bool &fs) override;

  /*!
   * @brief Returns energy and force between bond for \a no-fail region
   *
   * Peridynamic energy at point \f$ x \f$ is
   * \f[ e(x) = \frac{1}{|B_\epsilon(0)|} \int_{B_\epsilon(x)}
   * \frac{J^\epsilon(|y-x|)}{\epsilon} \psi(|y-x|S^2) dy \f]
   * and force at point x is
   * \f[ f(x) = \frac{4}{|B_\epsilon(0)|} \int_{B_\epsilon(x)}
   * \frac{J^\epsilon(|y-x|)}{\epsilon} \psi'(|y-x|S^2) S \frac{y-x}{|y-x|}
   * dy. \f]
   * where \f$\psi(r) = C(1-\exp(-\beta r))\f$.
   *
   * For material point \f$ x \f$ in \a no-fail region, we modify the function
   * \f$ \psi \f$ to \f[ \psi(r) = \psi'(0) r, \qquad \psi'(0) = C\beta . \f]
   *
   * For given initial bond length \f$ r \f$ and bond strain \f$ s\f$, this
   * function returns pair of
   * \f[ \hat{e} =  \frac{J^\epsilon(r)}{\epsilon |B_\epsilon(0)|} \psi(r s^2)
   * \f]
   * and
   * \f[ \hat{f} = \frac{4 J^\epsilon(r) s}{\epsilon |B_\epsilon(0)|}
   * \psi'(r s^2). \f]
   *
   * @param r Reference (initial) bond length
   * @param s Bond strain
   * @param J Influence function at r, i.e. \f$ J^\epsilon(r) \f$
   * @return pair Pair of energy and force
   */
  std::pair<double, double> getBondEFNoFail(const double &r, const double &s,
                                            const double &J) override;

  /*!
   * @brief Returns hydrostatic energy density
   *
   * Total hydrostatic energy is given by
   * \f[ \int_D \frac{g(\theta(x;u))}{\epsilon^2}  dx \f]
   * where \f$ g(r) = \bar{C}r^2/2 \f$.
   *
   * Given hydrostatic strain \f$ \theta \f$, this function returns \f$
   * \frac{g(\theta)}{\epsilon^2}\f$.
   *
   * @param theta Hydrostatic strain
   * @return energy Energy density
   */
  double getStateEnergy(const double &theta) override;

  /*!
   * @brief Returns hydrostatic force density
   *
   * Force at material point \f$ x\f$ due to state-based potential is given by
   * \f[ f_{state}(x) = \frac{1}{|B_\epsilon(x)|} \int_{B_\epsilon(x)}
   * \frac{J^\epsilon(|y-x|)}{\epsilon^2} [g'(\theta(x;u)) + g'(\theta(y;u))]
   * \frac{y-x}{|y-x|} dy.\f]
   * Here \f$ g(r) = \frac{\bar{C}}{2} r^2 \f$.
   *
   * Given hydrostatic strain \f$ \theta \f$, influence function \f$ J \f$ at
   * \f$ r\f$, this function returns
   * \f[ \frac{g'(\theta) J}{\epsilon^2|B_\epsilon(x)|} \f]
   *
   * @param theta Hydrostatic strain
   * @param J Influence function at r
   * @return force Force density
   */
  double getStateForce(const double &theta, const double &J) override;

  /*!
   * @brief Returns true if bond contributes to hydrostatic force
   * @param S Bond strain
   * @param r Reference bond length
   * @return bool True/false
   */
  bool doesBondContribToState(const double &S, const double &r) override;

  /*!
   * @brief Returns contribution of bond to hydrostatic strain
   *
   * Hydrostatic strain is given by
   * \f[ \theta(x;u) = \frac{1}{|H_\epsilon(x)|} \int_{H_\epsilon(x)}
   * J^\epsilon(|y-x|) S(y,x;u) |y-x| dy. \f]
   *
   * Given \f$ r\f$, strain \f$ s \f$, this function returns
   * \f[ \frac{J^\epsilon(r) s r}{|H_\epsilon(x)|}. \f]
   *
   * @param S Bond strain
   * @param r Reference (initial) bond length
   * @param J Influence function at r
   * @return strain Contribution to hydrostatic strain
   */
  double getBondContribToHydroStrain(const double &S, const double &r,
                                     const double &J) override;

  /*!
   * @brief Returns critical bond strain
   *
   * @param r Reference length of bond
   * @return strain Critical strain
   */
  double getSc(const double &r) override;

private:
  /*!
   * @brief Computes rnp material parameters from elastic constants
   *
   * Either Young's modulus E or bulk modulus K, and Poisson ratio \f$ \nu \f$,
   * and either critical energy release rate Gc or critical stress intensity
   * factor KIc are needed. We first compute lame
   * parameters \f$ \lambda, \mu \f$ from K (or E) and \f$ \nu \f$. With lame
   * parameters, we use following formula, see Equations (94), (95)
   * & (97) of [Lipton 2018](https://doi.org/10.1007/978-3-319-22977-5_33-1)
   * - if \f$ d=2 \f$
   * \f[ \mu = \frac{\psi'(0)}{4} M_2, \qquad \lambda = \mu + \frac{g''(0)}{2}
   * M_2^2, \qquad Gc = \frac{4\psi_{\infty}}{\pi} M_2. \f]
   * - if \f$ d=3\f$
   * \f[ \mu = \frac{\psi'(0)}{5} M_3, \qquad \lambda = \mu + \frac{g''(0)
   * }{2}M_3^2, \qquad Gc = \frac{3\psi_{\infty}}{2} M_3. \f]
   *
   * Where \f$M_2, M_3\f$ are defined by
   * \f[ M_2 =\int_0^1 r^2 J(r) dr, \qquad M_3 = \int_0^1 r^3 J(r) dr. \f]
   *
   * For potential function \f$ \psi(r) = c ( 1-\exp[-\beta r]), g(r) =
   * \bar{C} r^2/2\f$, we have \f$
   * \psi'(0) = c\beta, \psi_{\infty} = c, g''(0) = \bar{C} \f$. Thus, the
   * values of \f$ c, \beta, \bar{C} \f$ are given by
   * - if \f$ d=2 \f$
   * \f[ c = \frac{\pi G_c}{4} \frac{1}{M_2}, \qquad
   *     \beta = \frac{4 \mu}{c} \frac{1}{M_2}, \qquad \bar{C} = \frac{2
   *     (\lambda - \mu)}{M_2^2} .\f]
   * - if \f$ d=3 \f$
   * \f[ c = \frac{2 G_c}{3} \frac{1}{M_3}, \qquad
   *     \beta = \frac{5 \lambda}{c} \frac{1}{M_3}, \qquad \bar{C} = \frac{2
   *     (\lambda - \mu)}{M_3^2} .\f]
   *
   * @param deck Input material deck
   * @param M Moment of influence function
   */
  void computeParameters(inp::MaterialDeck *deck, const double &M);

  /*!
   * @brief Computes elastic and fracture properties from the rnp material
   * parameters
   *
   * This function is does opposite of computeParameters(). From peridynamic
   * material properties, it uses the relation between lame parameters and
   * peridynamic parameters (see computeParameters() description) to compute
   * lame parameters, and from lame parameters it computes the elastic
   * constants.
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

  /*! @brief Parameter \f$ \bar{C} \f$ */
  double d_barC;

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

  /*!
   * @brief Flag which indicates if the broken bond contributes to
   * state-based interaction
   */
  bool d_stateContributionFromBrokenBond;
};

} // namespace pd

} // namespace material

#endif // MATERIAL_PD_BONDMATERIAL_H
