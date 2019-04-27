// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef MATERIAL_PD_BONDMATERIAL_H
#define MATERIAL_PD_BONDMATERIAL_H

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
 * Provides method to compute energy and force using nonlinear bond-based
 * model introduced and studied in \b Lipton2016, \b JhaLipton208, \b
 * JhaLipton2019 **TODO** Add link to papers
 *
 * Given bond length \f$ r\f$ and bond strain \f$ s\f$, the pairwise energy
 * and force is given by
 * \f[ \psi(rs^2) = C(1 - \exp[-\beta r s^2]), \quad \psi'(rs^2) = C\beta
 * \exp[ -\beta rs^2] \f]
 * where \f$ C,\beta \f$ are material parameters, \f$\psi'(r) = d\psi(r)/dr
 * \f$.
 *
 * \f$ C,\beta\f$ are determined by the given material properties such as
 * Young's modulus and critical energy release rate.
 */
class RNPBond : public BaseMaterial {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   * @param dim Dimension
   * @param horizon Horizon
   * @param moment_inf Moment of influence function
   */
  RNPBond(inp::MaterialDeck *deck, const size_t &dim, const double &horizon,
          const double &moment_inf);

  /*!
   * @brief Returns energy and force between bond
   *
   * Peridynamic energy at point x is
   * \f[ e(x) = \frac{1}{\epsilon |B_\epsilon(0)|} \int_{B_\epsilon(x)}
   * J^\epsilon (|y-x|) \psi(|y-x|S^2) dy \f]
   * and force at point x is
   * \f[ f(x) = \frac{2}{\epsilon |B_\epsilon(0)|} \int_{B_\epsilon(x)}
   * J^\epsilon (|y-x|) \frac{\partial_S \psi(|y-x|S^2)}{|y-x|}
   * \frac{y-x}{|y-x|} dy. \f]
   * Above can be written as
   * \f[ f(x) = \frac{4}{\epsilon |B_\epsilon(0)|} \int_{B_\epsilon(x)}
   * J^\epsilon (|y-x|) \psi'(|y-x|S^2) S\frac{y-x}{|y-x|} dy \f]
   *
   * where \f$\psi(r) = C(1-\exp(-\beta r))\f$.
   *
   * This function returns pair of
   * \f[ \hat{e} =  \frac{1}{\epsilon |B_\epsilon(0)|} \psi(|y-x| S^2) \f]
   * and \f[ \hat{f} = \frac{4 S}{\epsilon |B_\epsilon(0)|} \psi'(|y-x|S^2). \f]
   *
   * @param r Reference (initial) bond length
   * @param s Bond strain
   * @param influence Value of influence function at r
   * @param fs Bond fracture state
   * @return Value Pair of energy and force
   */
  std::pair<double, double> getBondEF(const double &r, const double &s,
                                      const double &influence,
                                      bool &fs) override;

  /*!
   * @brief Returns energy and force between bond for \a no-fail region
   *
   * Peridynamic energy at point x is
   * \f[ e(x) = \frac{1}{\epsilon |B_\epsilon(0)|} \int_{B_\epsilon(x)}
   * J^\epsilon (|y-x|) \psi(|y-x|S^2) dy \f]
   * and force at point x is
   * \f[ f(x) = \frac{2}{\epsilon |B_\epsilon(0)|} \int_{B_\epsilon(x)}
   * J^\epsilon (|y-x|) \frac{\partial_S \psi(|y-x|S^2)}{|y-x|}
   * \frac{y-x}{|y-x|} dy. \f]
   * Above can be written as
   * \f[ f(x) = \frac{4}{\epsilon |B_\epsilon(0)|} \int_{B_\epsilon(x)}
   * J^\epsilon (|y-x|) \psi'(|y-x|S^2) S\frac{y-x}{|y-x|} dy \f]
   *
   * where \f$\psi(r) = C(1-\exp(-\beta r))\f$.
   *
   * For material point x in \a no-fail region, we modify the function \f$
   * \psi \f$ to \f[ \psi(r) = \psi'(0) r, \qquad \psi'(0) = C\beta . \f]
   *
   * This function returns pair of
   * \f[ \hat{e} =  \frac{1}{\epsilon |B_\epsilon(0)|} \psi(|y-x| S^2) \f]
   * and \f[ \hat{f} = \frac{4 S}{\epsilon |B_\epsilon(0)|} \psi'(|y-x|S^2) \f]
   * where \f$\psi(r)\f$ is modified as described above for points in \a
   * no-fail region.
   *
   * @param r Reference (initial) bond length
   * @param s Bond strain
   * @param influence Value of influence function at r
   * @return Value Pair of energy and force
   */
  std::pair<double, double> getBondEFNoFail(const double &r, const double &s,
                                            const double &influence) override;

  /*!
   * @brief Returns critical bond strain
   *
   * @param r Reference length of bond
   * @return Value Critical strain
   */
  double getSc(const double &r) override;

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
   * .8) of \b Lipton2016 **TODO** add links to paper
   * - if \f$ d=2 \f$
   * \f[ \lambda = \mu = \frac{f'(0)}{4} \int_0^1 r^2 J(r) dr f]
   * - if \f$ d=3\f$
   * \f[ \lambda = \mu = \frac{f'(0)}{5} \int_0^1 r^3 J(r) dr \f]
   * and
   * - if \f$ d=2 \f$
   * \f[ Gc = \frac{4 f_{\infty}}{\pi} \int_0^1 r^2 J(r) dr \f]
   * - if \f$ d=3 \f$
   * \f[ Gc = \frac{3 f_{\infty}}{2} \int_0^1 r^2 J(r) dr \f]
   *
   * For potential function \f$ f(r) = c ( 1-\exp[-\beta r])\f$, we have \f$
   * f'(0) = c\beta, f_{\infty} = c \f$. Thus, the values of \f$ c, \beta \f$
   * are given by
   * - if \f$ d=2 \f$
   * \f[ c = \frac{}{} \f]
   * @param deck Input material deck
   * @param moment_inf Moment of influence function
   */
  void computeParameters(inp::MaterialDeck *deck, const double &moment_inf);

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
};

} // namespace pd

} // namespace material

#endif // MATERIAL_PD_BONDMATERIAL_H
