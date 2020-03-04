////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef MATERIAL_PD_MATERIAL_H
#define MATERIAL_PD_MATERIAL_H

#include "inp/decks/materialDeck.h"
#include "material/pd/influenceFn.h"
#include "util/point.h" // definition of Point3
#include "util/compare.h"
#include "util/utilIO.h"
#include <limits>
#include <string>
#include <vector>

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
  explicit Material(std::string name = "") : d_name(name) {}

  /*!
   * @brief Returns name of the material
   * @return string Name
   */
  std::string name() { return d_name; }

  /*!
   * @brief Returns true if state-based potential is active
   * @return bool True/false
   */
  virtual bool isStateActive() const = 0;

  /*!
   * @brief Returns energy and force between bond due to pairwise interaction
   *
   * @param r Reference (initial) bond length
   * @param s Bond strain
   * @param fs Bond fracture state
   * @return value Pair of energy and force
   */
  virtual std::pair<double, double>
  getBondEF(const double &r, const double &s, bool &fs,
            const bool &break_bonds) const = 0;

  /*!
   * @brief Returns energy and force between bond due to state-based model
   *
   * @param r Reference (initial) bond length
   * @param s Bond strain
   * @param fs Bond fracture state
   * @param mx Weighted volume at node
   * @param thetax Dilation
   * @return value Pair of energy and force
   */
  virtual std::pair<double, double>
  getBondEF(const double &r, const double &s, bool &fs, const double
  &mx, const double &thetax) const = 0;

  virtual util::Point3 getBondForceDirection(const util::Point3 &dx,
                                             const util::Point3 &du) const = 0;

  /*!
   * @brief Returns the bond strain
   * @param dx Reference bond vector
   * @param du Difference of displacement
   * @return strain Bond strain \f$ S = \frac{du \cdot dx}{|dx|^2} \f$
   */
  virtual double getS(const util::Point3 &dx, const util::Point3 &du) const = 0;

  /*!
   * @brief Returns critical bond strain
   *
   * @param r Reference length of bond
   * @return strain Critical strain
   */
  virtual double getSc(const double &r) const = 0;

  /*!
   * @brief Returns the density of the material
   * @return density Density of the material
   */
  virtual double getDensity() const = 0;

  /*!
   * @brief Returns the value of influence function
   *
   * @param r Reference (initial) bond length
   * @return value Influence function at r
   */
  virtual double getInfFn(const double &r) const = 0;

  /*!
   * @brief Returns the moment of influence function
   *
   * If \f$ J(r) \f$ is the influence function for \f$ r\in [0,1)\f$ then \f$
   * i^{th}\f$ moment is given by \f[ M_i = \int_0^1 J(r) r^i dr. \f]
   *
   * @param i ith moment
   * @return value Moment
   */
  virtual double getMoment(const size_t &i) const = 0;

  /*!
   * @brief Returns the material deck
   * @return deck RnpMaterial deck
   */
  virtual const inp::MaterialDeck *getMaterialDeckP() const = 0;
  virtual inp::MaterialDeck *getMaterialDeckP() = 0;

  /*!
   * @brief Returns horizon
   *
   * @return horizon Horizon
   */
  virtual double getHorizon() const = 0;

  virtual std::string printStr(int nt, int lvl) const {

    auto tabS = util::io::getTabS(nt);
    std::ostringstream oss;
    oss << tabS << "------- particle::Material --------" << std::endl
        << std::endl;
    oss << tabS << "Abstract class of peridynamic materials" << std::endl;
    oss << tabS << "See RnpMaterial and PmbMaterial for implementation"
        << std::endl;
    oss << tabS << std::endl;

    return oss.str();
  }
  virtual void print(int nt, int lvl) const { std::cout << printStr(nt, lvl); }
  virtual void print() const { print(0, 0); }

private:
  /*! @brief Name of the material */
  std::string d_name;
};

/*!
 * @brief A class providing methods to compute energy density and force of
 * peridynamic material
 */
class RnpMaterial : public Material {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   * @param dim Dimension
   * @param horizon Horizon
   */
  RnpMaterial(inp::MaterialDeck *deck, const size_t &dim, const double &horizon)
      : Material("RNPBond"), d_deck_p(new inp::MaterialDeck(*deck)), d_horizon(horizon),
        d_dimension(dim), d_density(d_deck_p->d_density), d_infFn_p(nullptr),
        d_C(0.), d_beta(0.), d_rbar(0.), d_invFactor(0.),
        d_factorSc(d_deck_p->d_checkScFactor),
        d_irrevBondBreak(d_deck_p->d_irreversibleBondBreak) {

    // create influence function
    if (d_deck_p->d_influenceFnType == 0)
      d_infFn_p = new material::pd::ConstInfluenceFn(
          d_deck_p->d_influenceFnParams, d_dimension);
    else if (d_deck_p->d_influenceFnType == 1)
      d_infFn_p = new material::pd::LinearInfluenceFn(
          d_deck_p->d_influenceFnParams, d_dimension);
    else if (d_deck_p->d_influenceFnType == 2)
      d_infFn_p = new material::pd::GaussianInfluenceFn(
          d_deck_p->d_influenceFnParams, d_dimension);
    else {
      std::cerr << "Error: Influence function type unknown in input file.\n";
      exit(1);
    }

    if (dim == 1)
      d_invFactor = std::pow(horizon, 2) * 2.;
    else if (dim == 2)
      d_invFactor = std::pow(horizon, 3) * M_PI;
    else if (dim == 3)
      d_invFactor = std::pow(horizon, 4) * 4. * M_PI / 3.;

    // check if we need to compute the material parameters
    if (d_deck_p->d_computeParamsFromElastic)
      computeParameters(d_deck_p, d_infFn_p->getMoment(dim));
    else {
      d_C = d_deck_p->d_bondPotentialParams[0];
      d_beta = d_deck_p->d_bondPotentialParams[1];
      d_rbar = std::sqrt(0.5 / d_beta);
      computeMaterialProperties(d_deck_p, d_infFn_p->getMoment(dim));
    }
  };

  /*!
   * @brief Returns true if state-based potential is active
   * @return bool True/false
   */
  bool isStateActive() const override { return false; };

  /*!
   * @brief Returns energy and force between bond due to pairwise interaction
   *
   * Peridynamic energy at point \f$ x \f$ is
   * \f[ e(x) = \frac{1}{|B_\epsilon(0)|} \int_{B_\epsilon(x)}
   * \frac{J^\epsilon(|y-x|)}{\epsilon} \psi(|y-x|S^2) dy \f]
   * and force at point x is
   * \f[ f(x) = \frac{4}{|B_\epsilon(0)|} \int_{B_\epsilon(x)}
   * \frac{J^\epsilon(|y-x|)}{\epsilon} \psi'(|y-x|S^2) S \frac{y-x}{|y-x|}
   * dy, \f]
   * where \f$ \psi(r) = C(1-\exp(-\beta r))\f$.
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
   * @param fs Bond fracture state
   * @return value Pair of energy and force
   */
  std::pair<double, double> getBondEF(const double &r, const double &s,
                                      bool &fs,
                                      const bool &break_bonds) const override {

    if (break_bonds) {
      // check if fracture state of the bond need to be updated
      if (d_irrevBondBreak && !fs &&
          util::compare::definitelyGreaterThan(std::abs(s),
                                               d_factorSc * getSc(r)))
        fs = true;

      // if bond is not fractured, return energy and force from nonlinear
      // potential otherwise return energy of fractured bond, and zero force
      if (!fs)
        return std::make_pair(
            getInfFn(r) * d_C *
                (1. - std::exp(-d_beta * r * s * s)) / d_invFactor,
            getInfFn(r) * 4. * s * d_C * d_beta *
                std::exp(-d_beta * r * s * s) / d_invFactor);
      else
        return std::make_pair(getInfFn(r) * d_C / d_invFactor, 0.);
    } else {
      return std::make_pair(getInfFn(r) * d_C * d_beta *
                                r * s * s / d_invFactor,
                            getInfFn(r) * 4. * s * d_C *
                                d_beta / d_invFactor);
    }
  };

  /*!
   * @brief Returns energy and force between bond due to state-based model
   *
   * @param r Reference (initial) bond length
   * @param s Bond strain
   * @param fs Bond fracture state
   * @param mx Weighted volume at node
   * @param thetax Dilation
   * @return value Pair of energy and force
   */
  std::pair<double, double>
  getBondEF(const double &r, const double &s, bool &fs, const double
  &mx, const double &thetax) const override {

    return this->getBondEF(r, s, fs, true);
  };

  util::Point3 getBondForceDirection(const util::Point3 &dx,
                                     const util::Point3 &du) const override {
    return dx / dx.length();
  };

  /*!
   * @brief Returns the bond strain
   * @param dx Reference bond vector
   * @param du Difference of displacement
   * @return strain Bond strain \f$ S = \frac{du \cdot dx}{|dx|^2} \f$
   */
  double getS(const util::Point3 &dx, const util::Point3 &du) const override {
    return dx.dot(du) / dx.dot(dx);
  };

  /*!
   * @brief Returns critical bond strain
   *
   * @param r Reference length of bond
   * @return strain Critical strain
   */
  double getSc(const double &r) const override {
    return d_rbar / std::sqrt(r);
  };

  /*!
   * @brief Returns the density of the material
   * @return density Density of the material
   */
  double getDensity() const override { return d_density; };

  /*!
   * @brief Returns the value of influence function
   *
   * @param r Reference (initial) bond length
   * @return value Influence function at r
   */
  double getInfFn(const double &r) const override {
    return d_infFn_p->getInfFn(r / d_horizon);
  };

  /*!
   * @brief Returns the moment of influence function
   *
   * If \f$ J(r) \f$ is the influence function for \f$ r\in [0,1)\f$ then \f$
   * i^{th}\f$ moment is given by \f[ M_i = \int_0^1 J(r) r^i dr. \f]
   *
   * @param i ith moment
   * @return value Moment
   */
  double getMoment(const size_t &i) const override {
    return d_infFn_p->getMoment(i);
  };

  /*!
   * @brief Returns the material deck
   * @return deck RnpMaterial deck
   */
  const inp::MaterialDeck *getMaterialDeckP() const override {
    return d_deck_p;
  };
  inp::MaterialDeck *getMaterialDeckP() override { return d_deck_p; };

  /*!
   * @brief Returns horizon
   *
   * @return horizon Horizon
   */
  double getHorizon() const override { return d_horizon; };

  /*!
   * @brief Print information about this object
   */
  std::string printStr(int nt, int lvl) const override {

    auto tabS = util::io::getTabS(nt);
    std::ostringstream oss;
    oss << tabS << "------- particle::RnpMaterial --------" << std::endl
        << std::endl;
    oss << tabS << "Material deck address = " << d_deck_p << std::endl;
    oss << tabS << "Material deck info: " << std::endl;
    oss << d_deck_p->printStr(nt + 1, lvl);
    oss << tabS << "State active = " << 0 << std::endl;
    oss << tabS << "Horizon = " << d_horizon << std::endl;
    oss << tabS << "Dimension = " << d_dimension << std::endl;
    oss << tabS << "Influence fn address = " << d_infFn_p << std::endl;
    oss << tabS << "Influence fn info: " << std::endl;
    oss << d_infFn_p->printStr(nt + 1, lvl);
    oss << tabS << "Peridynamic parameters: " << std::endl;
    oss << tabS << "  C = " << d_C << std::endl;
    oss << tabS << "  beta = " << d_beta << std::endl;
    oss << tabS << "  r_bar = " << d_rbar << std::endl;
    oss << tabS << "  inv_factor = " << d_invFactor << std::endl;
    oss << tabS << "  factor_Sc = " << d_factorSc << std::endl;
    oss << tabS << "  irrev_bond_breaking = " << d_irrevBondBreak << std::endl;
    oss << tabS << std::endl;

    return oss.str();
  }

  void print(int nt, int lvl) const override {
    std::cout << printStr(nt, lvl);
  }
  void print() const override { print(0, 0); }

private:
  void computeParameters(inp::MaterialDeck *deck, const double &M) {
    //
    // Need following elastic and fracture properties
    // 1. E or K
    // 2. Gc or KIc
    // For bond-based, Poisson's ratio is fixed to 1/4
    //
    if (util::compare::definitelyLessThan(deck->d_matData.d_E, 0.) &&
        util::compare::definitelyLessThan(deck->d_matData.d_K, 0.)) {
      std::cerr << "Error: Require either Young's modulus E or Bulk modulus K"
                   " to compute the RNP bond-based peridynamic parameters.\n";
      exit(1);
    }
    if (util::compare::definitelyGreaterThan(deck->d_matData.d_E, 0.) &&
        util::compare::definitelyGreaterThan(deck->d_matData.d_K, 0.)) {
      std::cout << "Warning: Both Young's modulus E and Bulk modulus K are "
                   "provided.\n";
      std::cout << "Warning: To compute the RNP bond-based peridynamic "
                   "parameters, we only require one of those.\n";
      std::cout
          << "Warning: Selecting Young's modulus to compute parameters.\n";
    }

    if (util::compare::definitelyLessThan(deck->d_matData.d_Gc, 0.) &&
        util::compare::definitelyLessThan(deck->d_matData.d_KIc, 0.)) {
      std::cerr << "Error: Require either critical energy release rate Gc or "
                   "critical stress intensity factor KIc to compute the RNP "
                   "bond-based peridynamic parameters.\n";
      exit(1);
    } else if (util::compare::definitelyGreaterThan(deck->d_matData.d_Gc, 0.) &&
               util::compare::definitelyGreaterThan(deck->d_matData.d_KIc,
                                                    0.)) {
      std::cout << "Warning: Both critical energy release rate Gc and critical "
                   "stress intensity factor KIc are provided.\n";
      std::cout << "Warning: To compute the RNP bond-based peridynamic "
                   "parameters, we only require one of those.\n";
      std::cout << "Warning: Selecting critical energy release rate Gc to "
                   "compute parameters.\n";
    }

    // set Poisson's ratio to 1/4
    deck->d_matData.d_nu = 0.25;

    // compute E if not provided or K if not provided
    if (deck->d_matData.d_E > 0.)
      deck->d_matData.d_K =
          deck->d_matData.toK(deck->d_matData.d_E, deck->d_matData.d_nu);

    if (deck->d_matData.d_K > 0. && deck->d_matData.d_E < 0.)
      deck->d_matData.d_E =
          deck->d_matData.toE(deck->d_matData.d_K, deck->d_matData.d_nu);

    if (deck->d_matData.d_Gc > 0.)
      deck->d_matData.d_KIc = deck->d_matData.toKIc(
          deck->d_matData.d_Gc, deck->d_matData.d_nu, deck->d_matData.d_E);

    if (deck->d_matData.d_KIc > 0. && deck->d_matData.d_Gc < 0.)
      deck->d_matData.d_Gc = deck->d_matData.toGc(
          deck->d_matData.d_KIc, deck->d_matData.d_nu, deck->d_matData.d_E);

    // compute lame parameter
    deck->d_matData.d_lambda =
        deck->d_matData.toLambdaE(deck->d_matData.d_E, deck->d_matData.d_nu);
    deck->d_matData.d_G =
        deck->d_matData.toGE(deck->d_matData.d_E, deck->d_matData.d_nu);
    deck->d_matData.d_mu = deck->d_matData.d_G;

    // compute peridynamic parameters
    if (d_dimension == 2) {
      d_C = M_PI * deck->d_matData.d_Gc / (4. * M);
      d_beta = 4. * deck->d_matData.d_lambda / (d_C * M);
    } else if (d_dimension == 3) {
      d_C = 2. * deck->d_matData.d_Gc / (3. * M);
      d_beta = 5. * deck->d_matData.d_lambda / (d_C * M);
    }

    d_rbar = std::sqrt(0.5 / d_beta);
  };

  void computeMaterialProperties(inp::MaterialDeck *deck, const double &M) {
    // set Poisson's ratio to 1/4
    deck->d_matData.d_nu = 0.25;

    // compute peridynamic parameters
    if (d_dimension == 2) {
      deck->d_matData.d_Gc = 4. * M * d_C / M_PI;
      deck->d_matData.d_lambda = d_C * M * d_beta / 4.;
    } else if (d_dimension == 3) {
      deck->d_matData.d_Gc = 3. * M * d_C / 2.;
      deck->d_matData.d_lambda = d_C * M * d_beta / 5.;
    }
    deck->d_matData.d_mu = deck->d_matData.d_lambda;
    deck->d_matData.d_G = deck->d_matData.d_lambda;
    deck->d_matData.d_E = deck->d_matData.toELambda(deck->d_matData.d_lambda,
                                                    deck->d_matData.d_nu);
    deck->d_matData.d_K =
        deck->d_matData.toK(deck->d_matData.d_E, deck->d_matData.d_nu);
    deck->d_matData.d_KIc = deck->d_matData.toKIc(
        deck->d_matData.d_Gc, deck->d_matData.d_nu, deck->d_matData.d_E);
  };

private:
  /*! @brief Input material deck */
  inp::MaterialDeck *d_deck_p;

  /*! @brief Horizon */
  double d_horizon;

  /*! @brief Dimension */
  size_t d_dimension;

  /*! @brief Density */
  double d_density;

  /*! @brief Base object for influence function */
  material::pd::BaseInfluenceFn *d_infFn_p;

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

/*!
 * @brief A class providing methods to compute energy density and force of
 * peridynamic material
 */
class PmbMaterial : public Material {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   * @param dim Dimension
   * @param horizon Horizon
   */
  PmbMaterial(inp::MaterialDeck *deck, const size_t &dim, const double &horizon)
      : Material("PMBBond"), d_deck_p(new inp::MaterialDeck(*deck)), d_horizon(horizon),
        d_dimension(dim), d_density(d_deck_p->d_density), d_infFn_p(nullptr),
        d_c(0.), d_s0(0.) {

    // In case of RNP model, we compute certain parameters of influence
    // function so that d-order moment is 1
    // However, this is not expected for other models. So below we pass the
    // value which is considered for the standard influence functions
    // create influence function
    if (d_deck_p->d_influenceFnType == 0)
      d_infFn_p = new material::pd::ConstInfluenceFn(
          {1.}, d_dimension);
    else if (d_deck_p->d_influenceFnType == 1)
      d_infFn_p = new material::pd::LinearInfluenceFn(
          {1., -1.}, d_dimension);
    else if (d_deck_p->d_influenceFnType == 2)
      d_infFn_p = new material::pd::GaussianInfluenceFn(
          d_deck_p->d_influenceFnParams, d_dimension);
    else {
      std::cerr << "Error: Influence function type unknown in input file.\n";
      exit(1);
    }

    // check if we need to compute the material parameters
    if (d_deck_p->d_computeParamsFromElastic)
      computeParameters(d_deck_p, d_infFn_p->getMoment(dim));
    else {
      d_c = d_deck_p->d_bondPotentialParams[0];
      d_s0 = d_deck_p->d_bondPotentialParams[1];
    }
  };

  /*!
   * @brief Returns true if state-based potential is active
   * @return bool True/false
   */
  bool isStateActive() const override { return false; };

  /*!
   * @brief Returns energy and force between bond due to pairwise interaction
   *
   * @param r Reference (initial) bond length
   * @param s Bond strain
   * @param fs Bond fracture state
   * @return value Pair of energy and force
   */
  std::pair<double, double> getBondEF(const double &r, const double &s,
                                      bool &fs,
                                      const bool &break_bonds) const override {

    if (!break_bonds)
      return std::make_pair(getInfFn(r) * 0.5 * d_c * s *
                                s * r,
                            getInfFn(r) * d_c * s);

    // check if fracture state of the bond need to be updated
    if (!fs &&
        util::compare::definitelyGreaterThan(std::abs(s), d_s0 + 1.0e-10))
      fs = true;

    // if bond is not fractured, return energy and force from nonlinear
    // potential otherwise return energy of fractured bond, and zero force
    if (!fs)
      return std::make_pair(getInfFn(r) * 0.5 * d_c * s *
                                s * r,
                            getInfFn(r) * d_c * s);
    else
      return std::make_pair(
          getInfFn(r) * 0.5 * d_c * d_s0 * d_s0 * r, 0.);
  };

  /*!
   * @brief Returns energy and force between bond due to state-based model
   *
   * @param r Reference (initial) bond length
   * @param s Bond strain
   * @param fs Bond fracture state
   * @param mx Weighted volume at node
   * @param thetax Dilation
   * @return value Pair of energy and force
   */
  std::pair<double, double>
  getBondEF(const double &r, const double &s, bool &fs, const double
  &mx, const double &thetax) const override {

    return this->getBondEF(r, s, fs, true);
  };

  util::Point3 getBondForceDirection(const util::Point3 &dx,
                                     const util::Point3 &du) const override {
    return (dx + du) / (dx + du).length();
  };

  /*!
   * @brief Returns the bond strain
   * @param dx Reference bond vector
   * @param du Difference of displacement
   * @return strain Bond strain \f$ S = \frac{du \cdot dx}{|dx|^2} \f$
   */
  double getS(const util::Point3 &dx, const util::Point3 &du) const override {
    return ((dx + du).length() - dx.length()) / dx.length();
  };

  /*!
   * @brief Returns critical bond strain
   *
   * @param r Reference length of bond
   * @return strain Critical strain
   */
  double getSc(const double &r) const override { return d_s0; };

  /*!
   * @brief Returns the density of the material
   * @return density Density of the material
   */
  double getDensity() const override { return d_density; };

  /*!
   * @brief Returns the value of influence function
   *
   * @param r Reference (initial) bond length
   * @return value Influence function at r
   */
  double getInfFn(const double &r) const override {
    return d_infFn_p->getInfFn(r / d_horizon);
  };

  /*!
   * @brief Returns the moment of influence function
   *
   * If \f$ J(r) \f$ is the influence function for \f$ r\in [0,1)\f$ then \f$
   * i^{th}\f$ moment is given by \f[ M_i = \int_0^1 J(r) r^i dr. \f]
   *
   * @param i ith moment
   * @return value Moment
   */
  double getMoment(const size_t &i) const override {
    return d_infFn_p->getMoment(i);
  };

  /*!
   * @brief Returns the material deck
   * @return deck RnpMaterial deck
   */
  const inp::MaterialDeck *getMaterialDeckP() const override {
    return d_deck_p;
  };
  inp::MaterialDeck *getMaterialDeckP() override { return d_deck_p; };

  /*!
   * @brief Returns horizon
   *
   * @return horizon Horizon
   */
  double getHorizon() const override { return d_horizon; };

  /*!
   * @brief Print information about this object
   */
  std::string printStr(int nt, int lvl) const override {

    auto tabS = util::io::getTabS(nt);
    std::ostringstream oss;
    oss << tabS << "------- particle::PmbMaterial --------" << std::endl
        << std::endl;
    oss << tabS << "Material deck address = " << d_deck_p << std::endl;
    oss << tabS << "Material deck info: " << std::endl;
    oss << d_deck_p->printStr(nt + 1, lvl);
    oss << tabS << "State active = " << 0 << std::endl;
    oss << tabS << "Horizon = " << d_horizon << std::endl;
    oss << tabS << "Dimension = " << d_dimension << std::endl;
    oss << tabS << "Influence fn address = " << d_infFn_p << std::endl;
    oss << tabS << "Influence fn info: " << std::endl;
    oss << d_infFn_p->printStr(nt + 1, lvl);
    oss << tabS << "Peridynamic parameters: " << std::endl;
    oss << tabS << "  c = " << d_c << std::endl;
    oss << tabS << "  s0 = " << d_s0 << std::endl;
    oss << tabS << std::endl;

    return oss.str();
  };

  void print(int nt, int lvl) const override {
    std::cout << printStr(nt, lvl);
  };
  void print() const override { print(0, 0); };

private:
  void computeParameters(inp::MaterialDeck *deck, const double &M) {
    //
    // Need following elastic and fracture properties
    // 1. E or K
    // 2. Gc or KIc
    // For bond-based, Poisson's ratio is fixed to 1/4
    //
    if (util::compare::definitelyLessThan(deck->d_matData.d_E, 0.) &&
        util::compare::definitelyLessThan(deck->d_matData.d_K, 0.)) {
      std::cerr << "Error: Require either Young's modulus E or Bulk modulus K"
                   " to compute the RNP bond-based peridynamic parameters.\n";
      exit(1);
    }
    if (util::compare::definitelyGreaterThan(deck->d_matData.d_E, 0.) &&
        util::compare::definitelyGreaterThan(deck->d_matData.d_K, 0.)) {
      std::cout << "Warning: Both Young's modulus E and Bulk modulus K are "
                   "provided.\n";
      std::cout << "Warning: To compute the RNP bond-based peridynamic "
                   "parameters, we only require one of those.\n";
      std::cout
          << "Warning: Selecting Young's modulus to compute parameters.\n";
    }

    if (util::compare::definitelyLessThan(deck->d_matData.d_Gc, 0.) &&
        util::compare::definitelyLessThan(deck->d_matData.d_KIc, 0.)) {
      std::cerr << "Error: Require either critical energy release rate Gc or "
                   "critical stress intensity factor KIc to compute the RNP "
                   "bond-based peridynamic parameters.\n";
      exit(1);
    } else if (util::compare::definitelyGreaterThan(deck->d_matData.d_Gc, 0.) &&
               util::compare::definitelyGreaterThan(deck->d_matData.d_KIc,
                                                    0.)) {
      std::cout << "Warning: Both critical energy release rate Gc and critical "
                   "stress intensity factor KIc are provided.\n";
      std::cout << "Warning: To compute the RNP bond-based peridynamic "
                   "parameters, we only require one of those.\n";
      std::cout << "Warning: Selecting critical energy release rate Gc to "
                   "compute parameters.\n";
    }

    // set Poisson's ratio to 1/4
    deck->d_matData.d_nu = 0.25;

    // compute E if not provided or K if not provided
    if (deck->d_matData.d_E > 0.)
      deck->d_matData.d_K =
          deck->d_matData.toK(deck->d_matData.d_E, deck->d_matData.d_nu);

    if (deck->d_matData.d_K > 0. && deck->d_matData.d_E < 0.)
      deck->d_matData.d_E =
          deck->d_matData.toE(deck->d_matData.d_K, deck->d_matData.d_nu);

    if (deck->d_matData.d_Gc > 0.)
      deck->d_matData.d_KIc = deck->d_matData.toKIc(
          deck->d_matData.d_Gc, deck->d_matData.d_nu, deck->d_matData.d_E);

    if (deck->d_matData.d_KIc > 0. && deck->d_matData.d_Gc < 0.)
      deck->d_matData.d_Gc = deck->d_matData.toGc(
          deck->d_matData.d_KIc, deck->d_matData.d_nu, deck->d_matData.d_E);

    // compute lame parameter
    deck->d_matData.d_lambda =
        deck->d_matData.toLambdaE(deck->d_matData.d_E, deck->d_matData.d_nu);
    deck->d_matData.d_G =
        deck->d_matData.toGE(deck->d_matData.d_E, deck->d_matData.d_nu);
    deck->d_matData.d_mu = deck->d_matData.d_G;

    // compute peridynamic parameters
    if (d_dimension == 2) {
      d_c = 24.0 * deck->d_matData.d_lambda / (M_PI * std::pow(d_horizon, 3.0));
      d_s0 = std::sqrt(5.0 * M_PI * deck->d_matData.d_Gc /
                       (9.0 * deck->d_matData.d_E * d_horizon));
    } else if (d_dimension == 3) {
      d_c = 24.0 * deck->d_matData.d_lambda / (M_PI * std::pow(d_horizon, 3.0));
      d_s0 = std::sqrt(5.0 * M_PI * deck->d_matData.d_Gc /
                       (9.0 * deck->d_matData.d_E * d_horizon));
    }
  };

private:
  /*! @brief Input material deck */
  inp::MaterialDeck *d_deck_p;

  /*! @brief Horizon */
  double d_horizon;

  /*! @brief Dimension */
  size_t d_dimension;

  /*! @brief Density */
  double d_density;

  /*! @brief Base object for influence function */
  material::pd::BaseInfluenceFn *d_infFn_p;

  /**
   * @name Material parameters
   */
  /**@{*/

  /*! @brief Parameter C */
  double d_c;

  /*! @brief Parameter \f$ \beta \f$ */
  double d_s0;

  /** @}*/
};

/*!
 * @brief A class providing methods to compute energy density and force of
 * peridynamic material
 */
class PdElastic : public Material {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   * @param dim Dimension
   * @param horizon Horizon
   */
  PdElastic(inp::MaterialDeck *deck, const size_t &dim, const double &horizon)
      : Material("PDElasticBond"), d_deck_p(new inp::MaterialDeck(*deck)), d_horizon(horizon),
        d_dimension(dim), d_density(d_deck_p->d_density), d_infFn_p(nullptr),
        d_c(0.) {

    // create influence function
    if (d_deck_p->d_influenceFnType == 0)
      d_infFn_p = new material::pd::ConstInfluenceFn(
          {1.}, d_dimension);
    else if (d_deck_p->d_influenceFnType == 1)
      d_infFn_p = new material::pd::LinearInfluenceFn(
          {1., -1.}, d_dimension);
    else if (d_deck_p->d_influenceFnType == 2)
      d_infFn_p = new material::pd::GaussianInfluenceFn(
          d_deck_p->d_influenceFnParams, d_dimension);
    else {
      std::cerr << "Error: Influence function type unknown in input file.\n";
      exit(1);
    }

    // check if we need to compute the material parameters
    if (d_deck_p->d_computeParamsFromElastic)
      computeParameters(d_deck_p, d_infFn_p->getMoment(dim));
    else
      d_c = d_deck_p->d_bondPotentialParams[0];
  };

  /*!
   * @brief Returns true if state-based potential is active
   * @return bool True/false
   */
  bool isStateActive() const override { return false; };

  /*!
   * @brief Returns energy and force between bond due to pairwise interaction
   *
   * @param r Reference (initial) bond length
   * @param s Bond strain
   * @param fs Bond fracture state
   * @return value Pair of energy and force
   */
  std::pair<double, double> getBondEF(const double &r, const double &s,
                                      bool &fs,
                                      const bool &break_bonds) const override {

    return std::make_pair(getInfFn(r) * 0.5 * d_c * s *
                            s * r,
                          getInfFn(r) * d_c * s);
  };

  /*!
   * @brief Returns energy and force between bond due to state-based model
   *
   * @param r Reference (initial) bond length
   * @param s Bond strain
   * @param fs Bond fracture state
   * @param mx Weighted volume at node
   * @param thetax Dilation
   * @return value Pair of energy and force
   */
  std::pair<double, double>
  getBondEF(const double &r, const double &s, bool &fs, const double
  &mx, const double &thetax) const override {

    return this->getBondEF(r, s, fs, true);
  };

  util::Point3 getBondForceDirection(const util::Point3 &dx,
                                     const util::Point3 &du) const override {
    return (dx + du) / (dx + du).length();
  };

  /*!
   * @brief Returns the bond strain
   * @param dx Reference bond vector
   * @param du Difference of displacement
   * @return strain Bond strain \f$ S = \frac{du \cdot dx}{|dx|^2} \f$
   */
  double getS(const util::Point3 &dx, const util::Point3 &du) const override {
    return ((dx + du).length() - dx.length()) / dx.length();
  };

  /*!
   * @brief Returns critical bond strain
   *
   * @param r Reference length of bond
   * @return strain Critical strain
   */
  double getSc(const double &r) const override { return std::numeric_limits<double>::max(); };

  /*!
   * @brief Returns the density of the material
   * @return density Density of the material
   */
  double getDensity() const override { return d_density; };

  /*!
   * @brief Returns the value of influence function
   *
   * @param r Reference (initial) bond length
   * @return value Influence function at r
   */
  double getInfFn(const double &r) const override {
    return d_infFn_p->getInfFn(r / d_horizon);
  };

  /*!
   * @brief Returns the moment of influence function
   *
   * If \f$ J(r) \f$ is the influence function for \f$ r\in [0,1)\f$ then \f$
   * i^{th}\f$ moment is given by \f[ M_i = \int_0^1 J(r) r^i dr. \f]
   *
   * @param i ith moment
   * @return value Moment
   */
  double getMoment(const size_t &i) const override {
    return d_infFn_p->getMoment(i);
  };

  /*!
   * @brief Returns the material deck
   * @return deck RnpMaterial deck
   */
  const inp::MaterialDeck *getMaterialDeckP() const override {
    return d_deck_p;
  };
  inp::MaterialDeck *getMaterialDeckP() override { return d_deck_p; };

  /*!
   * @brief Returns horizon
   *
   * @return horizon Horizon
   */
  double getHorizon() const override { return d_horizon; };

  /*!
   * @brief Print information about this object
   */
  std::string printStr(int nt, int lvl) const override {

    auto tabS = util::io::getTabS(nt);
    std::ostringstream oss;
    oss << tabS << "------- particle::PdElastic --------" << std::endl
        << std::endl;
    oss << tabS << "Material deck address = " << d_deck_p << std::endl;
    oss << tabS << "Material deck info: " << std::endl;
    oss << d_deck_p->printStr(nt + 1, lvl);
    oss << tabS << "State active = " << 0 << std::endl;
    oss << tabS << "Horizon = " << d_horizon << std::endl;
    oss << tabS << "Dimension = " << d_dimension << std::endl;
    oss << tabS << "Influence fn address = " << d_infFn_p << std::endl;
    oss << tabS << "Influence fn info: " << std::endl;
    oss << d_infFn_p->printStr(nt + 1, lvl);
    oss << tabS << "Peridynamic parameters: " << std::endl;
    oss << tabS << "  c = " << d_c << std::endl;
    oss << tabS << std::endl;

    return oss.str();
  };

  void print(int nt, int lvl) const override {
    std::cout << printStr(nt, lvl);
  };
  void print() const override { print(0, 0); };

private:
  void computeParameters(inp::MaterialDeck *deck, const double &M) {
    //
    // Need following elastic and fracture properties
    // 1. E or K
    // For bond-based, Poisson's ratio is fixed to 1/4
    //
    if (util::compare::definitelyLessThan(deck->d_matData.d_E, 0.) &&
        util::compare::definitelyLessThan(deck->d_matData.d_K, 0.)) {
      std::cerr << "Error: Require either Young's modulus E or Bulk modulus K"
                   " to compute the RNP bond-based peridynamic parameters.\n";
      exit(1);
    }
    if (util::compare::definitelyGreaterThan(deck->d_matData.d_E, 0.) &&
        util::compare::definitelyGreaterThan(deck->d_matData.d_K, 0.)) {
      std::cout << "Warning: Both Young's modulus E and Bulk modulus K are "
                   "provided.\n";
      std::cout << "Warning: To compute the RNP bond-based peridynamic "
                   "parameters, we only require one of those.\n";
      std::cout
          << "Warning: Selecting Young's modulus to compute parameters.\n";
    }

    // set Poisson's ratio to 1/4
    deck->d_matData.d_nu = 0.25;

    // compute E if not provided or K if not provided
    if (deck->d_matData.d_E > 0.)
      deck->d_matData.d_K =
          deck->d_matData.toK(deck->d_matData.d_E, deck->d_matData.d_nu);

    if (deck->d_matData.d_K > 0. && deck->d_matData.d_E < 0.)
      deck->d_matData.d_E =
          deck->d_matData.toE(deck->d_matData.d_K, deck->d_matData.d_nu);

    // compute lame parameter
    deck->d_matData.d_lambda =
        deck->d_matData.toLambdaE(deck->d_matData.d_E, deck->d_matData.d_nu);
    deck->d_matData.d_G =
        deck->d_matData.toGE(deck->d_matData.d_E, deck->d_matData.d_nu);
    deck->d_matData.d_mu = deck->d_matData.d_G;

    // compute peridynamic parameters
    if (d_dimension == 2) {
      d_c = 24.0 * deck->d_matData.d_lambda / (M_PI * std::pow(d_horizon, 3.0));
    } else if (d_dimension == 3) {
      // TODO
      //  For PdElastic and PmbMaterial, implement correct formula in 3-d
      d_c = 24.0 * deck->d_matData.d_lambda / (M_PI * std::pow(d_horizon, 3.0));
    }
  };

private:
  /*! @brief Input material deck */
  inp::MaterialDeck *d_deck_p;

  /*! @brief Horizon */
  double d_horizon;

  /*! @brief Dimension */
  size_t d_dimension;

  /*! @brief Density */
  double d_density;

  /*! @brief Base object for influence function */
  material::pd::BaseInfluenceFn *d_infFn_p;

  /**
   * @name Material parameters
   */
  /**@{*/

  /*! @brief Parameter C */
  double d_c;

  /** @}*/
};

/*!
 * @brief A class providing methods to compute energy density and force of
 * peridynamic material
 */
class PdState : public Material {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   * @param dim Dimension
   * @param horizon Horizon
   */
  PdState(inp::MaterialDeck *deck, const size_t &dim, const double &horizon)
      : Material("PDState"), d_deck_p(new inp::MaterialDeck(*deck)), d_horizon
      (horizon),
        d_dimension(dim), d_density(d_deck_p->d_density), d_infFn_p(nullptr),
        d_K(0.), d_G(0.), d_s0(0.) {

    // create influence function
    if (d_deck_p->d_influenceFnType == 0)
      d_infFn_p = new material::pd::ConstInfluenceFn(
          {1.}, d_dimension);
    else if (d_deck_p->d_influenceFnType == 1)
      d_infFn_p = new material::pd::LinearInfluenceFn(
          {1., -1.}, d_dimension);
    else if (d_deck_p->d_influenceFnType == 2)
      d_infFn_p = new material::pd::GaussianInfluenceFn(
          d_deck_p->d_influenceFnParams, d_dimension);
    else {
      std::cerr << "Error: Influence function type unknown in input file.\n";
      exit(1);
    }

    // check if we need to compute the material parameters
    if (d_deck_p->d_computeParamsFromElastic)
      computeParameters(d_deck_p, d_infFn_p->getMoment(dim));
    else {
      d_K = d_deck_p->d_bondPotentialParams[0];
      d_G = d_deck_p->d_bondPotentialParams[1];
      d_s0 = d_deck_p->d_bondPotentialParams[2];
    }
  };

  /*!
   * @brief Returns true if state-based potential is active
   * @return bool True/false
   */
  bool isStateActive() const override { return true; };

  /*!
   * @brief Returns energy and force between bond due to pairwise interaction
   *
   * @param r Reference (initial) bond length
   * @param s Bond strain
   * @param fs Bond fracture state
   * @return value Pair of energy and force
   */
  std::pair<double, double> getBondEF(const double &r, const double &s,
                                      bool &fs,
                                      const bool &break_bonds) const override {

    return {0., 0.};
  };

  /*!
   * @brief Returns energy and force between bond due to state-based model
   *
   * @param r Reference (initial) bond length
   * @param s Bond strain
   * @param fs Bond fracture state
   * @param mx Weighted volume at node
   * @param thetax Dilation
   * @return value Pair of energy and force
   */
  std::pair<double, double>
  getBondEF(const double &r, const double &s, bool &fs, const double
  &mx, const double &thetax) const override {

    if (fs)
      return {0., 0.};

    double J = getInfFn(r);
    double change_length = s * r;

    double alpha = 15. * d_G / mx;
    double factor = (3. * d_K / mx) - alpha / 3.;

    return {0., J * (r * thetax * factor + change_length * alpha)};
  };

  util::Point3 getBondForceDirection(const util::Point3 &dx,
                                     const util::Point3 &du) const override {
    return (dx + du) / (dx + du).length();
  };

  /*!
   * @brief Returns the bond strain
   * @param dx Reference bond vector
   * @param du Difference of displacement
   * @return strain Bond strain \f$ S = \frac{du \cdot dx}{|dx|^2} \f$
   */
  double getS(const util::Point3 &dx, const util::Point3 &du) const override {
    return ((dx + du).length() - dx.length()) / dx.length();
  };

  /*!
   * @brief Returns critical bond strain
   *
   * @param r Reference length of bond
   * @return strain Critical strain
   */
  double getSc(const double &r) const override { return d_s0; };

  /*!
   * @brief Returns the density of the material
   * @return density Density of the material
   */
  double getDensity() const override { return d_density; };

  /*!
   * @brief Returns the value of influence function
   *
   * @param r Reference (initial) bond length
   * @return value Influence function at r
   */
  double getInfFn(const double &r) const override {
    return d_infFn_p->getInfFn(r / d_horizon);
  };

  /*!
   * @brief Returns the moment of influence function
   *
   * If \f$ J(r) \f$ is the influence function for \f$ r\in [0,1)\f$ then \f$
   * i^{th}\f$ moment is given by \f[ M_i = \int_0^1 J(r) r^i dr. \f]
   *
   * @param i ith moment
   * @return value Moment
   */
  double getMoment(const size_t &i) const override {
    return d_infFn_p->getMoment(i);
  };

  /*!
   * @brief Returns the material deck
   * @return deck RnpMaterial deck
   */
  const inp::MaterialDeck *getMaterialDeckP() const override {
    return d_deck_p;
  };
  inp::MaterialDeck *getMaterialDeckP() override { return d_deck_p; };

  /*!
   * @brief Returns horizon
   *
   * @return horizon Horizon
   */
  double getHorizon() const override { return d_horizon; };

  /*!
   * @brief Print information about this object
   */
  std::string printStr(int nt, int lvl) const override {

    auto tabS = util::io::getTabS(nt);
    std::ostringstream oss;
    oss << tabS << "------- particle::PdState --------" << std::endl
        << std::endl;
    oss << tabS << "Material deck address = " << d_deck_p << std::endl;
    oss << tabS << "Material deck info: " << std::endl;
    oss << d_deck_p->printStr(nt + 1, lvl);
    oss << tabS << "State active = " << 1 << std::endl;
    oss << tabS << "Horizon = " << d_horizon << std::endl;
    oss << tabS << "Dimension = " << d_dimension << std::endl;
    oss << tabS << "Influence fn address = " << d_infFn_p << std::endl;
    oss << tabS << "Influence fn info: " << std::endl;
    oss << d_infFn_p->printStr(nt + 1, lvl);
    oss << tabS << "Peridynamic parameters: " << std::endl;
    oss << tabS << "  K = " << d_K << std::endl;
    oss << tabS << "  G = " << d_G << std::endl;
    oss << tabS << "  s0 = " << d_s0 << std::endl;
    oss << tabS << std::endl;

    return oss.str();
  };

  void print(int nt, int lvl) const override {
    std::cout << printStr(nt, lvl);
  };
  void print() const override { print(0, 0); };

private:
  void computeParameters(inp::MaterialDeck *deck, const double &M) {
    //
    // Need following elastic and fracture properties
    // 1. E or K
    // 2. Poisson ratio or shear modulus
    // 3. Gc or KIc
    //
    bool found_E = false;
    bool found_K = false;
    bool found_G = false;
    bool found_nu = false;
    size_t num_props = 0;

    found_E = util::compare::definitelyGreaterThan(deck->d_matData.d_E, 0.);
    if (found_E)
      num_props++;

    found_K = util::compare::definitelyGreaterThan(deck->d_matData.d_K, 0.);
    if (found_K)
      num_props++;

    found_G = util::compare::definitelyGreaterThan(deck->d_matData.d_G, 0.);
    if (found_G)
      num_props++;

    found_nu = util::compare::definitelyGreaterThan(deck->d_matData.d_nu, 0.);
    if (found_nu)
      num_props++;

    if (num_props != 2) {
      std::ostringstream oss;
      oss << "Error: Require two different elastic properties for the "
                   "PdState material. Pairs supported are (E, K), (E, G), "
                   "(E, nu), (K, G).\n";
      oss << deck->printStr(0, 0) << "\n";
      std::cout << oss.str();
      exit(1);
    }

    if (util::compare::definitelyLessThan(deck->d_matData.d_Gc, 0.) &&
        util::compare::definitelyLessThan(deck->d_matData.d_KIc, 0.)) {
      std::cerr << "Error: Require either critical energy release rate Gc or "
                   "critical stress intensity factor KIc to compute the RNP "
                   "bond-based peridynamic parameters.\n";
      exit(1);
    } else if (util::compare::definitelyGreaterThan(deck->d_matData.d_Gc, 0.) &&
               util::compare::definitelyGreaterThan(deck->d_matData.d_KIc,
                                                    0.)) {
      std::cout << "Warning: Both critical energy release rate Gc and critical "
                   "stress intensity factor KIc are provided.\n";
      std::cout << "Warning: To compute the RNP bond-based peridynamic "
                   "parameters, we only require one of those.\n";
      std::cout << "Warning: Selecting critical energy release rate Gc to "
                   "compute parameters.\n";
    }

    // compute nu if not provided
    if (!found_nu) {
      if (found_E and found_G)
        deck->d_matData.d_nu =
            (0.5 * deck->d_matData.d_E / deck->d_matData.d_G) - 1.;

      if (found_E and found_K)
        deck->d_matData.d_nu =
            (3. * deck->d_matData.d_K - deck->d_matData.d_E) /
            (6. * deck->d_matData.d_K);

      if (found_G and found_K)
        deck->d_matData.d_nu =
            (3. * deck->d_matData.d_K - 2. * deck->d_matData.d_G) /
            (2. * (3. * deck->d_matData.d_K + deck->d_matData.d_G));
    }
    found_nu = true;

    // compute E if not provided
    if (!found_E) {
      if (found_K)
        deck->d_matData.d_E =
            deck->d_matData.toE(deck->d_matData.d_K, deck->d_matData.d_nu);

      if (found_G)
        deck->d_matData.d_E =
            2. * deck->d_matData.d_G * (1. + deck->d_matData.d_nu);
    }
    found_E = true;

    // compute K if not provided
    if (!found_K)
      deck->d_matData.d_K =
          deck->d_matData.toK(deck->d_matData.d_E, deck->d_matData.d_nu);

    found_K = true;


    // compute G if not provided
    if (!found_G)
      deck->d_matData.d_G =
          deck->d_matData.toGE(deck->d_matData.d_E, deck->d_matData.d_nu);

    found_G = true;

    // compute Gc (if not provided) and KIc (if not provided)
    if (deck->d_matData.d_Gc > 0.)
      deck->d_matData.d_KIc = deck->d_matData.toKIc(
          deck->d_matData.d_Gc, deck->d_matData.d_nu, deck->d_matData.d_E);

    if (deck->d_matData.d_KIc > 0. && deck->d_matData.d_Gc < 0.)
      deck->d_matData.d_Gc = deck->d_matData.toGc(
          deck->d_matData.d_KIc, deck->d_matData.d_nu, deck->d_matData.d_E);

    // compute lame parameter
    deck->d_matData.d_lambda =
        deck->d_matData.toLambdaE(deck->d_matData.d_E, deck->d_matData.d_nu);
    deck->d_matData.d_mu = deck->d_matData.d_G;

    // compute peridynamic parameters
    d_K = deck->d_matData.d_K;
    d_G = deck->d_matData.d_G;

    double d =
        (3. * d_G + std::pow(3. / 4., 4) * (d_K - 5. * d_G / 3.)) * d_horizon;
    d_s0 = std::sqrt(deck->d_matData.d_Gc / d);
  };

private:
  /*! @brief Input material deck */
  inp::MaterialDeck *d_deck_p;

  /*! @brief Horizon */
  double d_horizon;

  /*! @brief Dimension */
  size_t d_dimension;

  /*! @brief Density */
  double d_density;

  /*! @brief Base object for influence function */
  material::pd::BaseInfluenceFn *d_infFn_p;

  /**
   * @name Material parameters
   */
  /**@{*/

  /*! @brief Bulk modulus */
  double d_K;

  /*! @brief Shear modulus */
  double d_G;

  /*! @brief Critical stretch */
  double d_s0;

  /** @}*/
};

} // namespace particle

} // namespace material

#endif // MATERIAL_PD_MATERIAL_H
