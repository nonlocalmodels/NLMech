////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "rnpBond.h"

#include <iostream>

#include "inp/decks/materialDeck.h"
#include "util/compare.h"

material::pd::RNPBond::RNPBond(inp::MaterialDeck *deck, const size_t &dim,
                               const double &horizon, const double &M)
    : BaseMaterial(dim, horizon),
      d_C(0.),
      d_beta(0.),
      d_rbar(0.),
      d_invFactor(0.),
      d_factorSc(1.),
      d_irrevBondBreak(true) {
  d_irrevBondBreak = deck->d_irreversibleBondBreak;
  d_factorSc = deck->d_checkScFactor;
  if (dim == 1)
    d_invFactor = std::pow(horizon, 2) * 2.;
  else if (dim == 2)
    d_invFactor = std::pow(horizon, 3) * M_PI;
  else if (dim == 3)
    d_invFactor = std::pow(horizon, 4) * 4. * M_PI / 3.;

  // check if we need to compute the material parameters
  if (deck->d_computeParamsFromElastic)
    computeParameters(deck, M);
  else {
    d_C = deck->d_bondPotentialParams[0];
    d_beta = deck->d_bondPotentialParams[1];
    d_rbar = std::sqrt(0.5 / d_beta);
    computeMaterialProperties(deck, M);
  }
}

void material::pd::RNPBond::computeParameters(inp::MaterialDeck *deck,
                                              const double &M) {
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
    std::cout << "Warning: Selecting Young's modulus to compute parameters.\n";
  }

  if (util::compare::definitelyLessThan(deck->d_matData.d_Gc, 0.) &&
      util::compare::definitelyLessThan(deck->d_matData.d_KIc, 0.)) {
    std::cerr << "Error: Require either critical energy release rate Gc or "
                 "critical stress intensity factor KIc to compute the RNP "
                 "bond-based peridynamic parameters.\n";
    exit(1);
  } else if (util::compare::definitelyGreaterThan(deck->d_matData.d_Gc, 0.) &&
             util::compare::definitelyGreaterThan(deck->d_matData.d_KIc, 0.)) {
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
}

void material::pd::RNPBond::computeMaterialProperties(inp::MaterialDeck *deck,
                                                      const double &M) {
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
  deck->d_matData.d_E =
      deck->d_matData.toELambda(deck->d_matData.d_lambda, deck->d_matData.d_nu);
  deck->d_matData.d_K =
      deck->d_matData.toK(deck->d_matData.d_E, deck->d_matData.d_nu);
  deck->d_matData.d_KIc = deck->d_matData.toKIc(
      deck->d_matData.d_Gc, deck->d_matData.d_nu, deck->d_matData.d_E);
}

/*
std::pair<double, double>
material::pd::RNPBond::getBondEF(const double &r, const double &s,
                                 const double &influence, bool &fs) {

  // check if fracture state of the bond need to be updated
  if (d_irrevBondBreak && !fs &&
      util::compare::definitelyGreaterThan(std::abs(s), d_factorSc * getSc(r)))
    fs = true;

  // if bond is not fractured, return energy and force from nonlinear potential
  // otherwise return energy of fractured bond, and zero force
  if (!fs)
    return std::make_pair(
        influence * d_C * (1. - std::exp(-d_beta * r * s * s)) / d_invFactor,
        influence * 4. * s * d_C * d_beta * std::exp(-d_beta * r * s * s) /
            d_invFactor);
  else
    return std::make_pair(d_C / d_invFactor, 0.);
}
*/

std::pair<double, double> material::pd::RNPBond::getBondEFNoFail(
    const double &r, const double &s, const double &influence) {
  // Return force and energy for no-fail region bonds
  return std::make_pair(influence * d_C * d_beta * r * s * s / d_invFactor,
                        influence * 4. * s * d_C * d_beta / d_invFactor);
}

double material::pd::RNPBond::getSc(const double &r) {
  return d_rbar / std::sqrt(r);
}
