////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "rnpBond.h"

#include <iostream>

#include "influenceFn.h"
#include "data/DataManager.h"
#include "inp/decks/materialDeck.h"
#include "inp/decks/modelDeck.h"
#include "inp/decks/outputDeck.h"
#include "util/compare.h"
#include "geometry/fracture.h"
#include "geometry/neighbor.h"
#include "geometry/interiorFlags.h"

material::pd::RNPBond::RNPBond(inp::MaterialDeck *deck,
                               data::DataManager *dataManager)
    : BaseMaterial(0, 0.),
      d_C(0.),
      d_beta(0.),
      d_rbar(0.),
      d_invFactor(0.),
      d_factorSc(1.),
      d_irrevBondBreak(true),
      d_contact_Kn(0.),
      d_contact_Rc(0.),
      d_deck(nullptr),
      d_dataManager_p(dataManager),
      d_baseInfluenceFn_p(nullptr) {
  d_stateActive = false;
  d_name = "RNPBond";
  d_density = deck->d_density;

  d_horizon = dataManager->getModelDeckP()->d_horizon;

  d_dimension = dataManager->getModelDeckP()->d_dim;

  //  std::cout << "RNPBond \n" << std::flush;
  //  exit(0);

  // create influence function class
  if (deck->d_influenceFnType == 0)
    d_baseInfluenceFn_p = new material::pd::ConstInfluenceFn(
        deck->d_influenceFnParams, d_dimension);
  else if (deck->d_influenceFnType == 1)
    d_baseInfluenceFn_p = new material::pd::LinearInfluenceFn(
        deck->d_influenceFnParams, d_dimension);
  else if (deck->d_influenceFnType == 2)
    d_baseInfluenceFn_p = new material::pd::GaussianInfluenceFn(
        deck->d_influenceFnParams, d_dimension);
  else {
    std::cerr << "Error: Influence function type unknown in input file.\n";
    exit(1);
  }

  double influence_moment = d_baseInfluenceFn_p->getMoment(d_dimension);

  d_irrevBondBreak = deck->d_irreversibleBondBreak;
  d_factorSc = deck->d_checkScFactor;
  if (d_dimension == 1)
    d_invFactor = std::pow(d_horizon, 2) * 2.;
  else if (d_dimension == 2)
    d_invFactor = std::pow(d_horizon, 3) * M_PI;
  else if (d_dimension == 3)
    d_invFactor = std::pow(d_horizon, 4) * 4. * M_PI / 3.;

  // check if we need to compute the material parameters
  if (deck->d_computeParamsFromElastic)
    computeParameters(deck, influence_moment);
  else {
    d_C = deck->d_bondPotentialParams[0];
    d_beta = deck->d_bondPotentialParams[1];
    d_rbar = std::sqrt(0.5 / d_beta);
    computeMaterialProperties(deck, influence_moment);
  }

  // initialize contact forces between nodes with broken bonds
  if (deck->d_applyContact) {
    // set contact radius
    d_contact_Rc = 0.9 * d_dataManager_p->getMeshP()->getMeshSize();
    // set contact force coefficient
    d_contact_Kn = deck->d_matData.d_K * 18. / (M_PI * std::pow(d_horizon, 5));

    std::cout << "RNPBonc: Mesh size: "
              << d_dataManager_p->getMeshP()->getMeshSize()
              << ", Rc: " << d_contact_Rc << ", Kn: " << d_contact_Kn << "\n";
  }

  d_deck = deck;
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

  if (deck->d_isPlaneStrain == true)
    // set Poisson's ratio to 1/4 for plain strain
    deck->d_matData.d_nu = 0.25;
  else if (deck->d_isPlaneStrain == false)
    // set Poisson's ratio to 1/3 for plain stress
    deck->d_matData.d_nu = 1. / 3.;
  else {
    std::cerr << "Error: Please specifiy the Is_Plane_Strain attribute in the "
                 "Material section!"
              << std::endl;
    std::exit(1);
  }

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
  deck->d_matData.d_nu = 0.33;

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

// std::pair<double, double> material::pd::RNPBond::getBondEF(
//    const double &r, const double &s, const double &influence, bool &fs) {
std::pair<util::Point3, double> material::pd::RNPBond::getBondEF(size_t i,
                                                                 size_t j) {
  auto force = util::Point3();
  double energy = 0.;

  // get global id of j
  auto j_id = d_dataManager_p->getNeighborP()->getNeighbor(i, j);

  // get fracture state
  auto fs = d_dataManager_p->getFractureP()->getBondState(i, j);

  // get location of nodes
  auto xi = d_dataManager_p->getMeshP()->getNode(i);
  auto ui = (*d_dataManager_p->getDisplacementP())[i];

  auto xj = d_dataManager_p->getMeshP()->getNode(j_id);
  auto uj = (*d_dataManager_p->getDisplacementP())[j_id];

  // get interior flags (to enforce no-fail region method)
  auto node_i_interior =
      d_dataManager_p->getInteriorFlagsP()->getInteriorFlag(i, xi);
  auto node_j_interior =
      d_dataManager_p->getInteriorFlagsP()->getInteriorFlag(j_id, xj);
  bool break_bonds = true;
  if (!node_i_interior || !node_j_interior) break_bonds = false;

  // get distance between nodes and bond-strain
  auto rji = xj.dist(xi);
  auto Sji = this->getS(xj - xi, uj - ui);
  auto eij = this->getBondForceDirection(xj - xi, uj - ui);

  // upper and lower bound for volume correction
  auto h = d_dataManager_p->getMeshP()->getMeshSize();
  auto check_up = d_horizon + 0.5 * h;
  auto check_low = d_horizon - 0.5 * h;

  // get corrected volume of node j
  auto voli = d_dataManager_p->getMeshP()->getNodalVolume(i);
  auto volj = d_dataManager_p->getMeshP()->getNodalVolume(j_id);
  if (util::compare::definitelyGreaterThan(rji, check_low))
    volj *= (check_up - rji) / h;

  // get influence function
  auto influence = d_baseInfluenceFn_p->getInfFn(rji / d_horizon);

  if (break_bonds) {
    // check if fracture state of the bond need to be updated
    if (d_irrevBondBreak && !fs &&
        util::compare::definitelyGreaterThan(std::abs(Sji),
                                             d_factorSc * getSc(rji)))
      fs = true;

    // update bond-state
    d_dataManager_p->getFractureP()->setBondState(i, j, fs);

    // if bond is not fractured, return energy and force from nonlinear
    // potential otherwise return energy of fractured bond, and zero force
    if (!fs) {
      energy = (influence * d_C * (1. - std::exp(-d_beta * rji * Sji * Sji)) /
                d_invFactor) *
               volj;
      force = (influence * 4. * Sji * d_C * d_beta *
               std::exp(-d_beta * rji * Sji * Sji) / d_invFactor) *
              volj * eij;
      return {force, energy};
    } else {
      // energy
      energy = influence * d_C / d_invFactor * volj;

      // normal contact force between nodes of broken bond
      auto yji = xj + uj - (xi + ui);
      auto Rji = yji.length();
      auto scalar_f = d_contact_Kn * (voli * volj / (voli + volj)) *
                      (d_contact_Rc - Rji) / Rji;
      if (scalar_f < 0.) scalar_f = 0.;
      force += -scalar_f * yji;

      return std::make_pair(force, energy);
    }
  }  // if break_bonds
  else {
    energy = (influence * d_C * d_beta * rji * Sji * Sji / d_invFactor) * volj;
    force = (influence * 4. * Sji * d_C * d_beta / d_invFactor) * volj * eij;
    return {force, energy};
  }
}

double material::pd::RNPBond::getS(const util::Point3 &dx,
                                   const util::Point3 &du) {
  return dx.dot(du) / dx.dot(dx);
}

double material::pd::RNPBond::getS(size_t i, size_t j) {
  // get location of nodes
  auto xi = d_dataManager_p->getMeshP()->getNode(i);
  auto ui = (*d_dataManager_p->getDisplacementP())[i];

  auto xj = d_dataManager_p->getMeshP()->getNode(j);
  auto uj = (*d_dataManager_p->getDisplacementP())[j];

  return this->getS(xj - xi, uj - ui);
}

double material::pd::RNPBond::getSc(const double &r) {
  return d_rbar / std::sqrt(r);
}

double material::pd::RNPBond::getSc(size_t i, size_t j) {
  auto xi = d_dataManager_p->getMeshP()->getNode(i);
  auto xj = d_dataManager_p->getMeshP()->getNode(j);

  return d_rbar / std::sqrt(xj.dist(xi));
}

util::Point3 material::pd::RNPBond::getBondForceDirection(
    const util::Point3 &dx, const util::Point3 &du) const {
  return dx / dx.length();
}

double material::pd::RNPBond::getInfFn(const double &r) const {
  return d_baseInfluenceFn_p->getInfFn(r / d_horizon);
}