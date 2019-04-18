// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "pdMaterial.h"
#include "inp/decks/materialDeck.h"
#include "pd/influenceFn.h"
#include "pd/rnpBond.h"
#include "util/compare.h"
#include <iostream>

material::pd::Material::Material(inp::MaterialDeck *deck, const size_t &dim,
                                 const double &horizon)
    : d_deck_p(deck), d_stateActive(false), d_horizon(horizon),
      d_dimension(dim), d_density(deck->d_density), d_baseMaterial_p(nullptr),
      d_baseInfluenceFn_p(nullptr) {

  // create influence function class
  if (d_deck_p->d_influenceFnType == 0)
    d_baseInfluenceFn_p = new material::pd::ConstInfluenceFn(
        d_deck_p->d_influenceFnParams, d_dimension);
  else if (d_deck_p->d_influenceFnType == 1)
    d_baseInfluenceFn_p = new material::pd::LinearInfluenceFn(
        d_deck_p->d_influenceFnParams, d_dimension);
  else if (d_deck_p->d_influenceFnType == 2)
    d_baseInfluenceFn_p = new material::pd::GaussianInfluenceFn(
        d_deck_p->d_influenceFnParams, d_dimension);
  else {
    std::cerr << "Error: Influence function type unknown in input file.\n";
    exit(1);
  }

  // create material class
  if (d_deck_p->d_materialType == "PDBond") {
    d_stateActive = false;
    d_baseMaterial_p =
        new material::pd::RNPBond(d_deck_p, d_dimension, d_horizon,
                                  d_baseInfluenceFn_p->getMoment(d_dimension));
  } else {
    std::cerr << "Error: Material type = " << d_deck_p->d_materialType
              << ". Currently only PDBond is implemented.\n";
    exit(1);
  }
}

std::pair<double, double>
material::pd::Material::getBondEF(const double &r, const double &s, bool &fs,
                                  const bool &break_bonds) {
  if (break_bonds)
    return d_baseMaterial_p->getBondEF(
        r, s, d_baseInfluenceFn_p->getInfFn(r / d_horizon), fs);
  else
    return d_baseMaterial_p->getBondEFNoFail(
        r, s, d_baseInfluenceFn_p->getInfFn(r / d_horizon));
}

std::pair<double, double>
material::pd::Material::getStateEF(const double &theta) {
  return d_baseMaterial_p->getStateEF(theta);
}

double material::pd::Material::getInfFn(const double &r) {
  return d_baseInfluenceFn_p->getInfFn(r / d_horizon);
}

double material::pd::Material::getMoment(const size_t &i) {
  return d_baseInfluenceFn_p->getMoment(i);
}

double material::pd::Material::getS(const util::Point3 &dx,
                                    const util::Point3 &du) {
  return dx.dot(du) / dx.dot(dx);
}

double material::pd::Material::getDensity() { return d_density; }

bool material::pd::Material::isStateActive() { return d_stateActive; }

bool material::pd::Material::addBondContribToState(const double &S,
                                                   const double &r) {
  return d_deck_p->d_stateContributionFromBrokenBond ||
         util::compare::definitelyLessThan(S, d_deck_p->d_checkScFactor *
                                                  d_baseMaterial_p->getSc(r));
}

double material::pd::Material::getSc(const double &r) {
  return d_baseMaterial_p->getSc(r);
}

double material::pd::Material::getStateForce(const double &g_prime,
                                             const double &r) {
  return d_baseMaterial_p->getStateForce(g_prime, r);
}

double material::pd::Material::getBondContribToHydroStrain(const double &S,
                                                           const double &r) {
  return d_baseMaterial_p->getBondContribToHydroStrain(S, r);
}
