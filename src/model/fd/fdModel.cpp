// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "fdModel.h"
#include <iostream>

// include high level class declarations
#include "../../geometry/fracture.h"
#include "../../geometry/geometry.h"
#include "../../geometry/interiorFlags.h"
#include "../../geometry/neighbor.h"
#include "../../io/input.h"
#include "../../io/policy.h"
#include "../../loading/initialCondition.h"
#include "../../loading/loading.h"
#include "../../material/material.h"

model::FDModel::FDModel(io::Input *deck)  {

  std::cout<<"Here\n";

  // store pointer to the input data
  d_input_p = deck;

  // first initialize all the high level data
  initHObjects();

  // now initialize remaining data
  init();

  // setup and apply initial displacement and/or velocity and/or force boundary
  // condition
  setupBoundaryCondition();
}

size_t model::FDModel::currentStep() {
  return d_n;
}

float model::FDModel::getEnergy() {
  return d_te - d_tw + d_tk;
}

void model::FDModel::initHObjects() {
  d_geometry_p = new geometry::Geometry(d_input_p->getGeometryDeck());
  d_fracture_p = new geometry::Fracture(d_input_p->getFractureDeck());
  d_interiorFlags_p =
      new geometry::InteriorFlags(d_input_p->getInteriorFlagsDeck());
  d_neighbor_p = new geometry::Neighbor(d_input_p->getNeighborDeck());
  d_policy_p = new io::Policy(d_input_p->getPolicyDeck());
  d_initialCondition_p =
      new loading::InitialCondition(d_input_p->getInitialConditionDeck());
  d_loading_p = new loading::Loading(d_input_p->getLoadingDeck());
  d_material_p = new material::Material(d_input_p->getMaterialDeck());
}

void model::FDModel::init() {

  d_n = 0;

  // get number of nodes, total number of dofs (fixed and free together)
//  size_t nnodes = d_geometry_p->getNumNodes();
//  size_t ndofs = d_geometry_p->getNumDofs();
//  size_t dim = d_geometry_p->getDimension();
  size_t nnodes = 1;
  size_t ndofs = 1;
  size_t dim = 1;

  // initialize major simulation data
//  d_y = d_geometry_p->getNodes();
  d_v = std::vector<double>(ndofs, 0.0);
  d_f = std::vector<double>(ndofs, 0.0);
//  if (d_policy_p->populateData("FDModel_d_hS"))
//    d_hS = std::vector<double>(nnodes, 0.0);

  // initialize minor simulation data
//  if (d_policy_p->populateData("FDModel_d_e"))
//    d_e = std::vector<float>(nnodes, 0.0);
//  if (d_policy_p->populateData("FDModel_d_w"))
//    d_w = std::vector<float>(nnodes, 0.0);
//  if (d_policy_p->populateData("FDModel_d_phi"))
//    d_phi = std::vector<float>(nnodes, 0.0);
//  if (d_policy_p->populateData("FDModel_d_Z"))
//    d_Z = std::vector<float>(nnodes, 0.0);
//  if (d_policy_p->populateData("FDModel_d_eF"))
//    d_eF = std::vector<float>(nnodes, 0.0);
//  if (d_policy_p->populateData("FDModel_d_eFB"))
//    d_eFB = std::vector<float>(nnodes, 0.0);
  d_te = 0.0;
  d_tw = 0.0;
  d_tk = 0.0;
  d_teF = 0.0;
  d_teFB = 0.0;
}

void model::FDModel::integrate() {}

void model::FDModel::integrateCD() {}

void model::FDModel::integrateVerlet() {}

void model::FDModel::setupBoundaryCondition() {}

void model::FDModel::applyDisplacementBC() {}

void model::FDModel::applyForceBC() {}

void model::FDModel::applyInitialCondition() {}

void model::FDModel::output() {}

void model::FDModel::debug(const float e_old) {}