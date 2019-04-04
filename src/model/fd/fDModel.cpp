// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "fDModel.h"
#include <fstream>

// utils
#include "../../util/matrix.h"
#include "../../util/point.h"

// include high level class declarations
#include "../../fe/massMatrix.h"
#include "../../fe/mesh.h"
#include "../../fe/quadrature.h"
#include "../../geometry/fracture.h"
#include "../../geometry/interiorFlags.h"
#include "../../geometry/neighbor.h"
#include "../../inp/input.h"
#include "../../inp/policy.h"
#include "../../loading/initialCondition.h"
#include "../../loading/loading.h"
#include "../../material/material.h"

model::FDModel::FDModel(inp::Input *deck)
    : d_massMatrix_p(nullptr), d_mesh_p(nullptr), d_fracture_p(nullptr),
      d_neighbor_p(nullptr), d_interiorFlags_p(nullptr), d_input_p(deck),
      d_policy_p(nullptr), d_initialCondition_p(nullptr), d_loading_p(nullptr),
      d_material_p(nullptr) {

  std::cout << "Here\n";

  // first initialize all the high level data
  initHObjects();

  // now initialize remaining data
  init();

  // setup and apply initial displacement and/or velocity and/or force boundary
  // condition
  setupBoundaryCondition();
}

size_t model::FDModel::currentStep() { return d_n; }

float model::FDModel::getEnergy() { return d_te - d_tw + d_tk; }

void model::FDModel::initHObjects() {

  d_mesh_p = new fe::Mesh(d_input_p->getMeshDeck());

  //  d_quadrature_p = new fe::Quadrature(d_input_p->getQuadratureDeck(),
  //      d_mesh_p->getElementType());

  d_massMatrix_p = new fe::MassMatrix(d_input_p->getMassMatrixDeck());

  d_fracture_p = new geometry::Fracture(d_input_p->getFractureDeck());
  d_interiorFlags_p =
      new geometry::InteriorFlags(d_input_p->getInteriorFlagsDeck());
  d_neighbor_p = new geometry::Neighbor(d_input_p->getNeighborDeck());
  d_initialCondition_p =
      new loading::InitialCondition(d_input_p->getInitialConditionDeck());
  d_loading_p = new loading::Loading(d_input_p->getLoadingDeck());
  d_material_p = new material::Material(d_input_p->getMaterialDeck());

  // static class
  d_policy_p = inp::Policy::getInstance(d_input_p->getPolicyDeck());
}

void model::FDModel::init() {

  d_n = 0;

  // get number of nodes, total number of dofs (fixed and free together)
  size_t nnodes = d_mesh_p->getNumNodes();
  size_t ndofs = d_mesh_p->getNumDofs();
  size_t dim = d_mesh_p->getDimension();

  return;

  // initialize major simulation data
  d_u = std::vector<util::Point3>(nnodes, util::Point3());
  d_v = std::vector<util::Point3>(nnodes, util::Point3());
  d_f = std::vector<util::Point3>(nnodes, util::Point3());

  // check material type and if needed update policy
  if (!d_material_p->isStateActive()) {
    d_policy_p->addToTags(0, "Model_d_hS");
    d_policy_p->addToTags(0, "Model_d_eF");
  }

  // initialize minor simulation data
  if (d_policy_p->populateData("Model_d_e"))
    d_e = std::vector<float>(nnodes, 0.0);
  if (d_policy_p->populateData("Model_d_w"))
    d_w = std::vector<float>(nnodes, 0.0);
  if (d_policy_p->populateData("Model_d_phi"))
    d_phi = std::vector<float>(nnodes, 0.0);
  if (d_policy_p->populateData("Model_d_Z"))
    d_Z = std::vector<float>(nnodes, 0.0);
  if (d_policy_p->populateData("Model_d_eF"))
    d_eF = std::vector<float>(nnodes, 0.0);
  if (d_policy_p->populateData("Model_d_eFB"))
    d_eFB = std::vector<float>(nnodes, 0.0);
  if (d_policy_p->populateData("Model_d_hS"))
    d_hS = std::vector<double>(nnodes, 0.0);

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

void model::FDModel::debug(float e_old) {}