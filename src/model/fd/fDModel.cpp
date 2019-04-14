// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "fDModel.h"
#include <fstream>
#include <inp/decks/restartDeck.h>

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
#include "../../inp/decks/modelDeck.h"
#include "../../inp/decks/outputDeck.h"
#include "../../inp/decks/restartDeck.h"
#include "../../inp/input.h"
#include "../../inp/policy.h"
#include "../../loading/fLoading.h"
#include "../../loading/initialCondition.h"
#include "../../loading/uLoading.h"
#include "material/pdMaterial.h"

model::FDModel::FDModel(inp::Input *deck)
    : d_massMatrix_p(nullptr), d_mesh_p(nullptr), d_fracture_p(nullptr),
      d_neighbor_p(nullptr), d_interiorFlags_p(nullptr), d_input_p(deck),
      d_policy_p(nullptr), d_initialCondition_p(nullptr), d_uLoading_p(nullptr),
      d_fLoading_p(nullptr), d_material_p(nullptr) {

  std::cout << "Here\n";

  d_modelDeck_p = deck->getModelDeck();
  d_outputDeck_p = deck->getOutputDeck();

  if (d_modelDeck_p->d_isRestartActive)
    restart(deck);
  else
    run(deck);
}

size_t model::FDModel::currentStep() { return d_n; }

float model::FDModel::getEnergy() { return d_te - d_tw + d_tk; }

void model::FDModel::run(inp::Input *deck) {

  // first initialize all the high level data
  initHObjects();

  // now initialize remaining data
  init();

  // integrate in time
  integrate();
}

void model::FDModel::restart(inp::Input *deck) {

  d_restartDeck_p = deck->getRestartDeck();

  // first initialize all the high level data
  initHObjects();

  // now initialize remaining data
  init();

  // set time step to step specified in restart deck
  d_n = d_restartDeck_p->d_step;
  d_time = double(d_n) * d_modelDeck_p->d_dt;

  // read displacement and velocities from restart file
}

void model::FDModel::initHObjects() {

  // read mesh data
  d_mesh_p = new fe::Mesh(d_input_p->getMeshDeck());

  // create neighbor list
  d_neighbor_p = new geometry::Neighbor(d_modelDeck_p->d_horizon,
                                        d_input_p->getNeighborDeck(),
                                        d_mesh_p->getNodesP());

  // create fracture data
  d_fracture_p = new geometry::Fracture(d_input_p->getFractureDeck(),
                                        d_mesh_p->getNodesP(),
                                        d_neighbor_p->getNeighborsP());

  // create interior flags
  d_interiorFlags_p = new geometry::InteriorFlags(
      d_input_p->getInteriorFlagsDeck(), d_mesh_p->getNodesP(),
      d_mesh_p->getBoundingBox());

  // initialize initial condition class
  d_initialCondition_p =
      new loading::InitialCondition(d_input_p->getInitialConditionDeck());

  // initialize loading class
  d_uLoading_p = new loading::ULoading(d_input_p->getLoadingDeck(), d_mesh_p);
  d_fLoading_p = new loading::FLoading(d_input_p->getLoadingDeck(), d_mesh_p);

  // initialize material class
  d_material_p = new material::pd::Material(d_input_p->getMaterialDeck(),
                                            d_modelDeck_p->d_dim,
                                            d_modelDeck_p->d_horizon);

  // static class
  d_policy_p = inp::Policy::getInstance(d_input_p->getPolicyDeck());
}

void model::FDModel::init() {

  d_n = 0;
  d_time = 0.;

  // get number of nodes, total number of dofs (fixed and free together)
  size_t nnodes = d_mesh_p->getNumNodes();

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
  if (d_policy_p->enablePostProcessing()) {
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
  }
}

void model::FDModel::integrate() {

  // at the beginning compute forces and apply initial and boundary condition

  // initial condition
  if (d_n == 0)
    d_initialCondition_p->apply(&d_u, &d_v, d_mesh_p);

  // boundary condition
  d_uLoading_p->apply(d_time, &d_u, &d_v, d_mesh_p);
  d_fLoading_p->apply(d_time, &d_f, d_mesh_p);

  // internal forces
  computeForces();

  // perform output at the beginning
  if (d_n == 0) {
    if (d_policy_p->enablePostProcessing())
      computePostProcFields();

    output();
  }

  // start time integration
  size_t i = d_n;
  for (i; i < d_modelDeck_p->d_Nt; i++) {
    if (d_modelDeck_p->d_timeDiscretization == "central_difference")
      integrateCD();
    else if (d_modelDeck_p->d_timeDiscretization == "velocity_verlet")
      integrateVerlet();

    if (d_n % d_outputDeck_p->d_dtOut && d_n >= d_outputDeck_p->d_dtOut) {
      if (d_policy_p->enablePostProcessing())
        computePostProcFields();

      output();
    }
  } // loop over time steps
}

void model::FDModel::integrateCD() {

  // parallel for loop
  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_mesh_p->getNumNodes(), [this](boost::uint64_t i) {
        size_t dim = this->d_mesh_p->getDimension();
        double delta_t = this->d_modelDeck_p->d_dt;
        double density = this->d_material_p->getDensity();

        // modify dofs which are not marked fixed
        if (this->d_mesh_p->isNodeFree(i, 0)) {

          double u_old = this->d_u[i].d_x;

          this->d_u[i].d_x += delta_t * delta_t * this->d_f[i].d_x / density +
                              delta_t * this->d_v[i].d_x;

          this->d_v[i].d_x = (this->d_u[i].d_x - u_old) / delta_t;
        }

        if (dim > 1)
          if (this->d_mesh_p->isNodeFree(i, 1)) {

            double u_old = this->d_u[i].d_y;

            this->d_u[i].d_y += delta_t * delta_t * this->d_f[i].d_y / density +
                                delta_t * this->d_v[i].d_y;

            this->d_v[i].d_y = (this->d_u[i].d_y - u_old) / delta_t;
          }

        if (dim > 2)
          if (this->d_mesh_p->isNodeFree(i, 2)) {

            double u_old = this->d_u[i].d_z;

            this->d_u[i].d_z += delta_t * delta_t * this->d_f[i].d_z / density +
                                delta_t * this->d_v[i].d_z;

            this->d_v[i].d_z = (this->d_u[i].d_z - u_old) / delta_t;
          }

        // reset force
        this->d_f[i].d_x = 0.;
        this->d_f[i].d_y = 0.;
        this->d_f[i].d_z = 0.;
      }); // end of parallel for loop

  f.get();

  // compute forces and energy due to new displacement field (this will be
  // used in next time step)
  d_n++;
  d_time += d_modelDeck_p->d_dt;

  // boundary condition
  d_uLoading_p->apply(d_time, &d_u, &d_v, d_mesh_p);
  d_fLoading_p->apply(d_time, &d_f, d_mesh_p);

  // internal forces
  computeForces();
}

void model::FDModel::integrateVerlet() {

  // step 1 and 2 : Compute v_mid and u_new
  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_mesh_p->getNumNodes(), [this](boost::uint64_t i) {
        size_t dim = this->d_mesh_p->getDimension();
        double delta_t = this->d_modelDeck_p->d_dt;
        double density = this->d_material_p->getDensity();

        // modify dofs which are not marked fixed
        if (this->d_mesh_p->isNodeFree(i, 0)) {

          this->d_v[i].d_x += 0.5 * delta_t * this->d_f[i].d_x / density;

          this->d_u[i].d_x += delta_t * this->d_v[i].d_x;
        }

        if (dim > 1)
          if (this->d_mesh_p->isNodeFree(i, 1)) {

            this->d_v[i].d_y += 0.5 * delta_t * this->d_f[i].d_y / density;

            this->d_u[i].d_y += delta_t * this->d_v[i].d_y;
          }

        if (dim > 2)
          if (this->d_mesh_p->isNodeFree(i, 2)) {

            this->d_v[i].d_z += 0.5 * delta_t * this->d_f[i].d_z / density;

            this->d_u[i].d_z += delta_t * this->d_v[i].d_z;
          }

        // reset force
        this->d_f[i].d_x = 0.;
        this->d_f[i].d_y = 0.;
        this->d_f[i].d_z = 0.;
      }); // end of parallel for loop

  f.get();

  // compute forces and energy due to new displacement field (this will be
  // used in next time step)
  d_n++;
  d_time += d_modelDeck_p->d_dt;

  // boundary condition
  d_uLoading_p->apply(d_time, &d_u, &d_v, d_mesh_p);
  d_fLoading_p->apply(d_time, &d_f, d_mesh_p);

  // internal forces
  computeForces();

  // Step 3: Compute v_new
  f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_mesh_p->getNumNodes(), [this](boost::uint64_t i) {
        size_t dim = this->d_mesh_p->getDimension();
        double delta_t = this->d_modelDeck_p->d_dt;
        double density = this->d_material_p->getDensity();

        // modify dofs which are not marked fixed
        if (this->d_mesh_p->isNodeFree(i, 0)) {

          this->d_v[i].d_x += 0.5 * delta_t * this->d_f[i].d_x / density;

          this->d_u[i].d_x += delta_t * this->d_v[i].d_x;
        }

        if (dim > 1)
          if (this->d_mesh_p->isNodeFree(i, 1)) {

            this->d_v[i].d_y += 0.5 * delta_t * this->d_f[i].d_y / density;

            this->d_u[i].d_y += delta_t * this->d_v[i].d_y;
          }

        if (dim > 2)
          if (this->d_mesh_p->isNodeFree(i, 2)) {

            this->d_v[i].d_z += 0.5 * delta_t * this->d_f[i].d_z / density;

            this->d_u[i].d_z += delta_t * this->d_v[i].d_z;
          }
      }); // end of parallel for loop

  f.get();
}

void model::FDModel::output() {}

void model::FDModel::computeForces() {

  if (d_material_p->isStateActive())
    computeHydrostaticStrains();

  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_mesh_p->getNumNodes(),
      [this](boost::uint64_t i) {
        // local variable to hold force
        util::Point3 force_i = util::Point3();

        // reference coordinate and displacement at the node
        util::Point3 xi = this->d_mesh_p->getNode(i);
        util::Point3 ui = this->d_u[i];

        // get hydrostatic energy and force
        std::pair<double, double> gi;
        if (this->d_material_p->isStateActive())
          gi = d_material_p->getStateEF(this->d_hS[i]);

        // get interior flag
        auto node_i_interior = this->d_interiorFlags_p->getInteriorFlag(i, xi);

        // upper and lower bound for volume correction
        double h = d_mesh_p->getMeshSize();
        double check_up = d_modelDeck_p->d_horizon + 0.5 * h;
        double check_low = d_modelDeck_p->d_horizon - 0.5 * h;

        // inner loop over neighbors
        auto i_neighs = this->d_neighbor_p->getNeighbors(i);
        for (size_t j = 0; j < i_neighs.size(); j++) {

          size_t j_id = i_neighs[j];

          // there are two contributions to force at node i
          // 1. From bond j-i due to bond-based forces
          // 2. From hydrostatic strains at node j and node i due to
          // hydrostatic forces

          // compute bond-based contribution
          util::Point3 xj = this->d_mesh_p->getNode(j_id);
          util::Point3 uj = this->d_u[j_id];
          util::Point3 xji = xj - xi;
          double rji = xji.length();
          double Sji = this->d_material_p->getS(xji, uj - ui);

          // get corrected volume of node j
          double volj = this->d_mesh_p->getNodalVolume(j_id);
          if (util::compare::definitelyGreaterThan(rji, check_low))
            volj *= (check_up - rji) / h;

          // check if bond is in no-fail region
          bool break_bonds = true;
          if (!node_i_interior ||
              !this->d_interiorFlags_p->getInteriorFlag(j_id, xj))
            break_bonds = false;

          // get peridynamics force and energy density between bond i and j
          bool fs = this->d_fracture_p->getBondState(i, j);
          std::pair<double, double> ef =
              this->d_material_p->getBondEF(rji, Sji, fs, break_bonds);

          // update the fractured state of bond
          this->d_fracture_p->setBondState(i, j, fs);

          // compute the contribution of bond force to force at i
          double scalar_f = ef.second * volj / rji;

          force_i.d_x += scalar_f * xji.d_x;
          force_i.d_y += scalar_f * xji.d_y;
          force_i.d_z += scalar_f * xji.d_z;

          // compute state-based contribution
          if (d_material_p->isStateActive() &&
              d_material_p->addBondContribToState(Sji, rji)) {

            // Compute gj while noting that gi is already computed
            std::pair<double, double> gj =
                d_material_p->getStateEF(this->d_hS[j_id]);

            double scalar_g =
                d_material_p->getStateForce(gi.second + gj.second, rji) / rji;

            force_i.d_x += scalar_g * xji.d_x;
            force_i.d_y += scalar_g * xji.d_y;
            force_i.d_z += scalar_g * xji.d_z;
          }
        } // loop over neighboring nodes

        // update force and energy
        this->d_f[i] = force_i;
      } // loop over nodes

  ); // end of parallel for loop

  f.get();
}

void model::FDModel::computeHydrostaticStrains() {
  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_mesh_p->getNumNodes(),
      [this](boost::uint64_t i) {

        // local variable to hold strain
        double hydro_strain_i = 0.;

        // reference coordinate and displacement at the node
        util::Point3 xi = this->d_mesh_p->getNode(i);
        util::Point3 ui = this->d_u[i];

        // upper and lower bound for volume correction
        double h = d_mesh_p->getMeshSize();
        double check_up = d_modelDeck_p->d_horizon + 0.5 * h;
        double check_low = d_modelDeck_p->d_horizon - 0.5 * h;

        // inner loop over neighbors
        auto i_neighs = this->d_neighbor_p->getNeighbors(i);
        for (size_t j = 0; j < i_neighs.size(); j++) {

          size_t j_id = i_neighs[j];

          // compute bond-based contribution
          util::Point3 uji = this->d_u[j_id] - ui;
          util::Point3 xji = this->d_mesh_p->getNode(j_id) - xi;
          double rji = xji.length();
          double Sji = this->d_material_p->getS(xji, uji);

          // get corrected volume of node j
          double volj = this->d_mesh_p->getNodalVolume(j_id);
          if (util::compare::definitelyGreaterThan(rji, check_low))
            volj *= (check_up - rji) / h;

          // compute the contribution of bond to hydrostatic strain at i
          if (d_material_p->addBondContribToState(Sji, rji)) {

            hydro_strain_i += volj * d_material_p->getBondContribToHydroStrain
                (Sji, rji);
          }
        } // loop over neighboring nodes

        // update
        this->d_hS[i] = hydro_strain_i;
      } // loop over nodes

  ); // end of parallel for loop

  f.get();
}

void model::FDModel::computePostProcFields() {

  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_mesh_p->getNumNodes(),
      [this](boost::uint64_t i) {
        // local variable to hold force
        util::Point3 force_i = util::Point3();

        // reference coordinate and displacement at the node
        util::Point3 xi = this->d_mesh_p->getNode(i);
        util::Point3 ui = this->d_u[i];

        // get hydrostatic energy and force
        std::pair<double, double> gi;
        if (this->d_material_p->isStateActive())
          gi = d_material_p->getStateEF(this->d_hS[i]);

        // get interior flag
        auto node_i_interior = this->d_interiorFlags_p->getInteriorFlag(i, xi);

        // upper and lower bound for volume correction
        double h = d_mesh_p->getMeshSize();
        double check_up = d_modelDeck_p->d_horizon + 0.5 * h;
        double check_low = d_modelDeck_p->d_horizon - 0.5 * h;

        // inner loop over neighbors
        auto i_neighs = this->d_neighbor_p->getNeighbors(i);
        for (size_t j = 0; j < i_neighs.size(); j++) {

          size_t j_id = i_neighs[j];

          // there are two contributions to force at node i
          // 1. From bond j-i due to bond-based forces
          // 2. From hydrostatic strains at node j and node i due to
          // hydrostatic forces

          // compute bond-based contribution
          util::Point3 xj = this->d_mesh_p->getNode(j_id);
          util::Point3 uj = this->d_u[j_id];
          util::Point3 xji = xj - xi;
          double rji = xji.length();
          double Sji = this->d_material_p->getS(xji, uj - ui);

          // get corrected volume of node j
          double volj = this->d_mesh_p->getNodalVolume(j_id);
          if (util::compare::definitelyGreaterThan(rji, check_low))
            volj *= (check_up - rji) / h;

          // check if bond is in no-fail region
          bool break_bonds = true;
          if (!node_i_interior ||
              !this->d_interiorFlags_p->getInteriorFlag(j_id, xj))
            break_bonds = false;

          // get peridynamics force and energy density between bond i and j
          bool fs = this->d_fracture_p->getBondState(i, j);
          std::pair<double, double> ef =
              this->d_material_p->getBondEF(rji, Sji, fs, break_bonds);

          // update the fractured state of bond
          this->d_fracture_p->setBondState(i, j, fs);

          // compute the contribution of bond force to force at i
          double scalar_f = ef.second * volj / rji;

          force_i.d_x += scalar_f * xji.d_x;
          force_i.d_y += scalar_f * xji.d_y;
          force_i.d_z += scalar_f * xji.d_z;

          // compute state-based contribution
          if (d_material_p->isStateActive() &&
              d_material_p->addBondContribToState(Sji, rji)) {

            // Compute gj while noting that gi is already computed
            std::pair<double, double> gj =
                d_material_p->getStateEF(this->d_hS[j_id]);

            double scalar_g =
                d_material_p->getStateForce(gi.second + gj.second, rji) / rji;

            force_i.d_x += scalar_g * xji.d_x;
            force_i.d_y += scalar_g * xji.d_y;
            force_i.d_z += scalar_g * xji.d_z;
          }
        } // loop over neighboring nodes

        // update force and energy
        this->d_f[i] = force_i;
      } // loop over nodes

  ); // end of parallel for loop

  f.get();
}