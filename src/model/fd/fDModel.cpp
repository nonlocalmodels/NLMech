// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "fDModel.h"

// utils
#include "rw/reader.h"
#include "rw/writer.h"
#include "util/fastMethods.h"
#include "util/matrix.h"
#include "util/point.h"
#include "util/compare.h"

// include high level class declarations
#include "fe/massMatrix.h"
#include "fe/mesh.h"
#include "geometry/fracture.h"
#include "geometry/interiorFlags.h"
#include "geometry/neighbor.h"
#include "inp/decks/modelDeck.h"
#include "inp/decks/outputDeck.h"
#include "inp/decks/restartDeck.h"
#include "inp/input.h"
#include "inp/policy.h"
#include "loading/fLoading.h"
#include "loading/initialCondition.h"
#include "loading/uLoading.h"
#include "material/pdMaterial.h"

// standard lib
#include <fstream>

model::FDModel::FDModel(inp::Input *deck)
    : d_massMatrix_p(nullptr), d_mesh_p(nullptr), d_fracture_p(nullptr),
      d_neighbor_p(nullptr), d_interiorFlags_p(nullptr), d_input_p(deck),
      d_policy_p(nullptr), d_initialCondition_p(nullptr), d_uLoading_p(nullptr),
      d_fLoading_p(nullptr), d_material_p(nullptr) {

  d_modelDeck_p = deck->getModelDeck();
  d_outputDeck_p = deck->getOutputDeck();
  d_policy_p = inp::Policy::getInstance(d_input_p->getPolicyDeck());

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

  // read displacement and velocity from restart file
  rw::reader::readVtuFileRestart(d_restartDeck_p->d_file, &d_u, &d_v,
      d_mesh_p->getNodesP());

  // integrate in time
  integrate();
}

void model::FDModel::initHObjects() {
  std::cout << "FDModel: Initializing high level objects.\n";
  // read mesh data
  std::cout << "FDModel: Creating mesh.\n";
  d_mesh_p = new fe::Mesh(d_input_p->getMeshDeck());
  d_mesh_p->clearElementData();

  // create neighbor list
  std::cout << "FDModel: Creating neighbor list.\n";
  d_neighbor_p = new geometry::Neighbor(d_modelDeck_p->d_horizon,
                                        d_input_p->getNeighborDeck(),
                                        d_mesh_p->getNodesP());

  // create fracture data
  std::cout << "FDModel: Creating edge crack if any and modifying the "
               "fracture state of bonds.\n";
  d_fracture_p = new geometry::Fracture(d_input_p->getFractureDeck(),
                                        d_mesh_p->getNodesP(),
                                        d_neighbor_p->getNeighborsP());

  // create interior flags
  std::cout << "FDModel: Creating interior flags for nodes.\n";
  d_interiorFlags_p = new geometry::InteriorFlags(
      d_input_p->getInteriorFlagsDeck(), d_mesh_p->getNodesP(),
      d_mesh_p->getBoundingBox());

  // initialize initial condition class
  std::cout << "FDModel: Initializing initial condition object.\n";
  d_initialCondition_p =
      new loading::InitialCondition(d_input_p->getInitialConditionDeck());

  // initialize loading class
  std::cout << "FDModel: Initializing displacement loading object.\n";
  d_uLoading_p = new loading::ULoading(d_input_p->getLoadingDeck(), d_mesh_p);
  std::cout << "FDModel: Initializing force loading object.\n";
  d_fLoading_p = new loading::FLoading(d_input_p->getLoadingDeck(), d_mesh_p);

  // initialize material class
  std::cout << "FDModel: Initializing material object.\n";
  d_material_p = new material::pd::Material(d_input_p->getMaterialDeck(),
                                            d_modelDeck_p->d_dim,
                                            d_modelDeck_p->d_horizon);
}

void model::FDModel::init() {
  std::cout << "FDModel: Initializing basic datas.\n";

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
  if (d_policy_p->populateData("Model_d_hS"))
    d_hS = std::vector<double>(nnodes, 0.0);

  // initialize minor simulation data
  if (d_policy_p->enablePostProcessing()) {

    std::string tag = "Strain_Energy";
    if (d_outputDeck_p->isTagInOutput(tag)) {
      // this data is asked in output file
      // but check if policy allows its population
      if (d_policy_p->populateData("Model_d_e"))
        d_e = std::vector<float>(nnodes, 0.0);
    }
    else {
      // this data is not asked in output thus we disable it
      d_policy_p->addToTags(0, "Model_d_e");
    }

    tag = "Work_Done";
    if (d_outputDeck_p->isTagInOutput(tag)) {
      if (d_policy_p->populateData("Model_d_w"))
        d_w = std::vector<float>(nnodes, 0.0);
    }
    else
      d_policy_p->addToTags(0, "Model_d_w");

    tag = "Damage_Phi";
    if (d_outputDeck_p->isTagInOutput(tag)) {
      if (d_policy_p->populateData("Model_d_phi"))
        d_phi = std::vector<float>(nnodes, 0.0);
    }
    else
      d_policy_p->addToTags(0, "Model_d_phi");

    tag = "Damage_Z";
    if (d_outputDeck_p->isTagInOutput(tag)) {
      if (d_policy_p->populateData("Model_d_Z"))
        d_Z = std::vector<float>(nnodes, 0.0);
    }
    else
      d_policy_p->addToTags(0, "Model_d_Z");

    tag = "Fracture_Perienergy_Total";
    if (d_outputDeck_p->isTagInOutput(tag)) {
      if (d_policy_p->populateData("Model_d_eF"))
        d_eF = std::vector<float>(nnodes, 0.0);
    }
    else
      d_policy_p->addToTags(0, "Model_d_eF");

    tag = "Fracture_Perienergy_Bond";
    if (d_outputDeck_p->isTagInOutput(tag)) {
      if (d_policy_p->populateData("Model_d_eFB"))
        d_eFB = std::vector<float>(nnodes, 0.0);
    }
    else
      d_policy_p->addToTags(0, "Model_d_eFB");
  }
}

void model::FDModel::integrate() {

  // apply initial loading
  if (d_n == 0)
    d_initialCondition_p->apply(&d_u, &d_v, d_mesh_p);

  // apply loading
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

    // handle general output
    if ((d_n % d_outputDeck_p->d_dtOut == 0) &&
        (d_n >= d_outputDeck_p->d_dtOut)) {

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
        auto dim = this->d_mesh_p->getDimension();
        auto delta_t = this->d_modelDeck_p->d_dt;
        auto fact = delta_t * delta_t / this->d_material_p->getDensity();

        if (this->d_mesh_p->isNodeFree(i, 0)) {

          auto u_old = this->d_u[i].d_x;

          this->d_u[i].d_x +=
              fact * this->d_f[i].d_x + delta_t * this->d_v[i].d_x;

          this->d_v[i].d_x = (this->d_u[i].d_x - u_old) / delta_t;
        }

        if (dim > 1)
          if (this->d_mesh_p->isNodeFree(i, 1)) {

            auto u_old = this->d_u[i].d_y;

            this->d_u[i].d_y +=
                fact * this->d_f[i].d_y + delta_t * this->d_v[i].d_y;

            this->d_v[i].d_y = (this->d_u[i].d_y - u_old) / delta_t;
          }

        if (dim > 2)
          if (this->d_mesh_p->isNodeFree(i, 2)) {

            auto u_old = this->d_u[i].d_z;

            this->d_u[i].d_z +=
                fact * this->d_f[i].d_z + delta_t * this->d_v[i].d_z;

            this->d_v[i].d_z = (this->d_u[i].d_z - u_old) / delta_t;
          }

        // reset force
        this->d_f[i] = util::Point3();
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
        auto dim = this->d_mesh_p->getDimension();
        auto delta_t = this->d_modelDeck_p->d_dt;
        auto fact = 0.5 * delta_t / this->d_material_p->getDensity();

        // modify dofs which are not marked fixed
        if (this->d_mesh_p->isNodeFree(i, 0)) {
          this->d_v[i].d_x += fact * this->d_f[i].d_x;
          this->d_u[i].d_x += delta_t * this->d_v[i].d_x;
        }

        if (dim > 1)
          if (this->d_mesh_p->isNodeFree(i, 1)) {
            this->d_v[i].d_y += fact * this->d_f[i].d_y;
            this->d_u[i].d_y += delta_t * this->d_v[i].d_y;
          }

        if (dim > 2)
          if (this->d_mesh_p->isNodeFree(i, 2)) {
            this->d_v[i].d_z += fact * this->d_f[i].d_z;
            this->d_u[i].d_z += delta_t * this->d_v[i].d_z;
          }

        // reset force
        this->d_f[i] = util::Point3();
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
        auto dim = this->d_mesh_p->getDimension();
        auto fact =
            0.5 * this->d_modelDeck_p->d_dt / this->d_material_p->getDensity();

        // modify dofs which are not marked fixed
        if (this->d_mesh_p->isNodeFree(i, 0))
          this->d_v[i].d_x += fact * this->d_f[i].d_x;

        if (dim > 1)
          if (this->d_mesh_p->isNodeFree(i, 1))
            this->d_v[i].d_y += fact * this->d_f[i].d_y;

        if (dim > 2)
          if (this->d_mesh_p->isNodeFree(i, 2))
            this->d_v[i].d_z += fact * this->d_f[i].d_z;
      }); // end of parallel for loop

  f.get();
}

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
        auto xi = this->d_mesh_p->getNode(i);
        auto ui = this->d_u[i];

        // get interior flag
        auto node_i_interior = this->d_interiorFlags_p->getInteriorFlag(i, xi);

        // upper and lower bound for volume correction
        auto h = d_mesh_p->getMeshSize();
        auto check_up = d_modelDeck_p->d_horizon + 0.5 * h;
        auto check_low = d_modelDeck_p->d_horizon - 0.5 * h;

        // inner loop over neighbors
        auto i_neighs = this->d_neighbor_p->getNeighbors(i);
        for (size_t j = 0; j < i_neighs.size(); j++) {

          auto j_id = i_neighs[j];

          // there are two contributions to force at node i
          // 1. From bond j-i due to bond-based forces
          // 2. From hydrostatic strains at node j and node i due to
          // hydrostatic forces

          // compute bond-based contribution
          auto xj = this->d_mesh_p->getNode(j_id);
          auto uj = this->d_u[j_id];
          auto rji = xj.dist(xi);
          auto Sji = this->d_material_p->getS(xj - xi, uj - ui);

          // get corrected volume of node j
          auto volj = this->d_mesh_p->getNodalVolume(j_id);
          if (util::compare::definitelyGreaterThan(rji, check_low))
            volj *= (check_up - rji) / h;

          // check if bond is in no-fail region
          bool break_bonds = true;
          if (!node_i_interior ||
              !this->d_interiorFlags_p->getInteriorFlag(j_id, xj))
            break_bonds = false;

          // get peridynamics force and energy density between bond i and j
          auto fs = this->d_fracture_p->getBondState(i, j);
          auto ef = this->d_material_p->getBondEF(rji, Sji, fs, break_bonds);
          this->d_fracture_p->setBondState(i, j, fs);

          // compute the contribution of bond force to force at i
          auto scalar_f = ef.second * volj / rji;

          force_i.d_x += scalar_f * (xj.d_x - xi.d_x);
          force_i.d_y += scalar_f * (xj.d_y - xi.d_y);
          force_i.d_z += scalar_f * (xj.d_z - xi.d_z);

          // compute state-based contribution
          if (this->d_material_p->isStateActive() &&
              this->d_material_p->doesBondContribToState(Sji, rji)) {

            // Compute state force density
            auto dgi = this->d_material_p->getStateForce(this->d_hS[i], rji);
            auto dgj = this->d_material_p->getStateForce(this->d_hS[j_id], rji);

            auto scalar_g = (dgi + dgj) / rji;

            force_i.d_x += scalar_g * (xj.d_x - xi.d_x);
            force_i.d_y += scalar_g * (xj.d_y - xi.d_y);
            force_i.d_z += scalar_g * (xj.d_z - xi.d_z);
          }
        } // loop over neighboring nodes

        // update force and energy
        this->d_f[i] += force_i;
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
        auto xi = this->d_mesh_p->getNode(i);
        auto ui = this->d_u[i];

        // upper and lower bound for volume correction
        auto h = d_mesh_p->getMeshSize();
        auto check_up = d_modelDeck_p->d_horizon + 0.5 * h;
        auto check_low = d_modelDeck_p->d_horizon - 0.5 * h;

        // inner loop over neighbors
        auto i_neighs = this->d_neighbor_p->getNeighbors(i);
        for (size_t j = 0; j < i_neighs.size(); j++) {

          auto j_id = i_neighs[j];

          // compute bond-based contribution
          auto uji = this->d_u[j_id] - ui;
          auto xji = this->d_mesh_p->getNode(j_id) - xi;
          auto rji = xji.length();
          auto Sji = this->d_material_p->getS(xji, uji);

          // get corrected volume of node j
          auto volj = this->d_mesh_p->getNodalVolume(j_id);
          if (util::compare::definitelyGreaterThan(rji, check_low))
            volj *= (check_up - rji) / h;

          // compute the contribution of bond to hydrostatic strain at i
          if (this->d_material_p->doesBondContribToState(Sji, rji))
            hydro_strain_i +=
                volj * this->d_material_p->getBondContribToHydroStrain(Sji, rji);
        } // loop over neighboring nodes

        // update
        this->d_hS[i] = hydro_strain_i;
      } // loop over nodes

  ); // end of parallel for loop

  f.get();
}

void model::FDModel::computePostProcFields() {

  // if work done is to be computed, get the external forces
  std::vector<util::Point3> f_ext;
  if (d_policy_p->populateData("Model_d_w")) {
    f_ext = std::vector<util::Point3>(d_mesh_p->getNumNodes(), util::Point3());
    d_fLoading_p->apply(d_time, &f_ext, d_mesh_p);
  }

  // local data for kinetic energy
  std::vector<float> vec_ke;
  if (this->d_policy_p->populateData("Model_d_e"))
    vec_ke = std::vector<float>(d_mesh_p->getNumNodes(), 0.);

  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_mesh_p->getNumNodes(),
      [this, &f_ext, &vec_ke](boost::uint64_t i) {
        // local variable
        double energy_i = 0.0;
        double hydro_energy_i = 0.0;
        double a = 0.; // for damage
        double b = 0.; // for damage
        double z = 0.; // for damage

        // reference coordinate and displacement at the node
        auto xi = this->d_mesh_p->getNode(i);
        auto ui = this->d_u[i];

        // get volume of node i
        auto voli = this->d_mesh_p->getNodalVolume(i);

        // get interior flag
        auto node_i_interior = this->d_interiorFlags_p->getInteriorFlag(i, xi);

        // upper and lower bound for volume correction
        auto h = d_mesh_p->getMeshSize();
        auto check_up = d_modelDeck_p->d_horizon + 0.5 * h;
        auto check_low = d_modelDeck_p->d_horizon - 0.5 * h;

        // inner loop over neighbors
        auto i_neighs = this->d_neighbor_p->getNeighbors(i);
        for (size_t j = 0; j < i_neighs.size(); j++) {

          auto j_id = i_neighs[j];

          // there are two contributions to force at node i
          // 1. From bond j-i due to bond-based forces
          // 2. From hydrostatic strains at node j and node i due to
          // hydrostatic forces

          // compute bond-based contribution
          auto xj = this->d_mesh_p->getNode(j_id);
          auto uj = this->d_u[j_id];
          auto rji = xj.dist(xi);
          auto Sji = this->d_material_p->getS(xj - xi, uj - ui);

          // get corrected volume of node j
          auto volj = this->d_mesh_p->getNodalVolume(j_id);
          if (util::compare::definitelyGreaterThan(rji, check_low))
            volj *= (check_up - rji) / h;

          // check if bond is in no-fail region
          bool break_bonds = true;
          if (!node_i_interior ||
              !this->d_interiorFlags_p->getInteriorFlag(j_id, xj))
            break_bonds = false;

          // get peridynamics force and energy density between bond i and j
          auto fs = this->d_fracture_p->getBondState(i, j);
          auto ef = this->d_material_p->getBondEF(rji, Sji, fs, break_bonds);

          // energy
          energy_i += ef.first * volj;

          // parameters for damage function \phi
          if (!fs)
            a += volj;
          b += volj;

          // parameters for damage function Z
          double sr = 0.;
          if (util::compare::definitelyGreaterThan(rji, 1.0E-12))
            sr = std::abs(Sji) / this->d_material_p->getSc(rji);
          if (util::compare::definitelyLessThan(z, sr))
            z = sr;
        } // loop over neighboring nodes

        // compute hydrostatic energy
        if (this->d_material_p->isStateActive())
          hydro_energy_i = this->d_material_p->getStateEnergy(this->d_hS[i]);

        if (this->d_policy_p->populateData("Model_d_e"))
          this->d_e[i] = (energy_i + hydro_energy_i) * voli;

        if (this->d_policy_p->populateData("Model_d_w"))
          this->d_w[i] = ui.dot(f_ext[i]);

        if (this->d_policy_p->populateData("Model_d_eFB") &&
            util::compare::definitelyGreaterThan(z, 1.0 - 1.0E-10))
          this->d_eFB[i] = energy_i * voli;

        if (this->d_policy_p->populateData("Model_d_eF") &&
            util::compare::definitelyGreaterThan(z, 1.0 - 1.0E-10))
          this->d_eF[i] = (energy_i + hydro_energy_i) * voli;

        if (this->d_policy_p->populateData("Model_d_phi"))
          this->d_phi[i] = 1. - a / b;

        if (this->d_policy_p->populateData("Model_d_Z"))
          this->d_Z[i] = z;

        // compute kinetic energy
        if (this->d_policy_p->populateData("Model_d_e"))
          vec_ke[i] = 0.5 * this->d_material_p->getDensity() *
                    this->d_v[i].dot(this->d_v[i]) * voli;
      } // loop over nodes

  ); // end of parallel for loop

  f.get();

  // add energies to get total energy
  if (this->d_policy_p->populateData("Model_d_e"))
    d_te = util::methods::add(d_e);
  if (this->d_policy_p->populateData("Model_d_w"))
    d_tw = util::methods::add(d_w);
  if (this->d_policy_p->populateData("Model_d_eF"))
    d_teF = util::methods::add(d_eF);
  if (this->d_policy_p->populateData("Model_d_eFB"))
    d_teFB = util::methods::add(d_eFB);

  if (this->d_policy_p->populateData("Model_d_e"))
    d_tk = util::methods::add(vec_ke);
}

void model::FDModel::output() {

  std::cout << "Output: time step = " << d_n << "\n";

  // write out % completion of simulation at 10% interval
  {
    float p = float(d_n) * 100. / d_modelDeck_p->d_Nt;
    int m = std::max(1, int(d_modelDeck_p->d_Nt / 10));
    if (d_n % m == 0 && int(p) > 0)
      std::cout << "Message: Simulation " << int(p) << "% complete.\n";
  }

  // filename
  std::string filename = d_outputDeck_p->d_path + "output_" +
                         std::to_string(d_n / d_outputDeck_p->d_dtOut);

  // open
  auto writer =
      rw::writer::VtkWriterInterface(filename, d_outputDeck_p->d_compressType);

  // write mesh
  if (d_mesh_p->getNumElements() != 0 && d_outputDeck_p->d_performFEOut)
    writer.appendMesh(d_mesh_p->getNodesP(), d_mesh_p->getElementType(),
        d_mesh_p->getElementConnectivitiesP(), &d_u);
  else
    writer.appendNodes(d_mesh_p->getNodesP(), &d_u);

  //
  // major simulation data
  //
  std::string tag = "Displacement";
  if (d_outputDeck_p->isTagInOutput(tag))
    writer.appendPointData(tag, &d_u);

  tag = "Velocity";
  if (d_outputDeck_p->isTagInOutput(tag))
    writer.appendPointData(tag, &d_v);

  tag = "Force";
  if (d_outputDeck_p->isTagInOutput(tag)) {

    std::vector<util::Point3> force(d_mesh_p->getNumNodes(), util::Point3());

    for (size_t i = 0; i < d_f.size(); i++)
      force[i] = d_f[i] * d_mesh_p->getNodalVolume(i);

    writer.appendPointData(tag, &force);
  }

  tag = "time";
  writer.addTimeStep(d_time);

  //
  // minor simulation data
  //
  if (!d_policy_p->enablePostProcessing()) {
    writer.close();
    return;
  }

  tag = "Force_Density";
  if (d_outputDeck_p->isTagInOutput(tag))
    writer.appendPointData(tag, &d_f);

  tag = "Strain_Energy";
  if (d_outputDeck_p->isTagInOutput(tag) &&
      d_policy_p->populateData("Model_d_e"))
    writer.appendPointData(tag, &d_e);

  tag = "Work_Done";
  if (d_outputDeck_p->isTagInOutput(tag) &&
      d_policy_p->populateData("Model_d_w"))
    writer.appendPointData(tag, &d_w);

  tag = "Fixity";
  if (d_outputDeck_p->isTagInOutput(tag))
    writer.appendPointData(tag, d_mesh_p->getFixityP());

  tag = "Node_Volume";
  if (d_outputDeck_p->isTagInOutput(tag))
    writer.appendPointData(tag, d_mesh_p->getNodalVolumeP());

  tag = "Damage_Phi";
  if (d_outputDeck_p->isTagInOutput(tag) &&
      d_policy_p->populateData("Model_d_phi"))
    writer.appendPointData(tag, &d_phi);

  tag = "Damage_Z";
  if (d_outputDeck_p->isTagInOutput(tag) &&
      d_policy_p->populateData("Model_d_Z"))
    writer.appendPointData(tag, &d_Z);

  tag = "Fracture_Perienergy_Bond";
  if (d_outputDeck_p->isTagInOutput(tag) &&
      d_policy_p->populateData("Model_d_eFB"))
    writer.appendPointData(tag, &d_eFB);

  tag = "Fracture_Perienergy_Total";
  if (d_outputDeck_p->isTagInOutput(tag) &&
      d_policy_p->populateData("Model_d_eF"))
    writer.appendPointData(tag, &d_eF);

  tag = "Total_Energy";
  double te = d_te - d_tw + d_tk;
  if (d_outputDeck_p->isTagInOutput(tag) &&
      d_policy_p->populateData("Model_d_e"))
    writer.appendFieldData(tag, te);

  tag = "Total_Fracture_Perienergy_Bond";
  if (d_outputDeck_p->isTagInOutput(tag) &&
      d_policy_p->populateData("Model_d_eFB"))
    writer.appendFieldData(tag, d_teFB);

  tag = "Total_Fracture_Perienergy_Total";
  if (d_outputDeck_p->isTagInOutput(tag) &&
      d_policy_p->populateData("Model_d_eF"))
    writer.appendFieldData(tag, d_teF);

  writer.close();
}