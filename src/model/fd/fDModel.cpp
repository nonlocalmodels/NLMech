////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "fDModel.h"

#include "data/DataManager.h"

// utils
#include "rw/reader.h"
#include "rw/writer.h"
#include "util/compare.h"
#include "util/fastMethods.h"
#include "util/matrix.h"
#include "util/point.h"
#include "util/utilFunction.h"
#include "util/utilGeom.h"

// include high level class declarations
#include "fe/massMatrix.h"
#include "fe/mesh.h"
#include "geometry/dampingGeom.h"
#include "geometry/fracture.h"
#include "geometry/interiorFlags.h"
#include "geometry/neighbor.h"
#include "inp/decks/absborbingCondDeck.h"
#include "inp/decks/modelDeck.h"
#include "inp/decks/outputDeck.h"
#include "inp/decks/restartDeck.h"
#include "inp/decks/materialDeck.h"
#include "inp/input.h"
#include "inp/policy.h"
#include "loading/fLoading.h"
#include "loading/initialCondition.h"
#include "loading/uLoading.h"
#include "material/materialUtil.h"
#include "material/pdMaterial.h"
#include "material/materials.h"


// standard lib
#include <fstream>

model::FDModel::FDModel(inp::Input *deck) 
    : d_massMatrix_p(nullptr),
      d_fracture_p(nullptr),
      d_interiorFlags_p(nullptr),
      d_input_p(deck),
      d_policy_p(nullptr),
      d_initialCondition_p(nullptr),
      d_material_p(nullptr),
      d_dampingGeom_p(nullptr),
      d_stop(false) {
  d_modelDeck_p = deck->getModelDeck();
  d_outputDeck_p = deck->getOutputDeck();
  d_policy_p = inp::Policy::getInstance(d_input_p->getPolicyDeck());
  d_absorbingCondDeck_p = deck->getAbsorbingCondDeck();

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

  if (d_outputDeck_p->d_outFormat == "vtu")
    rw::reader::readVtuFileRestart(d_restartDeck_p->d_file, d_dataManager_p->getDisplacementP(), d_dataManager_p->getVelocityP(),
                                   d_dataManager_p->getMeshP()->getNodesP());
  else if (d_outputDeck_p->d_outFormat == "msh")
    rw::reader::readMshFileRestart(d_restartDeck_p->d_file, d_dataManager_p->getDisplacementP(), d_dataManager_p->getVelocityP(),
                                   d_dataManager_p->getMeshP()->getNodesP());

  // integrate in time
  integrate();
}

void model::FDModel::initHObjects() {
  std::cout << "FDModel: Initializing high level objects.\n";
  // read mesh data
  std::cout << "FDModel: Creating mesh.\n";

  d_dataManager_p->setMeshP(new fe::Mesh(d_input_p->getMeshDeck()));
  d_dataManager_p->getMeshP()->clearElementData();

  std::cout << "number of nodes = " << d_dataManager_p->getMeshP()->getNumNodes()
            << " number of elements = " << d_dataManager_p->getMeshP()->getNumElements() << "\n";

  // create neighbor list
  std::cout << "FDModel: Creating neighbor list.\n";

  d_dataManager_p->setNeighborP(new geometry::Neighbor(d_modelDeck_p->d_horizon,
                                        d_input_p->getNeighborDeck(),
                                        d_dataManager_p->getMeshP()->getNodesP()));

  // create fracture data
  std::cout << "FDModel: Creating edge crack if any and modifying the "
               "fracture state of bonds.\n";
  d_fracture_p = new geometry::Fracture(d_input_p->getFractureDeck(),
                                        d_dataManager_p->getMeshP()->getNodesP(),
                                        d_dataManager_p->getNeighborP()->getNeighborsListP());

  // create interior flags
  std::cout << "FDModel: Creating interior flags for nodes.\n";
  d_interiorFlags_p = new geometry::InteriorFlags(
      d_input_p->getInteriorFlagsDeck(), d_dataManager_p->getMeshP()->getNodesP(),
      d_dataManager_p->getMeshP()->getBoundingBox());

  // initialize initial condition class
  std::cout << "FDModel: Initializing initial condition object.\n";
  d_initialCondition_p =
      new loading::InitialCondition(d_input_p->getInitialConditionDeck());

  // initialize loading class
  std::cout << "FDModel: Initializing displacement loading object.\n";
  d_dataManager_p->setDisplacementLoadingP(new loading::ULoading(d_input_p->getLoadingDeck(), d_dataManager_p->getMeshP()));
  std::cout << "FDModel: Initializing force loading object.\n";
  d_dataManager_p->setForceLoadingP(new loading::FLoading(d_input_p->getLoadingDeck(), d_dataManager_p->getMeshP()));

  // initialize material class
  std::cout << "FDModel: Initializing material object.\n";
  d_material_p = new material::pd::Material(d_input_p->getMaterialDeck(),
                                            d_modelDeck_p->d_dim,
                                            d_modelDeck_p->d_horizon);

  // initialize damping geometry class
  std::cout << "FDModel: Initializing damping object.\n";
  d_dampingGeom_p = new geometry::DampingGeom(d_absorbingCondDeck_p, d_dataManager_p->getMeshP());
}

void model::FDModel::init() {
  std::cout << "FDModel: Initializing basic datas.\n";

  d_n = 0;
  d_time = 0.;

  // get number of nodes, total number of dofs (fixed and free together)
  size_t nnodes = d_dataManager_p->getMeshP()->getNumNodes();

  // initialize major simulation data
  d_dataManager_p->setDisplacementP(new std::vector<util::Point3>(nnodes, util::Point3()));
  d_dataManager_p->setVelocityP(new std::vector<util::Point3>(nnodes, util::Point3()));
  d_dataManager_p->setForceP(new std::vector<util::Point3>(nnodes, util::Point3()));

  // Allocate the reaction force vector
  if (d_outputDeck_p->isTagInOutput("Reaction_Force") or
      d_outputDeck_p->isTagInOutput("Total_Reaction_Force")) {

    d_dataManager_p->setReactionForceP(new std::vector<util::Point3>(nnodes, util::Point3()));
    d_dataManager_p->setTotalReactionForceP(new std::vector<double>(nnodes, 0.));
  }

  // check material type and if needed update policy
  if (!d_material_p->isStateActive()) {
    d_policy_p->addToTags(0, "Model_d_hS");
    d_policy_p->addToTags(0, "Model_d_eF");
  }

  // check if we need to compute and store hydrostatic strain
  if (d_material_p->isStateActive()) {
    d_thetaX = std::vector<double>(d_dataManager_p->getMeshP()->getNumNodes(), 0.);

    if (d_material_p->name() == "PDState") {
      d_mX = std::vector<double>(d_dataManager_p->getMeshP()->getNumNodes(), 0.);

      material::computeStateMx(
          d_dataManager_p->getMeshP()->getNodes(), d_dataManager_p->getMeshP()->getNodalVolumes(),
          d_dataManager_p->getNeighborP()->getNeighborsList(), d_dataManager_p->getMeshP()->getMeshSize(),
          d_material_p, d_mX, true);
    }
  }

  // initialize minor simulation data
  if (d_policy_p->enablePostProcessing()) {
    std::string tag = "Strain_Energy";
    if (d_outputDeck_p->isTagInOutput(tag)) {
      // this data is asked in output file
      // but check if policy allows its population
      if (d_policy_p->populateData("Model_d_e"))
        d_dataManager_p->setStrainEnergyP(new std::vector<float>(nnodes, 0.0));
    } else {
      // this data is not asked in output thus we disable it
      d_policy_p->addToTags(0, "Model_d_e");
    }

    tag = "Work_Done";
    if (d_outputDeck_p->isTagInOutput(tag)) {
      if (d_policy_p->populateData("Model_d_w"))
        d_dataManager_p->setWorkDoneP(new std::vector<float>(nnodes, 0.0));
    } else
      d_policy_p->addToTags(0, "Model_d_w");

    tag = "Damage_Phi";
    if (d_outputDeck_p->isTagInOutput(tag)) {
      if (d_policy_p->populateData("Model_d_phi"))
        d_dataManager_p->setPhiP(new std::vector<float>(nnodes, 0.0));
    } else
      d_policy_p->addToTags(0, "Model_d_phi");

    tag = "Damage_Z";
    if (d_outputDeck_p->isTagInOutput(tag)) {
      if (d_policy_p->populateData("Model_d_Z"))
      d_dataManager_p->setDamageFunctionP(new std::vector<float>(nnodes, 0.0));
    } else
      d_policy_p->addToTags(0, "Model_d_Z");

    tag = "Fracture_Perienergy_Total";
    if (d_outputDeck_p->isTagInOutput(tag)) {
      if (d_policy_p->populateData("Model_d_eF"))
        d_eF = std::vector<float>(nnodes, 0.0);
    } else
      d_policy_p->addToTags(0, "Model_d_eF");

    tag = "Fracture_Perienergy_Bond";
    if (d_outputDeck_p->isTagInOutput(tag)) {
      if (d_policy_p->populateData("Model_d_eFB"))
      d_dataManager_p->setFractureEnergyP(new std::vector<float>(nnodes, 0.0));
    } else
      d_policy_p->addToTags(0, "Model_d_eFB");
  }

  if (d_outputDeck_p->d_outCriteria == "max_Z" or
      d_outputDeck_p->d_outCriteria == "max_Z_stop") {
    size_t num_params = 1;
    if (d_outputDeck_p->d_outCriteria == "max_Z_stop") num_params = 6;

    if (d_outputDeck_p->d_outCriteriaParams.size() < num_params) {
      std::cerr << "Error: Output criteria " << d_outputDeck_p->d_outCriteria
                << " requires " << num_params << " parameters. \n";
      exit(1);
    }

    // issue warning when postprocessing is turned off as we need damage data
    // to compute output criteria
    if (!d_policy_p->enablePostProcessing()) {
      std::cout << "Warning: Output criteria " << d_outputDeck_p->d_outCriteria
                << " requires Damage data Z but either "
                   "postprocessing is set to off. "
                   "Therefore setting output criteria to "
                   "null.\n";
      d_outputDeck_p->d_outCriteria.clear();
    } else {
      // check if damage data is allocated
      if ((*d_dataManager_p->getDamageFunctionP()).size() != d_dataManager_p->getMeshP()->getNumNodes()) {
        // allocate data
        d_dataManager_p->setDamageFunctionP(new std::vector<float>(nnodes, 0.0));

        // check if damage data is allowed in policy class (if not, need to
        // allow it by removing the tag related to damage function Z)
        if (!d_policy_p->populateData("Model_d_Z"))
          d_policy_p->removeTag("Model_d_Z");
      }
    }
  }  // handle output criteria exceptions
}

void model::FDModel::integrate() {
  // apply initial loading
  if (d_n == 0) d_initialCondition_p->apply(d_dataManager_p->getDisplacementP(), d_dataManager_p->getVelocityP(), d_dataManager_p->getMeshP());

  // apply loading
  d_dataManager_p->getDisplacementLoadingP()->apply(d_time, d_dataManager_p->getDisplacementP(), d_dataManager_p->getVelocityP(), d_dataManager_p->getMeshP());
  d_dataManager_p->getForceLoadingP()->apply(d_time, d_dataManager_p->getForceP(), d_dataManager_p->getMeshP());

  // internal forces
  computeForces();

  // perform output at the beginning
  if (d_n == 0) {
    if (d_policy_p->enablePostProcessing()) computePostProcFields();

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
      if (d_policy_p->enablePostProcessing()) computePostProcFields();

      output();

      // exit early if output criteria has changed the d_stop flag to true
      if (d_stop) return;

      // check if we need to modify the output frequency
      checkOutputCriteria();
    }

    // check for crack application
    if (d_fracture_p->addCrack(d_time, d_dataManager_p->getMeshP()->getNodesP(),
                               d_dataManager_p->getNeighborP()->getNeighborsListP())) {
      // check if we need to modify the output frequency
      checkOutputCriteria();
    }
  }  // loop over time steps
}

void model::FDModel::integrateCD() {
  // parallel for loop
  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_dataManager_p->getMeshP()->getNumNodes(), [this](boost::uint64_t i) {
        auto dim = this->d_dataManager_p->getMeshP()->getDimension();
        auto delta_t = this->d_modelDeck_p->d_dt;
        auto fact = delta_t * delta_t / this->d_material_p->getDensity();

        if (this->d_dataManager_p->getMeshP()->isNodeFree(i, 0)) {
          auto u_old = (*d_dataManager_p->getDisplacementP())[i].d_x;

          (*d_dataManager_p->getDisplacementP())[i].d_x +=
              fact * (*this->d_dataManager_p->getForceP())[i].d_x + delta_t * (*d_dataManager_p->getVelocityP())[i].d_x;

          (*d_dataManager_p->getVelocityP())[i].d_x = ((*d_dataManager_p->getDisplacementP())[i].d_x - u_old) / delta_t;
        }

        if (dim > 1)
          if (this->d_dataManager_p->getMeshP()->isNodeFree(i, 1)) {
            auto u_old = (*d_dataManager_p->getDisplacementP())[i].d_y;

           (*d_dataManager_p->getDisplacementP())[i].d_y +=
                fact * (*this->d_dataManager_p->getForceP())[i].d_y + delta_t * (*d_dataManager_p->getVelocityP())[i].d_y;

            (*d_dataManager_p->getVelocityP())[i].d_y = ((*d_dataManager_p->getDisplacementP())[i].d_y - u_old) / delta_t;
          }

        if (dim > 2)
          if (this->d_dataManager_p->getMeshP()->isNodeFree(i, 2)) {
            auto u_old = (*d_dataManager_p->getDisplacementP())[i].d_z;

            (*d_dataManager_p->getDisplacementP())[i].d_z +=
                fact * (*this->d_dataManager_p->getForceP())[i].d_z + delta_t * (*d_dataManager_p->getVelocityP())[i].d_z;

            (*d_dataManager_p->getVelocityP())[i].d_z = ((*d_dataManager_p->getDisplacementP())[i].d_z - u_old) / delta_t;
          }

        // reset force
        (*this->d_dataManager_p->getForceP())[i] = util::Point3();
      });  // end of parallel for loop

  f.get();

  // compute forces and energy due to new displacement field (this will be
  // used in next time step)
  d_n++;
  d_time += d_modelDeck_p->d_dt;

  // boundary condition
  d_dataManager_p->getDisplacementLoadingP()->apply(d_time, d_dataManager_p->getDisplacementP(), d_dataManager_p->getVelocityP(), d_dataManager_p->getMeshP());
  d_dataManager_p->getForceLoadingP()->apply(d_time, d_dataManager_p->getForceP(), d_dataManager_p->getMeshP());

  // internal forces
  computeForces();
}

void model::FDModel::integrateVerlet() {
  // step 1 and 2 : Compute v_mid and u_new
  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_dataManager_p->getMeshP()->getNumNodes(), [this](boost::uint64_t i) {
        auto dim = this->d_dataManager_p->getMeshP()->getDimension();
        auto delta_t = this->d_modelDeck_p->d_dt;
        auto fact = 0.5 * delta_t / this->d_material_p->getDensity();

        // modify dofs which are not marked fixed
        if (this->d_dataManager_p->getMeshP()->isNodeFree(i, 0)) {
          (*d_dataManager_p->getVelocityP())[i].d_x += fact * (*this->d_dataManager_p->getForceP())[i].d_x;
          (*d_dataManager_p->getDisplacementP())[i].d_x += delta_t * (*d_dataManager_p->getVelocityP())[i].d_x;
        }

        if (dim > 1)
          if (this->d_dataManager_p->getMeshP()->isNodeFree(i, 1)) {
            (*d_dataManager_p->getVelocityP())[i].d_y += fact * (*this->d_dataManager_p->getForceP())[i].d_y;
            (*d_dataManager_p->getDisplacementP())[i].d_y += delta_t * (*d_dataManager_p->getVelocityP())[i].d_y;
          }

        if (dim > 2)
          if (this->d_dataManager_p->getMeshP()->isNodeFree(i, 2)) {
            (*d_dataManager_p->getVelocityP())[i].d_z += fact * (*this->d_dataManager_p->getForceP())[i].d_z;
            (*d_dataManager_p->getDisplacementP())[i].d_z += delta_t * (*d_dataManager_p->getVelocityP())[i].d_z;
          }

        // reset force
       (*this->d_dataManager_p->getForceP())[i] = util::Point3();
      });  // end of parallel for loop

  f.get();

  // compute forces and energy due to new displacement field (this will be
  // used in next time step)
  d_n++;
  d_time += d_modelDeck_p->d_dt;

  // boundary condition
  d_dataManager_p->getDisplacementLoadingP()->apply(d_time, d_dataManager_p->getDisplacementP(), d_dataManager_p->getVelocityP(), d_dataManager_p->getMeshP());
  d_dataManager_p->getForceLoadingP()->apply(d_time, d_dataManager_p->getForceP(), d_dataManager_p->getMeshP());

  // internal forces
  computeForces();

  // Step 3: Compute v_new
  f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_dataManager_p->getMeshP()->getNumNodes(), [this](boost::uint64_t i) {
        auto dim = this->d_dataManager_p->getMeshP()->getDimension();
        auto fact =
            0.5 * this->d_modelDeck_p->d_dt / this->d_material_p->getDensity();

        // modify dofs which are not marked fixed
        if (this->d_dataManager_p->getMeshP()->isNodeFree(i, 0))
          (*d_dataManager_p->getVelocityP())[i].d_x += fact * (*this->d_dataManager_p->getForceP())[i].d_x;

        if (dim > 1)
          if (this->d_dataManager_p->getMeshP()->isNodeFree(i, 1))
            (*d_dataManager_p->getVelocityP())[i].d_y += fact * (*this->d_dataManager_p->getForceP())[i].d_y;

        if (dim > 2)
          if (this->d_dataManager_p->getMeshP()->isNodeFree(i, 2))
            (*d_dataManager_p->getVelocityP())[i].d_z += fact * (*this->d_dataManager_p->getForceP())[i].d_z;
      });  // end of parallel for loop

  f.get();
}

void model::FDModel::computeForces() {
  const auto &nodes = d_dataManager_p->getMeshP()->getNodes();
  const auto &volumes = d_dataManager_p->getMeshP()->getNodalVolumes();

  // compute helper quantities for state-based model
  if (d_material_p->isStateActive()) {
    if (d_material_p->name() == "RNPState")
      material::computeHydrostaticStrain(
          nodes,(*d_dataManager_p->getDisplacementP()), volumes, d_dataManager_p->getNeighborP()->getNeighborsList(),
          d_dataManager_p->getMeshP()->getMeshSize(), d_material_p, d_fracture_p, d_thetaX,
          d_modelDeck_p->d_dim, true);
    else if (d_material_p->name() == "PDState") {
      // need to update the fracture state of bonds
      material::updateBondFractureData(nodes, (*d_dataManager_p->getDisplacementP()),
                                       d_dataManager_p->getNeighborP()->getNeighborsList(),
                                       d_material_p, d_fracture_p, true);

      material::computeStateThetax(nodes, (*d_dataManager_p->getDisplacementP()), volumes,
                                   d_dataManager_p->getNeighborP()->getNeighborsList(),
                                   d_dataManager_p->getMeshP()->getMeshSize(), d_material_p,
                                   d_fracture_p, d_mX, d_thetaX, true);
    }
  }

  if (d_material_p->isStateActive()) {
    auto f = hpx::parallel::for_loop(
        hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
        d_dataManager_p->getMeshP()->getNumNodes(),
        [this](boost::uint64_t i) {
          (*this->d_dataManager_p->getForceP())[i] += this->computeForceState(i).second;
        }  // loop over nodes
    );     // end of parallel for loop
    f.get();
  } else {
    auto f = hpx::parallel::for_loop(
        hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
        d_dataManager_p->getMeshP()->getNumNodes(),
        [this](boost::uint64_t i) {
          (*this->d_dataManager_p->getForceP())[i] += this->computeForce(i).second;
        }  // loop over nodes
    );     // end of parallel for loop
    f.get();
  }
}

std::pair<double, util::Point3> model::FDModel::computeForce(const size_t &i) {
  // local variable to hold force
  auto force_i = util::Point3();
  double energy_i = 0.;

  if (d_outputDeck_p->isTagInOutput("Reaction_Force") or
      d_outputDeck_p->isTagInOutput("Total_Reaction_Force")) {
    (*d_dataManager_p->getReactionForceP())[i] = util::Point3();
    (*d_dataManager_p->getTotalReactionForceP())[i] = 0.;
  }

  // reference coordinate and displacement at the node
  auto xi = this->d_dataManager_p->getMeshP()->getNode(i);
  auto ui = (*d_dataManager_p->getDisplacementP())[i];

  // get interior flag
  auto node_i_interior = this->d_interiorFlags_p->getInteriorFlag(i, xi);

  // upper and lower bound for volume correction
  auto h = d_dataManager_p->getMeshP()->getMeshSize();
  auto check_up = d_modelDeck_p->d_horizon + 0.5 * h;
  auto check_low = d_modelDeck_p->d_horizon - 0.5 * h;

  // inner loop over neighbors
  const auto &i_neighs = this->d_dataManager_p->getNeighborP()->getNeighbors(i);
  for (size_t j = 0; j < i_neighs.size(); j++) {
    auto j_id = i_neighs[j];

    // there are two contributions to force at node i
    // 1. From bond j-i due to bond-based forces
    // 2. From hydrostatic strains at node j and node i due to
    // hydrostatic forces

    // compute bond-based contribution
    auto xj = this->d_dataManager_p->getMeshP()->getNode(j_id);
    auto uj = (*d_dataManager_p->getDisplacementP())[j_id];
    auto rji = xj.dist(xi);
    auto Sji = this->d_material_p->getS(xj - xi, uj - ui);

    // get corrected volume of node j
    auto volj = this->d_dataManager_p->getMeshP()->getNodalVolume(j_id);
    if (util::compare::definitelyGreaterThan(rji, check_low))
      volj *= (check_up - rji) / h;

    // check if bond is in no-fail region
    bool break_bonds = true;
    if (!node_i_interior || !this->d_interiorFlags_p->getInteriorFlag(j_id, xj))
      break_bonds = false;

    // get peridynamics force and energy density between bond i and j
    auto fs = this->d_fracture_p->getBondState(i, j);
    auto ef = this->d_material_p->getBondEF(rji, Sji, fs, break_bonds);

    // update fracture state of bond
    this->d_fracture_p->setBondState(i, j, fs);

    // compute the contribution of bond force to force at i
    auto scalar_f = ef.second * volj;

    force_i +=
        scalar_f * this->d_material_p->getBondForceDirection(xj - xi, uj - ui);

    energy_i += ef.first * volj;

    // Todo: Add reaction force computation
    if (is_reaction_force(i, j_id) and
        (d_outputDeck_p->isTagInOutput("Reaction_Force") or
         d_outputDeck_p->isTagInOutput("Total_Reaction_Force")))
          (*d_dataManager_p->getReactionForceP())[i] +=
          (this->d_dataManager_p->getMeshP()->getNodalVolume(i) * scalar_f *
           this->d_material_p->getBondForceDirection(xj - xi, uj - ui));

  }  // loop over neighboring nodes

  if (d_outputDeck_p->isTagInOutput("Total_Reaction_Force"))
    (*d_dataManager_p->getTotalReactionForceP())[i] = (*d_dataManager_p->getReactionForceP())[i].length();

  return std::make_pair(energy_i, force_i);
}

std::pair<double, util::Point3> model::FDModel::computeForceState(
    const size_t &i) {
  // local variable to hold force
  auto force_i = util::Point3();
  double energy_i = 0.;

  // reference coordinate and displacement at the node
  auto xi = this->d_dataManager_p->getMeshP()->getNode(i);
  auto ui = (*d_dataManager_p->getDisplacementP())[i];
  auto thetai = d_thetaX[i];
  auto mi = d_mX[i];

  // get interior flag
  auto node_i_interior = this->d_interiorFlags_p->getInteriorFlag(i, xi);

  // upper and lower bound for volume correction
  auto h = d_dataManager_p->getMeshP()->getMeshSize();
  auto check_up = d_modelDeck_p->d_horizon + 0.5 * h;
  auto check_low = d_modelDeck_p->d_horizon - 0.5 * h;

  // inner loop over neighbors
  const auto &i_neighs = this->d_dataManager_p->getNeighborP()->getNeighbors(i);
  for (size_t j = 0; j < i_neighs.size(); j++) {
    auto j_id = i_neighs[j];

    // there are two contributions to force at node i
    // 1. From bond j-i due to bond-based forces
    // 2. From hydrostatic strains at node j and node i due to
    // hydrostatic forces

    // compute bond-based contribution
    auto xj = this->d_dataManager_p->getMeshP()->getNode(j_id);
    auto uj = (*d_dataManager_p->getDisplacementP())[j_id];
    auto thetaj = d_thetaX[j];
    auto mj = d_mX[j];
    auto rji = xj.dist(xi);
    auto Sji = this->d_material_p->getS(xj - xi, uj - ui);

    // get corrected volume of node j
    auto volj = this->d_dataManager_p->getMeshP()->getNodalVolume(j_id);
    if (util::compare::definitelyGreaterThan(rji, check_low))
      volj *= (check_up - rji) / h;

    // check if bond is in no-fail region
    bool break_bonds = true;
    if (!node_i_interior || !this->d_interiorFlags_p->getInteriorFlag(j_id, xj))
      break_bonds = false;

    // get peridynamics force and energy density between bond i and j
    auto fs = this->d_fracture_p->getBondState(i, j);
    // auto ef_i = this->d_material_p->getBondEF(rji, Sji, fs, mi, thetai);
    // auto ef_j = this->d_material_p->getBondEF(rji, Sji, fs, mj, thetaj);
    auto ef_i = this->d_material_p->getBondEF(rji, Sji, fs, true);
    auto ef_j = this->d_material_p->getBondEF(rji, Sji, fs, true);

    // for state-based, we do not update bond-state as it is already updated
    // this->d_fracture_p->setBondState(i, j, fs);

    // compute the contribution of bond force to force at i
    auto scalar_f = (ef_i.second + ef_j.second) * volj;

    force_i +=
        scalar_f * this->d_material_p->getBondForceDirection(xj - xi, uj - ui);
  }  // loop over neighboring nodes

  return std::make_pair(energy_i, force_i);
}

bool model::FDModel::is_reaction_force(size_t i, size_t j) {
  auto xi = this->d_dataManager_p->getMeshP()->getNode(i);
  auto xj = this->d_dataManager_p->getMeshP()->getNode(j);
  auto delta = d_modelDeck_p->d_horizon;
  auto min_x = this->d_dataManager_p->getMeshP()->getBoundingBox().first[0];
  auto max_x = this->d_dataManager_p->getMeshP()->getBoundingBox().second[0];
  auto max_y = this->d_dataManager_p->getMeshP()->getBoundingBox().second[1];
  auto min_y = this->d_dataManager_p->getMeshP()->getBoundingBox().first[1];
  auto eps = 1e-6;

  // Make sure that the node is not in the boundary layer and in
  // the L layer
  if (xi.d_y > max_y - 1.5 * delta && xi.d_y < max_y - delta) {
    util::Point3 A = util::Point3(min_x, max_y - delta + eps, 0.);
    util::Point3 B = util::Point3(max_x, max_y - delta + eps, 0.);

    if (util::geometry::doLinesIntersect(A, B, xi, xj)) return true;
  }

  return false;
}

void model::FDModel::computeDampingForces() {
  if (!d_dampingGeom_p->isDampingActive()) return;

  double delta_t = d_modelDeck_p->d_dt;
  auto dim = d_modelDeck_p->d_dim;
  bool is_viscous_damping = d_dampingGeom_p->isViscousDamping();

  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_dataManager_p->getMeshP()->getNumNodes(),
      [this, delta_t, is_viscous_damping, dim](boost::uint64_t i) {
        // local variable to hold force
        util::Point3 force_i = util::Point3();
        auto v_i = (*d_dataManager_p->getVelocityP())[i];

        // compute damping force
        // -alpha* |F| sign(vel)
        double coeff = this->d_dampingGeom_p->getCoefficient(i);
        if (is_viscous_damping) {
          coeff *= delta_t;
          force_i =
              util::Point3(coeff * v_i.d_x, coeff * v_i.d_y, coeff * v_i.d_z);
        } else {
          coeff *= delta_t * delta_t * (*this->d_dataManager_p->getForceP())[i].length();
          auto sv = util::function::signVector(v_i);

          force_i.d_x = -coeff * sv.d_x;

          if (dim > 1) force_i.d_y = -coeff * sv.d_y;

          if (dim > 2) force_i.d_z = -coeff * sv.d_z;
        }

        // add to the force
       (*this->d_dataManager_p->getForceP())[i] += force_i;
      }  // loop over nodes

  );  // end of parallel for loop

  f.get();
}

void model::FDModel::computePostProcFields() {
  std::cout << "Postprocessing\n";

  // if work done is to be computed, get the external forces
  std::vector<util::Point3> f_ext;
  if (d_policy_p->populateData("Model_d_w")) {
    f_ext = std::vector<util::Point3>(d_dataManager_p->getMeshP()->getNumNodes(), util::Point3());
    d_dataManager_p->getForceLoadingP()->apply(d_time, &f_ext, d_dataManager_p->getMeshP());
  }

  // local data for kinetic energy
  std::vector<float> vec_ke;
  if (this->d_policy_p->populateData("Model_d_e"))
    vec_ke = std::vector<float>(d_dataManager_p->getMeshP()->getNumNodes(), 0.);

  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_dataManager_p->getMeshP()->getNumNodes(),
      [this, &f_ext, &vec_ke](boost::uint64_t i) {
        // local variable
        double energy_i = 0.0;
        double hydro_energy_i = 0.0;
        double a = 0.;  // for damage
        double b = 0.;  // for damage
        double z = 0.;  // for damage

        // reference coordinate and displacement at the node
        auto xi = this->d_dataManager_p->getMeshP()->getNode(i);
        auto ui = (*d_dataManager_p->getDisplacementP())[i];

        // get volume of node i
        auto voli = this->d_dataManager_p->getMeshP()->getNodalVolume(i);

        // get interior flag
        auto node_i_interior = this->d_interiorFlags_p->getInteriorFlag(i, xi);

        // upper and lower bound for volume correction
        auto h = d_dataManager_p->getMeshP()->getMeshSize();
        auto check_up = d_modelDeck_p->d_horizon + 0.5 * h;
        auto check_low = d_modelDeck_p->d_horizon - 0.5 * h;

        // inner loop over neighbors
        const auto &i_neighs = this->d_dataManager_p->getNeighborP()->getNeighbors(i);
        for (size_t j = 0; j < i_neighs.size(); j++) {
          auto j_id = i_neighs[j];

          // there are two contributions to force at node i
          // 1. From bond j-i due to bond-based forces
          // 2. From hydrostatic strains at node j and node i due to
          // hydrostatic forces

          // compute bond-based contribution
          auto xj = this->d_dataManager_p->getMeshP()->getNode(j_id);
          auto uj = (*d_dataManager_p->getDisplacementP())[j_id];
          auto rji = xj.dist(xi);
          auto Sji = this->d_material_p->getS(xj - xi, uj - ui);

          // get corrected volume of node j
          auto volj = this->d_dataManager_p->getMeshP()->getNodalVolume(j_id);
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
          if (!fs) a += volj;
          b += volj;

          // parameters for damage function Z
          double sr = 0.;
          if (util::compare::definitelyGreaterThan(rji, 1.0E-12))
            sr = std::abs(Sji) / this->d_material_p->getSc(rji);
          if (util::compare::definitelyLessThan(z, sr)) z = sr;
        }  // loop over neighboring nodes

        // compute hydrostatic energy
        //        if (this->d_material_p->isStateActive())
        //          hydro_energy_i =
        //          this->d_material_p->getStateEnergy(this->d_hS[i]);

        if (this->d_policy_p->populateData("Model_d_e"))
          (*d_dataManager_p->getStrainEnergyP())[i] = (energy_i + hydro_energy_i) * voli;

        if (this->d_policy_p->populateData("Model_d_w"))
          (*d_dataManager_p->getWorkDoneP())[i] = ui.dot(f_ext[i]);

        if (this->d_policy_p->populateData("Model_d_eFB") &&
            util::compare::definitelyGreaterThan(z, 1.0 - 1.0E-10))
          (*d_dataManager_p->getFractureEnergyP())[i] = energy_i * voli;

        if (this->d_policy_p->populateData("Model_d_eF") &&
            util::compare::definitelyGreaterThan(z, 1.0 - 1.0E-10))
          this->d_eF[i] = (energy_i + hydro_energy_i) * voli;

        if (this->d_policy_p->populateData("Model_d_phi"))
          (*d_dataManager_p->getPhiP())[i] = 1. - a / b;

        if (this->d_policy_p->populateData("Model_d_Z")) (*d_dataManager_p->getDamageFunctionP())[i] = z;

        // compute kinetic energy
        if (this->d_policy_p->populateData("Model_d_e"))
          vec_ke[i] = 0.5 * this->d_material_p->getDensity() *
                      (*d_dataManager_p->getVelocityP())[i].dot((*d_dataManager_p->getVelocityP())[i]) * voli;
      }  // loop over nodes

  );  // end of parallel for loop

  f.get();

  // add energies to get total energy
  if (this->d_policy_p->populateData("Model_d_e"))
    d_te = util::methods::add((*d_dataManager_p->getStrainEnergyP()));
  if (this->d_policy_p->populateData("Model_d_w"))
    d_tw = util::methods::add((*d_dataManager_p->getWorkDoneP()));
  if (this->d_policy_p->populateData("Model_d_eF"))
    d_teF = util::methods::add(d_eF);
  if (this->d_policy_p->populateData("Model_d_eFB"))
    d_teFB = util::methods::add((*d_dataManager_p->getFractureEnergyP()));

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
  // use smaller dt_out as the tag for files
  size_t dt_out = d_outputDeck_p->d_dtOutCriteria;
  std::string filename =
      d_outputDeck_p->d_path + "output_" + std::to_string(d_n / dt_out);

  // open
  auto writer = rw::writer::Writer(filename, d_outputDeck_p->d_outFormat,
                                   d_outputDeck_p->d_compressType);

  // write mesh
  //if (d_dataManager_p->getMeshP()->getNumElements() != 0 && d_outputDeck_p->d_performFEOut)
   // writer.appendMesh(d_dataManager_p->getMeshP()->getNodesP(), d_dataManager_p->getMeshP()->getElementType(),
   //                   d_dataManager_p->getMeshP()->getElementConnectivitiesP(), &d_u);
  //else
  //  writer.appendNodes(d_dataManager_p->getMeshP()->getNodesP(), &d_u);

  //
  // major simulation data
  //
  std::string tag = "Displacement";
  //if (d_outputDeck_p->isTagInOutput(tag)) writer.appendPointData(tag, &d_u);

  //tag = "Velocity";
  //if (d_outputDeck_p->isTagInOutput(tag)) writer.appendPointData(tag, &d_v);

/*
  tag = "Force";
  if (d_outputDeck_p->isTagInOutput(tag)) {
    std::vector<util::Point3> force(d_dataManager_p->getMeshP()->getNumNodes(), util::Point3());

    for (size_t i = 0; i < (*this->d_dataManager_p->getForceP()).size(); i++)
      force[i] = (*this->d_dataManager_p->getForceP())[i] * d_dataManager_p->getMeshP()->getNodalVolume(i);

    writer.appendPointData(tag, &force);
  }
  */

  //tag = "Reaction_Force";
  //if (d_outputDeck_p->isTagInOutput(tag)) {
  //  writer.appendPointData(tag, &d_reaction_force);
  // }
  //tag = "Total_Reaction_Force";
  //if (d_outputDeck_p->isTagInOutput(tag)) {
  //  double sum = std::accumulate((*d_dataManager_p->getTotalReactionForceP()).begin(),
   //                              (*d_dataManager_p->getTotalReactionForceP()).end(), 0);

    // Computation of the area
    //auto delta = d_modelDeck_p->d_horizon;
    //auto min_x = this->d_dataManager_p->getMeshP()->getBoundingBox().first[0];
    //auto max_x = this->d_dataManager_p->getMeshP()->getBoundingBox().second[0];
    //auto max_y = this->d_dataManager_p->getMeshP()->getBoundingBox().second[1];
    //auto min_y = this->d_dataManager_p->getMeshP()->getBoundingBox().first[1];

    //double area =
    //    (std::abs(max_x - min_x) - delta) * (std::abs(max_y - min_y) - delta);

    //writer.appendFieldData("Total_Reaction_Force", sum * area);
  //}

  // //
  // // debug
  // //
  // if (d_dampingGeom_p->isDampingActive()) {
  //   tag = "Damping_Coefficients";
  //   writer.appendPointData(tag, d_dampingGeom_p->getCoefficientDataP());
  // }

  //tag = "time";
  //writer.addTimeStep(d_time);

  //
  // minor simulation data
  //
  if (!d_policy_p->enablePostProcessing()) {
    writer.close();
    return;
  }

  //tag = "Force_Density";
  //if (d_outputDeck_p->isTagInOutput(tag)) writer.appendPointData(tag, &d_f);

  /*
  tag = "Strain_Energy";
  if (d_outputDeck_p->isTagInOutput(tag) &&
      d_policy_p->populateData("Model_d_e"))
    writer.appendPointData(tag, &d_e);
  *

  tag = "Work_Done";
  if (d_outputDeck_p->isTagInOutput(tag) &&
      d_policy_p->populateData("Model_d_w"))
    writer.appendPointData(tag, &d_w);

  //tag = "Fixity";
  //if (d_outputDeck_p->isTagInOutput(tag))
  //  writer.appendPointData(tag, d_dataManager_p->getMeshP()->getFixityP());

  //tag = "Node_Volume";
  //if (d_outputDeck_p->isTagInOutput(tag))
  //  writer.appendPointData(tag, d_dataManager_p->getMeshP()->getNodalVolumesP());

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

  /*
  tag = "Neighbors";
  if (d_outputDeck_p->isTagInOutput(tag)) {
    std::vector<size_t> amountNeighbors;
    size_t nodes = d_dataManager_p->getMeshP()->getNumNodes();
    for (size_t i = 0; i < nodes; i++)
      amountNeighbors.push_back(d_dataManager_p->getNeighborP()->getNeighbors(i).size());
    writer.appendPointData(tag, &amountNeighbors);
  }
*/
  writer.close();
}

void model::FDModel::checkOutputCriteria() {
  // if output criteria is empty then we do nothing
  // if we two output frequency specified by user is same then we do nothing
  if (d_outputDeck_p->d_outCriteria.empty() ||
      d_outputDeck_p->d_dtOutOld == d_outputDeck_p->d_dtOutCriteria)
    return;

  // perform checks every dt large intervals
  if (d_n % d_outputDeck_p->d_dtOutOld != 0) return;

  if (d_outputDeck_p->d_outCriteria == "max_Z" ||
      d_outputDeck_p->d_outCriteria == "max_Z_stop") {
    // since damage function will not reduce once attaining desired maximum
    // value, we do not check if the criteria was met in the past
    if (d_outputDeck_p->d_outCriteria == "max_Z" &&
        d_outputDeck_p->d_dtOut == d_outputDeck_p->d_dtOutCriteria)
      return;

    // change from large interval to small interval should be done only once
    // to do this, we check below flag which will be set to true in first
    // change to small from small interval
    static bool changed_to_small = false;
    bool changed_to_small_at_current = false;
    if (d_outputDeck_p->d_dtOut > d_outputDeck_p->d_dtOutCriteria &&
        !changed_to_small) {
      // get maximum from the damage data
      auto max = util::methods::max((*d_dataManager_p->getDamageFunctionP()));

      // check if it is desired range and change output frequency
      if (util::compare::definitelyGreaterThan(
              max, d_outputDeck_p->d_outCriteriaParams[0])) {
        d_outputDeck_p->d_dtOut = d_outputDeck_p->d_dtOutCriteria;

        std::cout << "Message: Changing output interval to smaller value.\n";

        // dump this value
        std::ofstream fdump(d_outputDeck_p->d_path + "dt_out_change.info");
        fdump << "Large_Dt_To_Small_Dt:\n";
        fdump << "  N: " << d_n << "\n";
        fdump << "  dN: " << d_n / d_outputDeck_p->d_dtOutCriteria << "\n";
        fdump.close();

        changed_to_small = true;
        changed_to_small_at_current = true;
      }
    }  // if current dt out is larger

    // we now check (only if max_Z_stop is flag is provided) if we need to
    // revert to large time interval
    // we do this if time interval is small at present
    // Also do not change to large interval back if just now we changed it to
    // small interval

    // change from small to large interval should be done only once
    // to do this, we check below flag which will be set to true in first
    // change to large from small interval
    static bool changed_back_to_large = false;
    if (!changed_to_small_at_current && !changed_back_to_large &&
        d_outputDeck_p->d_outCriteria == "max_Z_stop" &&
        d_outputDeck_p->d_dtOut < d_outputDeck_p->d_dtOutOld) {
      // check maximum of Z function in the rectangle
      auto rect = std::make_pair(util::Point3(), util::Point3());
      auto ps = d_outputDeck_p->d_outCriteriaParams;
      auto refZ = ps[1];
      if (ps.size() == 6) {
        rect.first = util::Point3(ps[2], ps[3], 0.);
        rect.second = util::Point3(ps[4], ps[5], 0.);
      } else if (ps.size() == 8) {
        rect.first = util::Point3(ps[2], ps[3], ps[4]);
        rect.second = util::Point3(ps[5], ps[6], ps[7]);
      }

      // we divide the total number of nodes in parts and check in parallel
      std::vector<std::vector<size_t>> rect_ids;
      size_t N = 1000;
      if (N > d_dataManager_p->getMeshP()->getNumNodes()) N = d_dataManager_p->getMeshP()->getNumNodes();
      rect_ids.resize(N);

      auto f = hpx::parallel::for_loop(
          hpx::parallel::execution::par(hpx::parallel::execution::task), 0, N,
          [this, N, rect, refZ, &rect_ids](boost::uint64_t I) {
            size_t ibegin = I * N;
            size_t iend = (I + 1) * N;
            if (iend > d_dataManager_p->getMeshP()->getNumNodes()) iend = d_dataManager_p->getMeshP()->getNumNodes();
            for (size_t i = ibegin; i < iend; i++) {
              if (util::geometry::isPointInsideRectangle(
                      d_dataManager_p->getMeshP()->getNode(i), rect.first, rect.second))
                if (util::compare::definitelyGreaterThan((*d_dataManager_p->getDamageFunctionP())[i], refZ))
                  rect_ids[I].emplace_back(i);
            }  // loop over chunk of I
          }    // loop over chunks
      );
      f.get();

      bool found_valid_node = false;
      for (auto i : rect_ids) {
        if (found_valid_node) break;
        if (!i.empty()) found_valid_node = true;
      }
      if (found_valid_node) {
        d_outputDeck_p->d_dtOut = d_outputDeck_p->d_dtOutOld;

        std::cout << "Message: Changing output interval to larger value.\n";

        // dump this value
        std::ofstream fdump(d_outputDeck_p->d_path + "dt_out_change.info",
                            std::ios::app);
        fdump << "Small_Dt_To_Large_Dt:\n";
        fdump << "  N: " << d_n << "\n";
        fdump << "  dN: " << d_n / d_outputDeck_p->d_dtOutCriteria << "\n";
        fdump.close();

        changed_back_to_large = true;

        d_stop = true;
      }
    }  // if max_Z_stop and current dt out is smaller
  }    // if either max_Z or max_Z_stop criteria
}
