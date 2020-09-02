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
#include "inp/decks/loadingDeck.h"
#include "inp/policy.h"
#include "loading/fLoading.h"
#include "loading/initialCondition.h"
#include "loading/uLoading.h"
#include "material/materials.h"
#include "model/util.h"

// standard lib
#include <fstream>

template <class T>
model::FDModel<T>::FDModel(inp::Input *deck)
    : d_input_p(deck),
      d_policy_p(nullptr),
      d_initialCondition_p(nullptr),
      d_material_p(nullptr),
      d_dampingGeom_p(nullptr),
      d_stop(false) {

  d_dataManager_p = new data::DataManager();

  d_dataManager_p->setModelDeckP(deck->getModelDeck());
  d_dataManager_p->setOutputDeckP(deck->getOutputDeck());

  d_policy_p = inp::Policy::getInstance(d_input_p->getPolicyDeck());
  d_absorbingCondDeck_p = deck->getAbsorbingCondDeck();

  if (d_dataManager_p->getModelDeckP()->d_isRestartActive)
    restart(deck);
  else
    run(deck);
}

template <class T>
model::FDModel<T>::~FDModel() {
  delete d_dataManager_p->getMeshP();
  delete d_dataManager_p->getDisplacementLoadingP();
  delete d_dataManager_p->getForceLoadingP();
  delete d_dataManager_p->getNeighborP();
  delete d_dataManager_p->getFractureP();
  delete d_dataManager_p->getInteriorFlagsP();
  delete d_dataManager_p->getDisplacementP();
  delete d_dataManager_p->getVelocityP();
  delete d_dataManager_p->getForceP();
  delete d_dataManager_p->getModelDeckP();
  delete d_dataManager_p->getOutputDeckP();

  delete d_material_p;
  delete d_initialCondition_p;
  delete d_dampingGeom_p;

  delete d_dataManager_p;
}

template <class T>
void model::FDModel<T>::run(inp::Input *deck) {
  // first initialize all the high level data
  initHObjects();

  // now initialize remaining data
  init();

  // integrate in time
  integrate();
}

template <class T>
void model::FDModel<T>::restart(inp::Input *deck) {
  d_restartDeck_p = deck->getRestartDeck();

  // first initialize all the high level data
  initHObjects();

  // now initialize remaining data
  init();

  // set time step to step specified in restart deck
  d_n = d_restartDeck_p->d_step;
  d_time = double(d_n) * d_dataManager_p->getModelDeckP()->d_dt;

  // read displacement and velocity from restart file

  if (d_dataManager_p->getOutputDeckP()->d_outFormat == "vtu")
    rw::reader::readVtuFileRestart(d_restartDeck_p->d_file,
                                   d_dataManager_p->getDisplacementP(),
                                   d_dataManager_p->getVelocityP(),
                                   d_dataManager_p->getMeshP()->getNodesP());
  else if (d_dataManager_p->getOutputDeckP()->d_outFormat == "msh")
    rw::reader::readMshFileRestart(d_restartDeck_p->d_file,
                                   d_dataManager_p->getDisplacementP(),
                                   d_dataManager_p->getVelocityP(),
                                   d_dataManager_p->getMeshP()->getNodesP());

  // integrate in time
  integrate();
}

template <class T>
void model::FDModel<T>::initHObjects() {
  std::cout << "FDModel: Initializing high level objects.\n";
  // read mesh data
  std::cout << "FDModel: Creating mesh.\n";

  d_dataManager_p->setMeshP(new fe::Mesh(d_input_p->getMeshDeck()));
  d_dataManager_p->getMeshP()->clearElementData();

  std::cout << "number of nodes = "
            << d_dataManager_p->getMeshP()->getNumNodes()
            << " number of elements = "
            << d_dataManager_p->getMeshP()->getNumElements() << "\n";

  // create neighbor list
  std::cout << "FDModel: Creating neighbor list.\n";

  d_dataManager_p->setNeighborP(new geometry::Neighbor(
      d_dataManager_p->getModelDeckP()->d_horizon, d_input_p->getNeighborDeck(),
      d_dataManager_p->getMeshP()->getNodesP()));

  // create fracture data
  std::cout << "FDModel: Creating edge crack if any and modifying the "
               "fracture state of bonds.\n";
  d_dataManager_p->setFractureP(new geometry::Fracture(
      d_input_p->getFractureDeck(), d_dataManager_p->getMeshP()->getNodesP(),
      d_dataManager_p->getNeighborP()->getNeighborsListP()));

  // create interior flags
  std::cout << "FDModel: Creating interior flags for nodes.\n";
  d_dataManager_p->setInteriorFlagsP(new geometry::InteriorFlags(
      d_input_p->getInteriorFlagsDeck(),
      d_dataManager_p->getMeshP()->getNodesP(),
      d_dataManager_p->getMeshP()->getBoundingBox()));

  // initialize initial condition class
  std::cout << "FDModel: Initializing initial condition object.\n";
  d_initialCondition_p =
      new loading::InitialCondition(d_input_p->getInitialConditionDeck());

  // initialize loading class
  std::cout << "FDModel: Initializing displacement loading object.\n";
  d_dataManager_p->setDisplacementLoadingP(new loading::ULoading(
      d_input_p->getLoadingDeck(), d_dataManager_p->getMeshP()));
  std::cout << "FDModel: Initializing force loading object.\n";
  d_dataManager_p->setForceLoadingP(new loading::FLoading(
      d_input_p->getLoadingDeck(), d_dataManager_p->getMeshP()));

  // initialize material class
  std::cout << "FDModel: Initializing material object.\n";
  d_material_p = new T(d_input_p->getMaterialDeck(), d_dataManager_p);

  // initialize damping geometry class
  std::cout << "FDModel: Initializing damping object.\n";
  d_dampingGeom_p = new geometry::DampingGeom(d_absorbingCondDeck_p,
                                              d_dataManager_p->getMeshP());
}

template <class T>
void model::FDModel<T>::init() {
  std::cout << "FDModel: Initializing basic datas.\n";

  d_n = 0;
  d_time = 0.;

  // get number of nodes, total number of dofs (fixed and free together)
  size_t nnodes = d_dataManager_p->getMeshP()->getNumNodes();

  // initialize major simulation data
  d_dataManager_p->setDisplacementP(
      new std::vector<util::Point3>(nnodes, util::Point3()));
  d_dataManager_p->setVelocityP(
      new std::vector<util::Point3>(nnodes, util::Point3()));
  d_dataManager_p->setForceP(
      new std::vector<util::Point3>(nnodes, util::Point3()));

  // Allocate the reaction force vector
  if (d_dataManager_p->getOutputDeckP()->isTagInOutput("Reaction_Force") or
      d_dataManager_p->getOutputDeckP()->isTagInOutput("Total_Reaction_Force")) {
    d_dataManager_p->setReactionForceP(
        new std::vector<util::Point3>(nnodes, util::Point3()));
    d_dataManager_p->setTotalReactionForceP(
        new std::vector<double>(nnodes, 0.));
  }

  // initialize minor simulation data
  if (this->d_policy_p->populateData("Model_d_e"))
    d_dataManager_p->setKineticEnergyP(
        new std::vector<float>(d_dataManager_p->getMeshP()->getNumNodes(), 0.));

  if (d_policy_p->enablePostProcessing()) {
    std::string tag = "Strain_Energy";
    if (d_dataManager_p->getOutputDeckP()->isTagInOutput(tag)) {
      // this data is asked in output file
      // but check if policy allows its population
      if (d_policy_p->populateData("Model_d_e"))
        d_dataManager_p->setStrainEnergyP(new std::vector<float>(nnodes, 0.0));
    } else {
      // this data is not asked in output thus we disable it
      d_policy_p->addToTags(0, "Model_d_e");
    }

    tag = "Work_Done";
    if (d_dataManager_p->getOutputDeckP()->isTagInOutput(tag)) {
      if (d_policy_p->populateData("Model_d_w"))
        d_dataManager_p->setWorkDoneP(new std::vector<float>(nnodes, 0.0));
    } else
      d_policy_p->addToTags(0, "Model_d_w");

    tag = "Damage_Phi";
    if (d_dataManager_p->getOutputDeckP()->isTagInOutput(tag)) {
      if (d_policy_p->populateData("Model_d_phi"))
        d_dataManager_p->setPhiP(new std::vector<float>(nnodes, 0.0));
    } else
      d_policy_p->addToTags(0, "Model_d_phi");

    tag = "Damage_Z";
    if (d_dataManager_p->getOutputDeckP()->isTagInOutput(tag)) {
      if (d_policy_p->populateData("Model_d_Z"))
        d_dataManager_p->setDamageFunctionP(
            new std::vector<float>(nnodes, 0.0));
    } else
      d_policy_p->addToTags(0, "Model_d_Z");

    tag = "Fracture_Perienergy_Total";
    if (d_dataManager_p->getOutputDeckP()->isTagInOutput(tag)) {
      if (d_policy_p->populateData("Model_d_eF"))

        d_dataManager_p->setFractureEnergyP(
            new std::vector<float>(nnodes, 0.0));
    } else
      d_policy_p->addToTags(0, "Model_d_eF");

    tag = "Fracture_Perienergy_Bond";
    if (d_dataManager_p->getOutputDeckP()->isTagInOutput(tag)) {
      if (d_policy_p->populateData("Model_d_eFB"))
        d_dataManager_p->setBBFractureEnergyP(
            new std::vector<float>(nnodes, 0.0));
    } else
      d_policy_p->addToTags(0, "Model_d_eFB");
  }

  if (d_dataManager_p->getOutputDeckP()->d_outCriteria == "max_Z" or
      d_dataManager_p->getOutputDeckP()->d_outCriteria == "max_Z_stop") {
    size_t num_params = 1;
    if (d_dataManager_p->getOutputDeckP()->d_outCriteria == "max_Z_stop") num_params = 6;

    if (d_dataManager_p->getOutputDeckP()->d_outCriteriaParams.size() < num_params) {
      std::cerr << "Error: Output criteria " << d_dataManager_p->getOutputDeckP()->d_outCriteria
                << " requires " << num_params << " parameters. \n";
      exit(1);
    }

    // issue warning when postprocessing is turned off as we need damage data
    // to compute output criteria
    if (!d_policy_p->enablePostProcessing()) {
      std::cout << "Warning: Output criteria " << d_dataManager_p->getOutputDeckP()->d_outCriteria
                << " requires Damage data Z but either "
                   "postprocessing is set to off. "
                   "Therefore setting output criteria to "
                   "null.\n";
      d_dataManager_p->getOutputDeckP()->d_outCriteria.clear();
    } else {
      // check if damage data is allocated
      if ((*d_dataManager_p->getDamageFunctionP()).size() !=
          d_dataManager_p->getMeshP()->getNumNodes()) {
        // allocate data
        d_dataManager_p->setDamageFunctionP(
            new std::vector<float>(nnodes, 0.0));

        // check if damage data is allowed in policy class (if not, need to
        // allow it by removing the tag related to damage function Z)
        if (!d_policy_p->populateData("Model_d_Z"))
          d_policy_p->removeTag("Model_d_Z");
      }
    }
  }  // handle output criteria exceptions
}

template <class T>
void model::FDModel<T>::integrate() {
  // apply initial loading
  if (d_n == 0)
    d_initialCondition_p->apply(d_dataManager_p->getDisplacementP(),
                                d_dataManager_p->getVelocityP(),
                                d_dataManager_p->getMeshP());

  // apply loading
  d_dataManager_p->getDisplacementLoadingP()->apply(
      d_time, d_dataManager_p->getDisplacementP(),
      d_dataManager_p->getVelocityP(), d_dataManager_p->getMeshP());
  d_dataManager_p->getForceLoadingP()->apply(
      d_time, d_dataManager_p->getForceP(), d_dataManager_p->getMeshP());

  // internal forces
  computeForces();

  // perform output at the beginning
  if (d_n == 0) {
    if (d_policy_p->enablePostProcessing()) computePostProcFields();

    model::Output(d_input_p, d_dataManager_p, d_n, d_time);
  }

  // start time integration
  size_t i = d_n;
  for (i; i < d_dataManager_p->getModelDeckP()->d_Nt; i++) {
    if (d_dataManager_p->getModelDeckP()->d_timeDiscretization == "central_difference")
      integrateCD();
    else if (d_dataManager_p->getModelDeckP()->d_timeDiscretization == "velocity_verlet")
      integrateVerlet();

    // handle general output
    if ((d_n % d_dataManager_p->getOutputDeckP()->d_dtOut == 0) &&
        (d_n >= d_dataManager_p->getOutputDeckP()->d_dtOut)) {
      if (d_policy_p->enablePostProcessing()) computePostProcFields();

      model::Output(d_input_p, d_dataManager_p, d_n, d_time);

      // exit early if output criteria has changed the d_stop flag to true
      if (d_stop) return;

      // check if we need to modify the output frequency
      checkOutputCriteria();
    }

    // check for crack application
    if (d_dataManager_p->getFractureP()->addCrack(
            d_time, d_dataManager_p->getMeshP()->getNodesP(),
            d_dataManager_p->getNeighborP()->getNeighborsListP())) {
      // check if we need to modify the output frequency
      checkOutputCriteria();
    }
  }  // loop over time steps
}

template <class T>
void model::FDModel<T>::integrateCD() {
  // parallel for loop
  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_dataManager_p->getMeshP()->getNumNodes(), [this](boost::uint64_t i) {
        auto dim = this->d_dataManager_p->getMeshP()->getDimension();
        auto delta_t = this->d_dataManager_p->getModelDeckP()->d_dt;
        auto fact = delta_t * delta_t / this->d_material_p->getDensity();

        if (this->d_dataManager_p->getMeshP()->isNodeFree(i, 0)) {
          auto u_old = (*d_dataManager_p->getDisplacementP())[i].d_x;

          (*d_dataManager_p->getDisplacementP())[i].d_x +=
              fact * (*this->d_dataManager_p->getForceP())[i].d_x +
              delta_t * (*d_dataManager_p->getVelocityP())[i].d_x;

          (*d_dataManager_p->getVelocityP())[i].d_x =
              ((*d_dataManager_p->getDisplacementP())[i].d_x - u_old) / delta_t;
        }

        if (dim > 1)
          if (this->d_dataManager_p->getMeshP()->isNodeFree(i, 1)) {
            auto u_old = (*d_dataManager_p->getDisplacementP())[i].d_y;

            (*d_dataManager_p->getDisplacementP())[i].d_y +=
                fact * (*this->d_dataManager_p->getForceP())[i].d_y +
                delta_t * (*d_dataManager_p->getVelocityP())[i].d_y;

            (*d_dataManager_p->getVelocityP())[i].d_y =
                ((*d_dataManager_p->getDisplacementP())[i].d_y - u_old) /
                delta_t;
          }

        if (dim > 2)
          if (this->d_dataManager_p->getMeshP()->isNodeFree(i, 2)) {
            auto u_old = (*d_dataManager_p->getDisplacementP())[i].d_z;

            (*d_dataManager_p->getDisplacementP())[i].d_z +=
                fact * (*this->d_dataManager_p->getForceP())[i].d_z +
                delta_t * (*d_dataManager_p->getVelocityP())[i].d_z;

            (*d_dataManager_p->getVelocityP())[i].d_z =
                ((*d_dataManager_p->getDisplacementP())[i].d_z - u_old) /
                delta_t;
          }

        // reset force
        (*this->d_dataManager_p->getForceP())[i] = util::Point3();
      });  // end of parallel for loop

  f.get();

  // compute forces and energy due to new displacement field (this will be
  // used in next time step)
  d_n++;
  d_time += d_dataManager_p->getModelDeckP()->d_dt;

  // boundary condition
  d_dataManager_p->getDisplacementLoadingP()->apply(
      d_time, d_dataManager_p->getDisplacementP(),
      d_dataManager_p->getVelocityP(), d_dataManager_p->getMeshP());
  d_dataManager_p->getForceLoadingP()->apply(
      d_time, d_dataManager_p->getForceP(), d_dataManager_p->getMeshP());

  // internal forces
  computeForces();
}

template <class T>
void model::FDModel<T>::integrateVerlet() {

  // step 1 and 2 : Compute v_mid and u_new
  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_dataManager_p->getMeshP()->getNumNodes(), [this](boost::uint64_t i) {
        auto dim = this->d_dataManager_p->getMeshP()->getDimension();
        auto delta_t = this->d_dataManager_p->getModelDeckP()->d_dt;
        auto fact = 0.5 * delta_t / this->d_material_p->getDensity();

        // modify dofs which are not marked fixed
        if (this->d_dataManager_p->getMeshP()->isNodeFree(i, 0)) {
          (*d_dataManager_p->getVelocityP())[i].d_x +=
              fact * (*this->d_dataManager_p->getForceP())[i].d_x;
          (*d_dataManager_p->getDisplacementP())[i].d_x +=
              delta_t * (*d_dataManager_p->getVelocityP())[i].d_x;
        }

        if (dim > 1)
          if (this->d_dataManager_p->getMeshP()->isNodeFree(i, 1)) {
            (*d_dataManager_p->getVelocityP())[i].d_y +=
                fact * (*this->d_dataManager_p->getForceP())[i].d_y;
            (*d_dataManager_p->getDisplacementP())[i].d_y +=
                delta_t * (*d_dataManager_p->getVelocityP())[i].d_y;
          }

        if (dim > 2)
          if (this->d_dataManager_p->getMeshP()->isNodeFree(i, 2)) {
            (*d_dataManager_p->getVelocityP())[i].d_z +=
                fact * (*this->d_dataManager_p->getForceP())[i].d_z;
            (*d_dataManager_p->getDisplacementP())[i].d_z +=
                delta_t * (*d_dataManager_p->getVelocityP())[i].d_z;
          }

        // reset force
        (*this->d_dataManager_p->getForceP())[i] = util::Point3();
      });  // end of parallel for loop

  f.get();

  // compute forces and energy due to new displacement field (this will be
  // used in next time step)
  d_n++;
  d_time += d_dataManager_p->getModelDeckP()->d_dt;

  // boundary condition
  d_dataManager_p->getDisplacementLoadingP()->apply(
      d_time, d_dataManager_p->getDisplacementP(),
      d_dataManager_p->getVelocityP(), d_dataManager_p->getMeshP());
  d_dataManager_p->getForceLoadingP()->apply(
      d_time, d_dataManager_p->getForceP(), d_dataManager_p->getMeshP());

  // internal forces
  computeForces();

  // Step 3: Compute v_new
  f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_dataManager_p->getMeshP()->getNumNodes(), [this](boost::uint64_t i) {
        auto dim = this->d_dataManager_p->getMeshP()->getDimension();
        auto fact =
            0.5 * this->d_dataManager_p->getModelDeckP()->d_dt / this->d_material_p->getDensity();

        // modify dofs which are not marked fixed
        if (this->d_dataManager_p->getMeshP()->isNodeFree(i, 0))
          (*d_dataManager_p->getVelocityP())[i].d_x +=
              fact * (*this->d_dataManager_p->getForceP())[i].d_x;

        if (dim > 1)
          if (this->d_dataManager_p->getMeshP()->isNodeFree(i, 1))
            (*d_dataManager_p->getVelocityP())[i].d_y +=
                fact * (*this->d_dataManager_p->getForceP())[i].d_y;

        if (dim > 2)
          if (this->d_dataManager_p->getMeshP()->isNodeFree(i, 2))
            (*d_dataManager_p->getVelocityP())[i].d_z +=
                fact * (*this->d_dataManager_p->getForceP())[i].d_z;
      });  // end of parallel for loop

  f.get();
}

template <class T>
void model::FDModel<T>::computeForces() {
  const auto &nodes = d_dataManager_p->getMeshP()->getNodes();
  const auto &volumes = d_dataManager_p->getMeshP()->getNodalVolumes();

  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_dataManager_p->getMeshP()->getNumNodes(),
      [this](boost::uint64_t i) {
        (*this->d_dataManager_p->getForceP())[i] +=
            this->computeForce(i).second;
      }  // loop over nodes
  );     // end of parallel for loop
  f.get();
}

template <class T>
std::pair<double, util::Point3> model::FDModel<T>::computeForce(const size_t &i) {
  // local variable to hold force
  auto force_i = util::Point3();
  double energy_i = 0.;

  if (d_dataManager_p->getOutputDeckP()->isTagInOutput("Reaction_Force") or
      d_dataManager_p->getOutputDeckP()->isTagInOutput("Total_Reaction_Force")) {
    (*d_dataManager_p->getReactionForceP())[i] = util::Point3();
    (*d_dataManager_p->getTotalReactionForceP())[i] = 0.;
  }

  // inner loop over neighbors
  const auto &i_neighs = this->d_dataManager_p->getNeighborP()->getNeighbors(i);
  for (size_t j = 0; j < i_neighs.size(); j++) {
    auto j_id = i_neighs[j];

    auto fe_pair = d_material_p->getBondEF(i, j);
    force_i += fe_pair.first;
    energy_i += fe_pair.second;

    // Todo: Add reaction force computation
    if (is_reaction_force(i, j_id) and
        (d_dataManager_p->getOutputDeckP()->isTagInOutput("Reaction_Force") or
         d_dataManager_p->getOutputDeckP()->isTagInOutput("Total_Reaction_Force")))
      (*d_dataManager_p->getReactionForceP())[i] +=
          (this->d_dataManager_p->getMeshP()->getNodalVolume(i) *
            fe_pair.first);

  }  // loop over neighboring nodes

  if (d_dataManager_p->getOutputDeckP()->isTagInOutput("Total_Reaction_Force"))
    (*d_dataManager_p->getTotalReactionForceP())[i] =
        (*d_dataManager_p->getReactionForceP())[i].length();

  return std::make_pair(energy_i, force_i);
}

template <class T>
bool model::FDModel<T>::is_reaction_force(size_t i, size_t j) {
  auto xi = this->d_dataManager_p->getMeshP()->getNode(i);
  auto xj = this->d_dataManager_p->getMeshP()->getNode(j);
  auto delta = d_dataManager_p->getModelDeckP()->d_horizon;
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

template <class T>
void model::FDModel<T>::computeDampingForces() {
  if (!d_dampingGeom_p->isDampingActive()) return;

  double delta_t = d_dataManager_p->getModelDeckP()->d_dt;
  auto dim = d_dataManager_p->getModelDeckP()->d_dim;
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
          coeff *= delta_t * delta_t *
                   (*this->d_dataManager_p->getForceP())[i].length();
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

template <class T>
void model::FDModel<T>::computePostProcFields() {
  std::cout << "Postprocessing\n";

  // if work done is to be computed, get the external forces
  std::vector<util::Point3> f_ext;
  if (d_policy_p->populateData("Model_d_w")) {
    f_ext = std::vector<util::Point3>(
        d_dataManager_p->getMeshP()->getNumNodes(), util::Point3());
    d_dataManager_p->getForceLoadingP()->apply(d_time, &f_ext,
                                               d_dataManager_p->getMeshP());
  }

  // local data for kinetic energy
  std::vector<float> vec_ke;
  if (this->d_policy_p->populateData("Model_d_e"))
    vec_ke = (*d_dataManager_p->getKineticEnergyP());

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

        // upper and lower bound for volume correction
        auto h = d_dataManager_p->getMeshP()->getMeshSize();
        auto horizon = d_material_p->getHorizon();
        auto check_up = horizon + 0.5 * h;
        auto check_low = horizon - 0.5 * h;

        // get volume of node i
        auto voli = this->d_dataManager_p->getMeshP()->getNodalVolume(i);

        // inner loop over neighbors
        const auto &i_neighs =
            this->d_dataManager_p->getNeighborP()->getNeighbors(i);
        for (size_t j = 0; j < i_neighs.size(); j++) {
          auto j_id = i_neighs[j];

          auto fe_pair = this->d_material_p->getBondEF(i, j);
          auto fs = this->d_dataManager_p->getFractureP()->getBondState(i, j);

          // energy
          energy_i += fe_pair.second;

          auto xj = d_dataManager_p->getMeshP()->getNode(j_id);
          auto uj = (*d_dataManager_p->getDisplacementP())[j_id];
          auto rji = xj.dist(xi);
          auto Sji = this->d_material_p->getS(xj - xi, uj - ui);

          // get corrected volume of node j
          auto volj = d_dataManager_p->getMeshP()->getNodalVolume(j_id);
          if (util::compare::definitelyGreaterThan(rji, check_low))
            volj *= (check_up - rji) / h;

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
          (*d_dataManager_p->getStrainEnergyP())[i] =
              (energy_i + hydro_energy_i) * voli;

        if (this->d_policy_p->populateData("Model_d_w"))
          (*d_dataManager_p->getWorkDoneP())[i] = ui.dot(f_ext[i]);

        if (this->d_policy_p->populateData("Model_d_eFB") &&
            util::compare::definitelyGreaterThan(z, 1.0 - 1.0E-10))
          (*d_dataManager_p->getBBFractureEnergyP())[i] = energy_i * voli;

        if (this->d_policy_p->populateData("Model_d_eF") &&
            util::compare::definitelyGreaterThan(z, 1.0 - 1.0E-10))
          (*d_dataManager_p->getFractureEnergyP())[i] =
              (energy_i + hydro_energy_i) * voli;

        if (this->d_policy_p->populateData("Model_d_phi"))
          (*d_dataManager_p->getPhiP())[i] = 1. - a / b;

        if (this->d_policy_p->populateData("Model_d_Z"))
          (*d_dataManager_p->getDamageFunctionP())[i] = z;

        // compute kinetic energy
        if (this->d_policy_p->populateData("Model_d_e"))
          (*d_dataManager_p->getKineticEnergyP())[i] =
              0.5 * this->d_material_p->getDensity() *
              (*d_dataManager_p->getVelocityP())[i].dot(
                  (*d_dataManager_p->getVelocityP())[i]) *
              voli;
      }  // loop over nodes

  );  // end of parallel for loop

  f.get();

  // add energies to get total energy
  if (this->d_policy_p->populateData("Model_d_e"))
    d_te = util::methods::add((*d_dataManager_p->getStrainEnergyP()));
  if (this->d_policy_p->populateData("Model_d_w"))
    d_tw = util::methods::add((*d_dataManager_p->getWorkDoneP()));
  if (this->d_policy_p->populateData("Model_d_eF"))
    d_teF = util::methods::add((*d_dataManager_p->getFractureEnergyP()));
  if (this->d_policy_p->populateData("Model_d_eFB"))
    d_teFB = util::methods::add((*d_dataManager_p->getBBFractureEnergyP()));

  if (this->d_policy_p->populateData("Model_d_e"))
    d_tk = util::methods::add((*d_dataManager_p->getKineticEnergyP()));
}

template <class T>
void model::FDModel<T>::checkOutputCriteria() {
  // if output criteria is empty then we do nothing
  // if we two output frequency specified by user is same then we do nothing
  if (d_dataManager_p->getOutputDeckP()->d_outCriteria.empty() ||
      d_dataManager_p->getOutputDeckP()->d_dtOutOld == d_dataManager_p->getOutputDeckP()->d_dtOutCriteria)
    return;

  // perform checks every dt large intervals
  if (d_n % d_dataManager_p->getOutputDeckP()->d_dtOutOld != 0) return;

  if (d_dataManager_p->getOutputDeckP()->d_outCriteria == "max_Z" ||
      d_dataManager_p->getOutputDeckP()->d_outCriteria == "max_Z_stop") {
    // since damage function will not reduce once attaining desired maximum
    // value, we do not check if the criteria was met in the past
    if (d_dataManager_p->getOutputDeckP()->d_outCriteria == "max_Z" &&
        d_dataManager_p->getOutputDeckP()->d_dtOut == d_dataManager_p->getOutputDeckP()->d_dtOutCriteria)
      return;

    // change from large interval to small interval should be done only once
    // to do this, we check below flag which will be set to true in first
    // change to small from small interval
    static bool changed_to_small = false;
    bool changed_to_small_at_current = false;
    if (d_dataManager_p->getOutputDeckP()->d_dtOut > d_dataManager_p->getOutputDeckP()->d_dtOutCriteria &&
        !changed_to_small) {
      // get maximum from the damage data
      auto max = util::methods::max((*d_dataManager_p->getDamageFunctionP()));

      // check if it is desired range and change output frequency
      if (util::compare::definitelyGreaterThan(
              max, d_dataManager_p->getOutputDeckP()->d_outCriteriaParams[0])) {
        d_dataManager_p->getOutputDeckP()->d_dtOut = d_dataManager_p->getOutputDeckP()->d_dtOutCriteria;

        std::cout << "Message: Changing output interval to smaller value.\n";

        // dump this value
        std::ofstream fdump(d_dataManager_p->getOutputDeckP()->d_path + "dt_out_change.info");
        fdump << "Large_Dt_To_Small_Dt:\n";
        fdump << "  N: " << d_n << "\n";
        fdump << "  dN: " << d_n / d_dataManager_p->getOutputDeckP()->d_dtOutCriteria << "\n";
        fdump.close();

        changed_to_small = true;
        changed_to_small_at_current = true;
      }
    }  // if current dt out is larger

    // we now check (only if max_Z_stop flag is provided) if we need to
    // revert to large time interval
    // we do this if time interval is small at present
    // Also do not change to large interval back if just now we changed it to
    // small interval

    // change from small to large interval should be done only once
    // to do this, we check below flag which will be set to true in first
    // change to large from small interval
    static bool changed_back_to_large = false;
    if (!changed_to_small_at_current && !changed_back_to_large &&
        d_dataManager_p->getOutputDeckP()->d_outCriteria == "max_Z_stop" &&
        d_dataManager_p->getOutputDeckP()->d_dtOut < d_dataManager_p->getOutputDeckP()->d_dtOutOld) {
      // check maximum of Z function in the rectangle
      auto rect = std::make_pair(util::Point3(), util::Point3());
      auto ps = d_dataManager_p->getOutputDeckP()->d_outCriteriaParams;
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
      if (N > d_dataManager_p->getMeshP()->getNumNodes())
        N = d_dataManager_p->getMeshP()->getNumNodes();
      rect_ids.resize(N);

      auto f = hpx::parallel::for_loop(
          hpx::parallel::execution::par(hpx::parallel::execution::task), 0, N,
          [this, N, rect, refZ, &rect_ids](boost::uint64_t I) {
            size_t ibegin = I * N;
            size_t iend = (I + 1) * N;
            if (iend > d_dataManager_p->getMeshP()->getNumNodes())
              iend = d_dataManager_p->getMeshP()->getNumNodes();
            for (size_t i = ibegin; i < iend; i++) {
              if (util::geometry::isPointInsideRectangle(
                      d_dataManager_p->getMeshP()->getNode(i), rect.first,
                      rect.second))
                if (util::compare::definitelyGreaterThan(
                        (*d_dataManager_p->getDamageFunctionP())[i], refZ))
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
        d_dataManager_p->getOutputDeckP()->d_dtOut = d_dataManager_p->getOutputDeckP()->d_dtOutOld;

        std::cout << "Message: Changing output interval to larger value.\n";

        // dump this value
        std::ofstream fdump(d_dataManager_p->getOutputDeckP()->d_path + "dt_out_change.info",
                            std::ios::app);
        fdump << "Small_Dt_To_Large_Dt:\n";
        fdump << "  N: " << d_n << "\n";
        fdump << "  dN: " << d_n / d_dataManager_p->getOutputDeckP()->d_dtOutCriteria << "\n";
        fdump.close();

        changed_back_to_large = true;

        d_stop = true;
      }
    }  // if max_Z_stop and current dt out is smaller
  }    // if either max_Z or max_Z_stop criteria
}

template class model::FDModel<material::pd::RNPBond>;
