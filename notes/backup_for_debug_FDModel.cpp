// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "fDModel.h"
#include <fstream>
#include <inp/decks/restartDeck.h>

// utils
#include "../../rw/reader.h"
#include "../../rw/writer.h"
#include "../../util/fastMethods.h"
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

namespace {

const double C = 392.699;
const double beta = 1.52789e+08;
const double rbar = std::sqrt(0.5 / beta);
const double factor_bond_fracture = 10.;
const double horizon = 0.008;
const double h = 0.008/4.;
const double vol_ball = M_PI * horizon * horizon;

static std::vector<std::vector<size_t>> neighbors;
static std::vector<util::Point3> nodes;
static std::vector<std::vector<int>> fracture;

void create_nodes(std::vector<util::Point3> *nodes) {

  nodes->clear();

  for (size_t i=0; i<=size_t(0.1/h); i++)
    for (size_t j=0; j<=size_t(0.1/h); j++)
      nodes->push_back((util::Point3(i*h, j*h, 0.)));
}

void create_neighbors(const std::vector<util::Point3> *nds,
                      std::vector<std::vector<size_t>> *ngs) {

  //  std::vector<std::vector<size_t>> neighbors;
  ngs->resize(nds->size());
  for (size_t i = 0; i < nds->size(); i++) {

    std::vector<size_t> ns;
    for (size_t j = 0; j < nds->size(); j++) {

      if (i == j)
        continue;

      auto pij = (*nds)[i] - (*nds)[j];
      double r = std::sqrt(pij.d_x * pij.d_x + pij.d_y * pij.d_y);
      if (util::compare::definitelyLessThan(r, horizon))
        ns.push_back(j);
    }
    (*ngs)[i] = ns;
  }
}

void create_fracture(const std::vector<util::Point3> *nds,
                     const std::vector<std::vector<size_t>> *ngs,
                     std::vector<std::vector<int>> *fs){

  fs->resize(nds->size());
  for (size_t i=0; i < nds->size(); i++) {
    std::vector<int> a((*ngs)[i].size(), 0);

    for (size_t j=0; j<a.size(); j++) {

      auto pi = nds->at(i);
      auto pj = nds->at((*ngs)[i][j]);
      auto pji = pj - pi;

      bool intersect = false;
      if (util::compare::definitelyLessThan(pi.d_x, 0.05 + 1.0E-8) &&
          util::compare::definitelyLessThan(pj.d_x, 0.05 + 1.0E-8))
        continue;
      else if (util::compare::definitelyGreaterThan(pi.d_x, 0.05 + 1.0E-8) &&
          util::compare::definitelyGreaterThan(pj.d_x, 0.05 + 1.0E-8))
        continue;
      else {

        if (util::compare::definitelyLessThan(pi.d_y, 0.02 + 1.0E-8) &&
            util::compare::definitelyLessThan(pj.d_y, 0.02 + 1.0E-8))
          a[j] = 1;
      }
    }

    (*fs)[i] = a;
  }
}

double get_critical_bond_strain(double r) { return rbar / std::sqrt(r); }

double getInfFn(double r) { return  12. * (1. - r / horizon); }

std::pair<double, double>
getFandFprimeNew(double r, double s, int &bond_fractured, bool
break_bonds) {

  if (break_bonds && bond_fractured == 0)
    if (util::compare::definitelyGreaterThan(
        std::abs(s), factor_bond_fracture * get_critical_bond_strain(r)))
      bond_fractured = 1;

  double argument = -beta * r * s * s;

  if (bond_fractured == 0)
    return std::make_pair(C * (1 - std::exp(argument)) / (horizon * vol_ball),
                          C * beta * std::exp(argument) * 4.0 * s /
                              (r * horizon * vol_ball));
  else
    return std::make_pair(C / (horizon * vol_ball), 0.0);
}

void apply_u_bc(double time, std::vector<util::Point3> *u,
                std::vector<util::Point3> *v, const std::vector<util::Point3> *nds) {

  for (size_t i = 0; i < nds->size(); i++) {

    auto p = (*nds)[i];

    if (util::compare::definitelyGreaterThan(p.d_y, 0.092)) {
      (*u)[i].d_x = 0.;
      (*u)[i].d_y = 0.;

      (*v)[i].d_x = 0.;
      (*v)[i].d_y = 0.;
    }

    if (util::compare::definitelyLessThan(p.d_y, horizon)) {
      if (util::compare::definitelyLessThan(p.d_x, 0.05 + 1.0E-10)) {
        (*u)[i].d_x = -0.5 * time;

        (*v)[i].d_x = -0.5;
      } else {
        (*u)[i].d_x = 0.5 * time;

        (*v)[i].d_x = 0.5;
      }
    }
  }
}
} // namespace

model::FDModel::FDModel(inp::Input *deck)
    : d_massMatrix_p(nullptr), d_mesh_p(nullptr), d_fracture_p(nullptr),
      d_neighbor_p(nullptr), d_interiorFlags_p(nullptr), d_input_p(deck),
      d_policy_p(nullptr), d_initialCondition_p(nullptr), d_uLoading_p(nullptr),
      d_fLoading_p(nullptr), d_material_p(nullptr) {

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

  // read displacement and velocity from restart file
  rw::reader::readVtuFileRestart(d_restartDeck_p->d_file, &d_u, &d_v);

  // integrate in time
  integrate();
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
  if (d_policy_p->populateData("Model_d_hS"))
    d_hS = std::vector<double>(nnodes, 0.0);

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
  }
}

void model::FDModel::integrate() {

  create_nodes(&nodes);
  create_neighbors(&nodes, &neighbors);
  create_fracture(&nodes, &neighbors, &fracture);
  apply_u_bc(d_time, &d_u, &d_v, &nodes);

  // internal forces
  computeForces();

  // perform output at the beginning
  if (d_n == 0)
    output();

  // start time integration
  size_t i = d_n;
  for (i; i < d_modelDeck_p->d_Nt; i++) {

    integrateCD();

    if ((d_n % d_outputDeck_p->d_dtOut == 0) &&
        (d_n >= d_outputDeck_p->d_dtOut))
      output();
  } // loop over time steps
}

void model::FDModel::integrateCD() {

  double factor = 1.;
  if (d_n == 0)
    factor = 0.5;

  //   parallel for loop
//  auto f = hpx::parallel::for_loop(
//      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
//      d_mesh_p->getNumNodes(), [this, factor](boost::uint64_t i) {

  for (size_t i=0; i < nodes.size(); i++) {
    size_t dim = 2;                   // this->d_mesh_p->getDimension();
    double delta_t = 0.0001 / 25000.; // this->d_modelDeck_p->d_dt;
    double density = 1200.;           // this->d_material_p->getDensity();
    double fact = delta_t * delta_t / density;

    auto p = nodes[i];

    bool process_x = true;
    if (util::compare::definitelyGreaterThan(p.d_y, 0.092 - 1.0E-8))
      process_x = false;
    if (util::compare::definitelyLessThan(p.d_y, horizon + 1.0E-8))
      process_x = false;

    bool process_y = true;
    if (util::compare::definitelyGreaterThan(p.d_y, 0.092 - 1.0E-8))
      process_y = false;

    if (process_x) {

      double u_old = this->d_u[i].d_x;

      this->d_u[i].d_x +=
          factor * fact * this->d_f[i].d_x + delta_t * this->d_v[i].d_x;

      this->d_v[i].d_x = (this->d_u[i].d_x - u_old) / delta_t;
    }


    if (process_y) {

      double u_old = this->d_u[i].d_y;

      this->d_u[i].d_y +=
          factor * fact * this->d_f[i].d_y + delta_t * this->d_v[i].d_y;

      this->d_v[i].d_y = (this->d_u[i].d_y - u_old) / delta_t;
    }
  }

//); // end of parallel for loop
//
//  f.get();

  // compute forces and energy due to new displacement field (this will be
  // used in next time step)
  d_n += 1;
  d_time += d_modelDeck_p->d_dt;

  // boundary condition
  apply_u_bc(d_time, &d_u, &d_v, &nodes);

  // internal forces
  computeForces();
}

void model::FDModel::computeForces() {
//  auto f = hpx::parallel::for_loop(
//      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
//      d_mesh_p->getNumNodes(),
//      [this](boost::uint64_t i) {

  for (size_t i=0; i < nodes.size(); i++) {
    // local variable to hold force
    util::Point3 force_i = util::Point3();

    // reference coordinate and displacement at the node
    util::Point3 xi = nodes[i];
    util::Point3 ui = this->d_u[i];

    // upper and lower bound for volume correction
    double check_up = horizon + 0.5 * h;
    double check_low = horizon - 0.5 * h;

    // inner loop over neighbors
    for (size_t j = 0; j < neighbors[i].size(); j++) {
      size_t j_id = neighbors[i][j];

      // compute bond-based contribution
      util::Point3 xj = nodes[j_id];
      util::Point3 uj = this->d_u[j_id];
      util::Point3 xji = xj - xi;
      util::Point3 uji = uj - ui;
      double rji = xji.length();
      double Sji = xji.dot(uji) / xji.dot(xji);

      // get corrected volume of node j
      double volj = h * h;
      if (util::compare::definitelyGreaterThan(rji, check_low))
        volj *= (check_up - rji) / h;

      auto bond_fractured = fracture[i][j];
      auto ef = getFandFprimeNew(rji, Sji, bond_fractured, true);

      // update the fractured state of bond
      fracture[i][j] = bond_fractured;

      // compute influence function
      auto wji = getInfFn(rji);

      // compute the contribution of bond force to force at i
      double scalar_f = wji * ef.second * volj;

      if (bond_fractured == 0) {
        force_i.d_x += scalar_f * xji.d_x;
        force_i.d_y += scalar_f * xji.d_y;
        force_i.d_z += scalar_f * xji.d_z;
      }
    } // loop over neighboring nodes

    // update force and energy
    this->d_f[i] = force_i;
  } // loop over nodes

//  ); // end of parallel for loop
//
//  f.get();
}

void model::FDModel::integrateVerlet() {}

void model::FDModel::computeHydrostaticStrains() {}

void model::FDModel::computePostProcFields() {}

void model::FDModel::output() {

  std::cout << "Output: time step = " << this->d_n << "\n";

  //
  // List of data that we need to output
  //
  // Major simulation data
  // 1. Displacement (vector)
  // 2. Velocity (vector)
  // 3. Force (vector)
  // 4. time (single data)
  //
  // Minor simulation data (if these are computed)
  // 1. Force_Density (vector)
  // 2. Strain_Energy (scalar)
  // 3. Energy (scalar)
  // 4. Work_Done (scalar)
  // 5. Fixity (scalar)
  // 6. Node_Volume (scalar)
  // 7. Damage (scalar)
  // 8. Damage Z (scalar)
  // 9. Fracture_Perienergy_Bond (scalar)
  // 10. Fracture_Perienergy_Total (scalar) if StateFracture material
  // 11. Total energy (single data)
  // 12. Total fracture energy bond (single data)
  // 13. Total fracture energy total if StateFracture material (single data)
  //

  // filename
  std::string filename = d_outputDeck_p->d_path + "output_" +
      std::to_string(d_n / d_outputDeck_p->d_dtOut);

  // open
  auto writer = rw::writer::VtkWriterInterface(filename);

  // write mesh
//  writer.appendNodes(d_mesh_p->getNodesP(), &d_u);
  writer.appendNodes(&nodes, &d_u);

  //
  // major simulation data
  //
  std::string tag = "Displacement";
  writer.appendPointData(tag, &d_u);

  tag = "Velocity";
  writer.appendPointData(tag, &d_v);

  tag = "Force";
  if (d_outputDeck_p->isTagInOutput(tag)) {

//    std::vector<util::Point3> force(d_mesh_p->getNumNodes(), util::Point3());
    std::vector<util::Point3> force(nodes.size(), util::Point3());

    for (size_t i = 0; i < d_f.size(); i++)
      force[i] = d_f[i] * h * h;

    writer.appendPointData(tag, &force);
  }

  tag = "time";
  writer.addTimeStep(d_time);

  //
  // minor simulation data
  //
  //  if (!d_policy_p->enablePostProcessing()) {
  writer.close();
  return;
  //  }

  tag = "Force_Density";
  if (d_outputDeck_p->isTagInOutput(tag))
    writer.appendPointData(tag, &d_f);

  tag = "Strain_Energy";
  if (d_outputDeck_p->isTagInOutput(tag) &&
      d_policy_p->populateData("Model_d_e"))
    writer.appendPointData(tag, &d_e);

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
    writer.appendPointData(tag, d_mesh_p->getFixityP());

  tag = "Damage";
  if (d_outputDeck_p->isTagInOutput(tag) &&
      d_policy_p->populateData("Model_d_phi"))
    writer.appendPointData(tag, &d_phi);

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

  tag = "Total_Strain_Energy";
  if (d_outputDeck_p->isTagInOutput(tag) &&
      d_policy_p->populateData("Model_d_e"))
    writer.appendFieldData(tag, d_te);

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