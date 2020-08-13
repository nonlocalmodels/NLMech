////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "materialUtil.h"
#include "util/compare.h"
#include <iostream>

// hpx lib
#include <hpx/include/parallel_algorithm.hpp>

namespace {

double computeStateMxI(size_t i, const std::vector<util::Point3> &nodes,
                       const std::vector<double> &nodal_vol,
                       const std::vector<std::vector<size_t>> &neighbor_list,
                       const double &mesh_size,
                       material::pd::Material *material) {

  double horizon = material->getHorizon();
  const auto &xi = nodes[i];
  double m = 0.;

  // upper and lower bound for volume correction
  auto check_up = horizon + 0.5 * mesh_size;
  auto check_low = horizon - 0.5 * mesh_size;

  const auto &i_neighs = neighbor_list[i];
  for (size_t j = 0; j < i_neighs.size(); j++) {

    auto j_id = i_neighs[j];

    const auto &xj = nodes[j_id];
    double rji = (xj - xi).length();

    if (util::compare::definitelyGreaterThan(rji, horizon) or j_id == i)
      continue;

    // get corrected volume of node j
    auto volj = nodal_vol[j_id];

    if (util::compare::definitelyGreaterThan(rji, check_low))
      volj *= (check_up - rji) / mesh_size;

    m += std::pow(rji, 2) * material->getInfFn(rji) * volj;
  }

  if (util::compare::definitelyLessThan(m, 1.0E-18)) {
    std::ostringstream oss;
    oss << "Error: Weighted nodal volume = " << m
              << " should not be too close to zero.\n";

    oss << "Mesh size = " << mesh_size << "\n";
    oss << "J = " << material->getInfFn(1.) << "\n";
    oss << material->printStr(0, 0);

    std::cout << oss.str();

    exit(1);
  }
  return m;
}

double computeStateThetaxI(size_t i, const std::vector<util::Point3> &nodes,
                           const std::vector<util::Point3> &nodes_disp,
                           const std::vector<double> &nodal_vol,
                           const std::vector<std::vector<size_t>> &neighbor_list,
                           const double &mesh_size,
                           material::pd::Material *material,
                           const geometry::Fracture *fracture,
                           const std::vector<double> &mx) {

  double horizon = material->getHorizon();
  const auto &xi = nodes[i];
  const auto &ui = nodes_disp[i];
  double m = mx[i];
  double theta = 0.;

  // upper and lower bound for volume correction
  auto check_up = horizon + 0.5 * mesh_size;
  auto check_low = horizon - 0.5 * mesh_size;

  const auto &i_neighs = neighbor_list[i];
  for (size_t j = 0; j < i_neighs.size(); j++) {

    auto j_id = i_neighs[j];

    const auto &xj = nodes[j_id];
    const auto &uj = nodes_disp[j_id];
    double rji = (xj - xi).length();

    if (util::compare::definitelyGreaterThan(rji, horizon) or j_id == i)
      continue;

    // get corrected volume of node j
    auto volj = nodal_vol[j_id];

    if (util::compare::definitelyGreaterThan(rji, check_low))
      volj *= (check_up - rji) / mesh_size;

    // get bond state
    double bond_state = fracture->getBondState(i, j) ? 0. : 1.;

    // get change in bond length
    auto yi = xi + ui;
    auto yj = xj + uj;
    double change_length = (yj - yi).length() - rji;

    theta += bond_state * rji * change_length *
             material->getInfFn(rji) * volj;
  }

  return 3. * theta / m;
}

double computeHydrostaticStrainI(size_t i, const std::vector<util::Point3> &nodes,
                           const std::vector<util::Point3> &nodes_disp,
                           const std::vector<double> &nodal_vol,
                                 const std::vector<std::vector<size_t>> &neighbor_list,
                           const double &mesh_size,
                           material::pd::Material *material,
                           const geometry::Fracture *fracture,
                           size_t dim) {

  double horizon = material->getHorizon();
  const auto &xi = nodes[i];
  const auto &ui = nodes_disp[i];
  double theta = 0.;

  // upper and lower bound for volume correction
  auto check_up = horizon + 0.5 * mesh_size;
  auto check_low = horizon - 0.5 * mesh_size;

  // get volume of ball
  double vol_ball = std::pow(horizon, 2) * M_PI;
  if (dim == 3)
    vol_ball *= horizon * 4. / 3.;

  const auto &i_neighs = neighbor_list[i];
  for (size_t j = 0; j < i_neighs.size(); j++) {

    auto j_id = i_neighs[j];

    const auto &xj = nodes[j_id];
    const auto &uj = nodes_disp[j_id];
    double rji = (xj - xi).length();

    if (util::compare::definitelyGreaterThan(rji, horizon) or j_id == i)
      continue;

    // get corrected volume of node j
    auto volj = nodal_vol[j_id];

    if (util::compare::definitelyGreaterThan(rji, check_low))
      volj *= (check_up - rji) / mesh_size;

    // get bond state
    double bond_state = fracture->getBondState(i, j) ? 0. : 1.;

    // get bond strain
    double Sji = material->getS(xj - xi, uj - ui);

    theta += bond_state * rji * Sji *
             material->getInfFn(rji) * volj / vol_ball;
  }

  return theta;
}

void updateBondFractureDataI(size_t i, const std::vector<util::Point3> &nodes,
                                 const std::vector<util::Point3> &nodes_disp,
                             const std::vector<std::vector<size_t>> &neighbor_list,
                                 material::pd::Material *material,
                                 geometry::Fracture *fracture) {

  const auto &i_neighs = neighbor_list[i];
  for (size_t j = 0; j < i_neighs.size(); j++) {

    auto j_id = i_neighs[j];

    double s = material->getS(nodes[j_id] - nodes[i], nodes_disp[j_id] -
    nodes_disp[i]);
    double sc = material->getSc((nodes[j_id] - nodes[i]).length());

    // get fracture state, modify, and set
    auto fs = fracture->getBondState(i, j);
    if (!fs &&
        util::compare::definitelyGreaterThan(std::abs(s), sc + 1.0e-10)) {
      fs = true;
      fracture->setBondState(i, j, fs);
    }
  }
}

} // anonymous namespace

void material::computeStateMx(const std::vector<util::Point3> &nodes,
                              const std::vector<double> &nodal_vol,
                              const std::vector<std::vector<size_t>> &neighbor_list,
                              const double &mesh_size,
                              material::pd::Material *material,
                              std::vector<double> &mx,
                              bool compute_in_parallel) {

  mx.resize(nodes.size());
  if (!compute_in_parallel) {
    for (size_t i = 0; i < nodes.size(); i++)
      mx[i] = computeStateMxI(i, nodes, nodal_vol, neighbor_list, mesh_size, material);
  } else {

    auto f = hpx::parallel::for_loop(
        hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
        nodes.size(),
        [nodes, nodal_vol, neighbor_list, mesh_size, material, &mx](boost::uint64_t i) {
          mx[i] = computeStateMxI(i, nodes, nodal_vol, neighbor_list, mesh_size,
              material);
        } // loop over nodes
    );    // end of parallel for loop
    f.get();
  }
}

void material::computeStateThetax(
    const std::vector<util::Point3> &nodes,
    const std::vector<util::Point3> &nodes_disp,
    const std::vector<double> &nodal_vol,
    const std::vector<std::vector<size_t>> &neighbor_list,
    const double &mesh_size,
    material::pd::Material *material,
    const geometry::Fracture *fracture, const std::vector<double> &mx,
    std::vector<double> &thetax, bool compute_in_parallel) {

  thetax.resize(nodes.size());
  if (!compute_in_parallel) {
    for (size_t i = 0; i < nodes.size(); i++)
      thetax[i] =
          computeStateThetaxI(i, nodes, nodes_disp, nodal_vol, neighbor_list,
                              mesh_size, material, fracture, mx);
  } else {

    auto f = hpx::parallel::for_loop(
        hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
        nodes.size(),
        [nodes, nodes_disp, nodal_vol, neighbor_list, mesh_size, material, fracture, mx,
         &thetax](boost::uint64_t i) {
          thetax[i] = computeStateThetaxI(i, nodes, nodes_disp, nodal_vol,
                                          neighbor_list, mesh_size, material,
                                          fracture, mx);
        } // loop over nodes
    );    // end of parallel for loop
    f.get();
  }
}

void material::computeHydrostaticStrain(const std::vector<util::Point3> &nodes,
                               const std::vector<util::Point3> &nodes_disp,
                               const std::vector<double> &nodal_vol,
                                        const std::vector<std::vector<size_t>> &neighbor_list,
                               const double &mesh_size,
                               material::pd::Material *material,
                               const geometry::Fracture *fracture,
                               std::vector<double> &thetax, size_t dim,
                               bool compute_in_parallel) {

  thetax.resize(nodes.size());
  if (!compute_in_parallel) {
    for (size_t i = 0; i < nodes.size(); i++)
      thetax[i] = computeHydrostaticStrainI(i, nodes, nodes_disp, nodal_vol,
                                            neighbor_list, mesh_size, material,
                                            fracture, dim);
  } else {

    auto f = hpx::parallel::for_loop(
        hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
        nodes.size(),
        [nodes, nodes_disp, nodal_vol, neighbor_list, mesh_size, material, fracture, dim,
         &thetax](boost::uint64_t i) {
          thetax[i] = computeHydrostaticStrainI(i, nodes, nodes_disp, nodal_vol,
                                                neighbor_list, mesh_size,
                                                material, fracture, dim);
        } // loop over nodes
    );    // end of parallel for loop
    f.get();
  }

}

void material::updateBondFractureData(const std::vector<util::Point3> &nodes,
                            const std::vector<util::Point3> &nodes_disp,
                                      const std::vector<std::vector<size_t>> &neighbor_list,
                            material::pd::Material *material,
                            geometry::Fracture *fracture, bool
                            compute_in_parallel) {

  if (!compute_in_parallel) {
    for (size_t i = 0; i < nodes.size(); i++)
      updateBondFractureDataI(i, nodes, nodes_disp, neighbor_list, material, fracture);
  } else {

    auto f = hpx::parallel::for_loop(
        hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
        nodes.size(),
        [nodes, nodes_disp, neighbor_list, material, &fracture](boost::uint64_t i) {
          updateBondFractureDataI(i, nodes, nodes_disp, neighbor_list, material,
              fracture);
        } // loop over nodes
    );    // end of parallel for loop
    f.get();
  }
}