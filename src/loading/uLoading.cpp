// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "uLoading.h"
#include "../inp/decks/loadingDeck.h"
#include "fe/mesh.h"
#include "util/utilGeom.h"

loading::ULoading::ULoading(inp::LoadingDeck *deck, fe::Mesh *mesh) {

  d_bcData = deck->d_uBCData;

  // fill the list of nodes where bc is applied and also set fixity of these
  // nodes
  for (const auto &bc : d_bcData) {

    // check bc first
    if (bc.d_regionType != "rectangle" and
        bc.d_regionType != "angled_rectangle") {
      std::cerr
          << "Error: Displacement bc region type = " << bc.d_regionType
          << " not recognised. Should be rectangle or angled_rectangle. \n";
      exit(1);
    }

    if (bc.d_spatialFnType != "constant") {
      std::cerr << "Error: Displacement bc space function type = "
                << bc.d_spatialFnType << " not recognised. "
                << "Currently only constant function is implemented. \n";
      exit(1);
    }

    if (bc.d_timeFnType != "constant" and bc.d_timeFnType != "linear" and
        bc.d_timeFnType != "quadratic" and bc.d_timeFnType != "sin") {

      std::cerr << "Error: Displacement bc space function type = "
                << bc.d_timeFnType << " not recognised. "
                << "Currently constant, linear, quadratic, sin functions are "
                   "implemented. \n";
      exit(1);
    }

    size_t time_num_params = 1;
    if (bc.d_timeFnType == "quadratic" or bc.d_timeFnType == "sin")
      time_num_params = 2;
    if (bc.d_timeFnParams.size() != time_num_params) {
      std::cerr << "Error: Displacement bc insufficient parameters for time "
                   "function. Need "
                << time_num_params << " parameters but only "
                << bc.d_timeFnParams.size() << " parameters provided.\n";
      exit(1);
    }

    // compute list of nodes which are marked fixed
    std::vector<size_t> fix_nodes;

    // now loop over nodes
    for (size_t i = 0; i < mesh->getNumNodes(); i++) {

      if (bc.d_regionType == "rectangle" &&
          util::geometry::isPointInsideRectangle(mesh->getNode(i), bc.d_x1,
                                                 bc.d_x2, bc.d_y1, bc.d_y2)) {
        // loop over direction and set fixity
        for (auto dof : bc.d_direction) {
          // pass 0 for x, 1 for y, and 2 for z dof
          mesh->setFixity(i, dof - 1, true);
        }

        // store the id of this node
        fix_nodes.push_back(i);
      } else if (bc.d_regionType == "angled_rectangle" &&
                 util::geometry::isPointInsideAngledRectangle(
                     mesh->getNode(i), bc.d_x1, bc.d_x2, bc.d_y1, bc.d_y2,
                     bc.d_theta)) {
        // loop over direction and set fixity
        for (auto dof : bc.d_direction) {
          // pass 0 for x, 1 for y, and 2 for z dof
          mesh->setFixity(i, dof - 1, true);
        }
        // store the id of this node
        fix_nodes.push_back(i);
      }

    } // loop over nodes

    // add computed list of nodes to the data
    d_bcNodes.push_back(fix_nodes);
  } // loop over bc sets
}

void loading::ULoading::apply(const double &time, std::vector<util::Point3> *u,
                              std::vector<util::Point3> *v, fe::Mesh *mesh) {

  for (size_t s = 0; s < d_bcData.size(); s++) {

    inp::BCData bc = d_bcData[s];
    for (auto i : d_bcNodes[s]) {

      double umax = bc.d_timeFnParams[0];
      double du = 0.;
      double dv = 0.;

      if (bc.d_timeFnType == "constant")
        du = umax;
      else if (bc.d_timeFnType == "linear") {
        du = umax * time;
        dv = umax;
      } else if (bc.d_timeFnType == "quadratic") {
        du = umax * time + bc.d_timeFnParams[1] * time * time;
        dv = umax + bc.d_timeFnParams[1] * time;
      } else if (bc.d_timeFnType == "sin") {
        double a = M_PI * bc.d_timeFnParams[1];
        du = umax * std::sin(a * time);
        dv = umax * a * std::cos(a * time);
      }

      for (auto d : bc.d_direction) {
        if (d == 1) {
          (*u)[i].d_x = du;
          (*v)[i].d_x = dv;
        } else if (d == 2) {
          (*u)[i].d_y = du;
          (*v)[i].d_y = dv;
        } else if (d == 3) {
          (*u)[i].d_z = du;
          (*v)[i].d_z = dv;
        }
      }
    } // loop over nodes
  }   // loop over bc sets
}
