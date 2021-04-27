////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "uLoading.h"

#include "../inp/decks/loadingDeck.h"
#include "fe/mesh.h"
#include "util/compare.h"
#include "util/utilFunction.h"
#include "util/utilGeom.h"
#include "util/utilIO.h"

loading::ULoading::ULoading(inp::LoadingDeck *deck, fe::Mesh *mesh) {
  d_bcData = deck->d_uBCData;

  // fill the list of nodes where bc is applied and also set fixity of these
  // nodes

  
  for (const auto &bc : d_bcData) {
    // check bc first
   
    if (bc.d_regionType != "rectangle" and
        bc.d_regionType != "angled_rectangle" and
        bc.d_regionType != "circle" and bc.d_regionType != "torus" and
        bc.d_regionType != "line" and bc.d_regionType != "displacement_from_pum")  {
      std::cerr << "Error: Displacement bc region type = " << bc.d_regionType
                << " not recognised. Should be rectangle or angled_rectangle "
                   "or circle or torus or displacement_from_pum. \n";
      exit(1);
    }

    if (bc.d_spatialFnType != "constant" and bc.d_spatialFnType != "sin_x" and
        bc.d_spatialFnType != "sin_y" and bc.d_spatialFnType != "sin_z" and
        bc.d_spatialFnType != "linear_x" and
        bc.d_spatialFnType != "linear_y" and bc.d_spatialFnType != "linear_z") {
      std::cerr << "Error: Displacement bc space function type = "
                << bc.d_spatialFnType << " not recognized. "
                << "Currently only constant function is implemented. \n";
      exit(1);
    }

    if (bc.d_timeFnType != "constant" and bc.d_timeFnType != "linear" and
        bc.d_timeFnType != "quadratic" and bc.d_timeFnType != "linear_step" and
        bc.d_timeFnType != "linear_slow_fast" and bc.d_timeFnType != "sin") {
      std::cerr << "Error: Displacement bc space function type = "
                << bc.d_timeFnType << " not recognised. "
                << "Currently constant, linear, quadratic, sin functions are "
                   "implemented. \n";
      exit(1);
    }

    size_t time_num_params = 1;
    if (bc.d_timeFnType == "quadratic" or bc.d_timeFnType == "sin")
      time_num_params = 2;
    if (bc.d_timeFnType == "linear_step")
      time_num_params = 3;
    else if (bc.d_timeFnType == "linear_slow_fast")
      time_num_params = 4;
    if (bc.d_timeFnParams.size() != time_num_params) {
      std::cerr << "Error: Displacement bc insufficient parameters for time "
                   "function. Need "
                << time_num_params << " parameters but only "
                << bc.d_timeFnParams.size() << " parameters provided.\n";
      exit(1);
    }

    size_t spatial_num_params = 0;
    if (bc.d_spatialFnType == "sin_x" or bc.d_spatialFnType == "sin_y" or
        bc.d_spatialFnType == "sin_z" or bc.d_spatialFnType == "linear_x" or
        bc.d_spatialFnType == "linear_y" or bc.d_spatialFnType == "linear_y")
      spatial_num_params = 1;
    if (bc.d_spatialFnParams.size() < spatial_num_params) {
      std::cerr << "Error: Force bc insufficient parameters for spatial "
                   "function. Need "
                << spatial_num_params << " parameters but only "
                << bc.d_spatialFnParams.size() << " parameters provided.\n";
      exit(1);
    }

    // compute list of nodes which are marked fixed
    std::vector<size_t> fix_nodes;

    // now loop over nodes
    for (size_t i = 0; i < mesh->getNumNodes(); i++) {
      bool node_fixed = false;

      if (bc.d_regionType == "rectangle" &&
          util::geometry::isPointInsideRectangle(mesh->getNode(i), bc.d_x1,
                                                 bc.d_x2, bc.d_y1, bc.d_y2)) {
        node_fixed = true;

        // loop over direction and set fixity
        for (auto dof : bc.d_direction) {
          // pass 0 for x, 1 for y, and 2 for z dof
          mesh->setFixity(i, dof - 1, true);
        }

      } else if (bc.d_regionType == "angled_rectangle" &&
                 util::geometry::isPointInsideAngledRectangle(
                     mesh->getNode(i), bc.d_x1, bc.d_x2, bc.d_y1, bc.d_y2,
                     bc.d_theta)) {
        node_fixed = true;

        // loop over direction and set fixity
        for (auto dof : bc.d_direction) {
          // pass 0 for x, 1 for y, and 2 for z dof
          mesh->setFixity(i, dof - 1, true);
        }

      } else if (bc.d_regionType == "circle" &&
                 util::geometry::isPointinCircle(
                     mesh->getNode(i), util::Point3(bc.d_x1, bc.d_y1, 0.),
                     bc.d_r1)) {
        node_fixed = true;

        // loop over direction and set fixity
        for (auto dof : bc.d_direction) {
          // pass 0 for x, 1 for y, and 2 for z dof
          mesh->setFixity(i, dof - 1, true);
        }
      } else if (bc.d_regionType == "torus" &&
                 util::geometry::isPointinCircle(
                     mesh->getNode(i), util::Point3(bc.d_x1, bc.d_y1, 0.),
                     bc.d_r1) &&
                 !util::geometry::isPointinCircle(
                     mesh->getNode(i), util::Point3(bc.d_x1, bc.d_y1, 0.),
                     bc.d_r2)) {
        node_fixed = true;

        // loop over direction and set fixity
        for (auto dof : bc.d_direction) {
          // pass 0 for x, 1 for y, and 2 for z dof
          mesh->setFixity(i, dof - 1, true);
        }
      } else if (bc.d_regionType == "line" &&
                 util::compare::definitelyGreaterThan(mesh->getNode(i).d_x,
                                                      bc.d_x1) &&
                 util::compare::definitelyLessThan(mesh->getNode(i).d_x,
                                                   bc.d_x2)) {
        node_fixed = true;

        // loop over direction and set fixity
        for (auto dof : bc.d_direction) {
          // pass 0 for x, 1 for y, and 2 for z dof
          mesh->setFixity(i, dof - 1, true);
        }
      }
      else if (bc.d_regionType == "displacement_from_pum"){

        if (mesh->getPrescribedNodes()[i]==1)
          node_fixed = true;
      }

      // store the id of this node
      if (node_fixed) fix_nodes.push_back(i);
    }  // loop over nodes

    // add computed list of nodes to the data
    d_bcNodes.push_back(fix_nodes);
  }  // loop over bc sets
}

void loading::ULoading::apply(const double &time, std::vector<util::Point3> *u,
                              std::vector<util::Point3> *v, fe::Mesh *mesh) {
  for (size_t s = 0; s < d_bcData.size(); s++) {
    inp::BCData bc = d_bcData[s];
    for (auto i : d_bcNodes[s]) {
      util::Point3 x = mesh->getNode(i);
      double umax = bc.d_timeFnParams[0];
      double du = 0.;
      double dv = 0.;

      if (bc.d_regionType == "displacement_from_pum")
      {

        if (bc.d_direction.size() != 1)
        {
          std::cerr << "The region type: displacement_from_pum support only one direction per set. One set per direction is needed if multiple directions are used for the coupling." << std::endl;
          exit(1);
        }
        umax = mesh->getPrescribedValues()[i][bc.d_direction[0]-1] / 0.001 ;
        //std::cout << i << " " << umax << " " << bc.d_direction[0] << std::endl;

      }

      // apply spatial function
      if (bc.d_spatialFnType == "sin_x") {
        double a = M_PI * bc.d_spatialFnParams[0];
        umax = umax * std::sin(a * x.d_x);
      } else if (bc.d_spatialFnType == "sin_y") {
        double a = M_PI * bc.d_spatialFnParams[0];
        umax = umax * std::sin(a * x.d_y);
      } else if (bc.d_spatialFnType == "sin_z") {
        double a = M_PI * bc.d_spatialFnParams[0];
        umax = umax * std::sin(a * x.d_z);
      } else if (bc.d_spatialFnType == "linear_x") {
        double a = bc.d_spatialFnParams[0];
        umax = umax * a * x.d_x;
      } else if (bc.d_spatialFnType == "linear_y") {
        double a = bc.d_spatialFnParams[0];
        umax = umax * a * x.d_y;
      } else if (bc.d_spatialFnType == "linear_z") {
        double a = bc.d_spatialFnParams[0];
        umax = umax * a * x.d_z;
      }

      // apply time function
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
      } else if (bc.d_timeFnType == "linear_step") {
        du = umax * util::function::linearStepFunc(time, bc.d_timeFnParams[1],
                                                   bc.d_timeFnParams[2]);
        dv = umax * util::function::derLinearStepFunc(
                        time, bc.d_timeFnParams[1], bc.d_timeFnParams[2]);
      } else if (bc.d_timeFnType == "linear_slow_fast") {
        if (util::compare::definitelyGreaterThan(time, bc.d_timeFnParams[1])) {
          du = umax * bc.d_timeFnParams[3] * time;
          dv = umax * bc.d_timeFnParams[3];
        } else {
          du = umax * bc.d_timeFnParams[2] * time;
          dv = umax * bc.d_timeFnParams[2];
        }
      }

      for (auto d : bc.d_direction) {
        if (d == 1) {
          (*u)[i].d_x = du;
          (*v)[i].d_x = dv;
        } else if (d == 2) {
          (*u)[i].d_y = du ;
          (*v)[i].d_y = dv ;
        } else if (d == 3) {
          (*u)[i].d_z = du;
          (*v)[i].d_z = dv;
        }
      }
    }  // loop over nodes
  }    // loop over bc sets
}

std::string loading::ULoading::printStr(int nt, int lvl) const {
  auto tabS = util::io::getTabS(nt);
  std::ostringstream oss;
  oss << tabS << "------- FLoading --------" << std::endl << std::endl;
  oss << tabS << "Number of loading data = " << d_bcData.size() << std::endl;
  for (size_t i = 0; i < d_bcData.size(); i++) {
    oss << tabS << "Loading data " << i + 1 << " information" << std::endl;
    oss << d_bcData[i].printStr(nt + 1) << std::endl;
    oss << "Number of nodes affected = " << d_bcNodes[i].size() << std::endl;
  }
  oss << tabS << std::endl;

  return oss.str();
}
