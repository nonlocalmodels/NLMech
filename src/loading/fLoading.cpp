////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "fLoading.h"
#include "../inp/decks/loadingDeck.h"
#include "fe/mesh.h"
#include "util/compare.h"
#include "util/utilFunction.h"
#include "util/utilGeom.h"

loading::FLoading::FLoading(inp::LoadingDeck *deck, fe::Mesh *mesh) {
  d_bcData = deck->d_fBCData;

  // fill the list of nodes where bc is applied and also set fixity of these
  // nodes
  for (const auto &bc : d_bcData) {
    // check bc first
    if (bc.d_regionType != "rectangle" and
        bc.d_regionType != "angled_rectangle") {
      std::cerr
          << "Error: Force bc region type = " << bc.d_regionType
          << " not recognised. Should be rectangle or angled_rectangle. \n";
      exit(1);
    }

    if (bc.d_spatialFnType != "constant" and bc.d_spatialFnType != "hat_x" and
        bc.d_spatialFnType != "hat_y" and bc.d_spatialFnType != "sin_x" and
        bc.d_spatialFnType != "sin_y" and bc.d_spatialFnType != "linear_x" and
        bc.d_spatialFnType != "linear_y") {
      std::cerr << "Error: Force bc space function type = "
                << bc.d_spatialFnType << " not recognised. "
                << "Currently only constant, hat_x, hat_y function is "
                   "implemented. \n";
      exit(1);
    }

    if (bc.d_timeFnType != "constant" and bc.d_timeFnType != "linear" and
        bc.d_timeFnType != "linear_step" and
        bc.d_timeFnType != "linear_slow_fast" and bc.d_timeFnType != "sin") {
      std::cerr << "Error: Force bc space function type = " << bc.d_timeFnType
                << " not recognised. "
                << "Currently constant, linear, linear_step, linear_slow_fast"
                   " are implemented. \n";
      exit(1);
    }

    size_t time_num_params = 1;
    if (bc.d_timeFnType == "linear_step")
      time_num_params = 3;
    else if (bc.d_timeFnType == "linear_slow_fast")
      time_num_params = 4;
    else if (bc.d_timeFnType == "sin")
      time_num_params = 2;
    if (bc.d_timeFnParams.size() != time_num_params) {
      std::cerr << "Error: Force bc insufficient parameters for time function. "
                   "Need "
                << time_num_params << " parameters but only "
                << bc.d_timeFnParams.size() << " parameters provided.\n";
      exit(1);
    }

    size_t spatial_num_params = 0;
    if (bc.d_spatialFnType == "hat_x" or bc.d_spatialFnType == "hat_y" or
        bc.d_spatialFnType == "sin_x" or bc.d_spatialFnType == "sin_y" or
        bc.d_spatialFnType == "linear_x" or bc.d_spatialFnType == "linear_y")
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
      if (bc.d_regionType == "rectangle" &&
          util::geometry::isPointInsideRectangle(mesh->getNode(i), bc.d_x1,
                                                 bc.d_x2, bc.d_y1, bc.d_y2))
        fix_nodes.push_back(i);
      else if (bc.d_regionType == "angled_rectangle" &&
               util::geometry::isPointInsideAngledRectangle(
                   mesh->getNode(i), bc.d_x1, bc.d_x2, bc.d_y1, bc.d_y2,
                   bc.d_theta))
        fix_nodes.push_back(i);
    }  // loop over nodes

    // add computed list of nodes to the data
    d_bcNodes.push_back(fix_nodes);
  }  // loop over bc sets
}

void loading::FLoading::apply(const double &time, std::vector<util::Point3> *f,
                              fe::Mesh *mesh) {
  for (size_t s = 0; s < d_bcData.size(); s++) {
    inp::BCData bc = d_bcData[s];
    for (auto i : d_bcNodes[s]) {
      util::Point3 x = mesh->getNode(i);
      double fmax = 1.0;

      // apply spatial function
      if (bc.d_spatialFnType == "hat_x") {
        // Hat function
        //
        //     f ^
        //       |
        //       |
        // f_max o
        //       |           /|\
        //       |         /  |  \
        //       |       /    |    \
        //       |     /      |      \
        //       |   /        |        \
        //       | /          |          \
	//       o____________o____________o______\ x
        //                                        /
        //    loc_x_min                 loc_x_max
        //
        fmax = bc.d_spatialFnParams[0] *
               util::function::hatFunction(x.d_x, bc.d_x1, bc.d_x2);
      } else if (bc.d_spatialFnType == "hat_y") {
        fmax = bc.d_spatialFnParams[0] *
               util::function::hatFunction(x.d_y, bc.d_y1, bc.d_y2);
      } else if (bc.d_spatialFnType == "sin_x") {
        double a = M_PI * bc.d_spatialFnParams[0];
        fmax = bc.d_spatialFnParams[0] * std::sin(a * x.d_x);
      } else if (bc.d_spatialFnType == "sin_y") {
        double a = M_PI * bc.d_spatialFnParams[0];
        fmax = bc.d_spatialFnParams[0] * std::sin(a * x.d_y);
      } else if (bc.d_spatialFnType == "linear_x") {
        double a = bc.d_spatialFnParams[0];
        fmax = bc.d_spatialFnParams[0] * a * x.d_x;
      } else if (bc.d_spatialFnType == "linear_y") {
        double a = bc.d_spatialFnParams[0];
        fmax = bc.d_spatialFnParams[0] * a * x.d_y;
      }

      // apply time function
      if (bc.d_timeFnType == "linear")
        fmax *= time;
      else if (bc.d_timeFnType == "linear_step")
        fmax *= util::function::linearStepFunc(time, bc.d_timeFnParams[1],
                                               bc.d_timeFnParams[2]);
      else if (bc.d_timeFnType == "linear_slow_fast") {
        if (util::compare::definitelyGreaterThan(time, bc.d_timeFnParams[1]))
          fmax *= bc.d_timeFnParams[3] * time;
        else
          fmax *= bc.d_timeFnParams[2] * time;
      } else if (bc.d_timeFnType == "sin") {
        double a = M_PI * bc.d_timeFnParams[1];
        fmax *= std::sin(a * time);
      }

      // multiply by the slope
      fmax *= bc.d_timeFnParams[0];

      for (auto d : bc.d_direction) {
        if (d == 1)
          (*f)[i].d_x += fmax;
        else if (d == 2)
          (*f)[i].d_y += fmax;
        else if (d == 3)
          (*f)[i].d_z += fmax;
      }
    }  // loop over nodes
  }    // loop over bc sets
}
