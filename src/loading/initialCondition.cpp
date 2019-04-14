// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "initialCondition.h"
#include "fe/mesh.h"
#include "util/utilFunction.h"
#include "../inp/decks/initialConditionDeck.h"
#include <iostream>

loading::InitialCondition::InitialCondition(inp::InitialConditionDeck *deck)
:d_deck_p(deck) {}

void loading::InitialCondition::apply(std::vector<util::Point3> *u,
                                      std::vector<util::Point3> *v,
                                      fe::Mesh *mesh) {

  // process displacement
  if (!d_deck_p->d_uICData.d_type.empty() && d_deck_p->d_uICData.d_type !=
  "file") {
    for (size_t i=0; i<mesh->getNumNodes(); i++) {
      if (mesh->isNodeFree(i, 0))
        (*u)[i].d_x = getICFormula(d_deck_p->d_uICData.d_type,
                                   d_deck_p->d_uICData.d_params,
                                   mesh->getNode(i), 0, mesh->getDimension());
      if (mesh->getDimension() > 1)
        if (mesh->isNodeFree(i, 1))
          (*u)[i].d_y = getICFormula(d_deck_p->d_uICData.d_type,
                                   d_deck_p->d_uICData.d_params,
                                   mesh->getNode(i), 1, mesh->getDimension());

      if (mesh->getDimension() > 2)
        if (mesh->isNodeFree(i, 2))
          (*u)[i].d_z = getICFormula(d_deck_p->d_uICData.d_type,
                                     d_deck_p->d_uICData.d_params,
                                     mesh->getNode(i), 2, mesh->getDimension());
    }
  }

  // process velocity
  if (!d_deck_p->d_vICData.d_type.empty() && d_deck_p->d_vICData.d_type !=
      "file") {
    for (size_t i=0; i<mesh->getNumNodes(); i++) {
      if (mesh->isNodeFree(i, 0))
        (*v)[i].d_x = getICFormula(d_deck_p->d_vICData.d_type,
                                   d_deck_p->d_vICData.d_params,
                                   mesh->getNode(i), 0, mesh->getDimension());
      if (mesh->getDimension() > 1)
        if (mesh->isNodeFree(i, 1))
          (*v)[i].d_y = getICFormula(d_deck_p->d_vICData.d_type,
                                     d_deck_p->d_vICData.d_params,
                                     mesh->getNode(i), 1, mesh->getDimension());

      if (mesh->getDimension() > 2)
        if (mesh->isNodeFree(i, 2))
          (*v)[i].d_z = getICFormula(d_deck_p->d_vICData.d_type,
                                     d_deck_p->d_vICData.d_params,
                                     mesh->getNode(i), 2, mesh->getDimension());
    }
  }
}

double loading::InitialCondition::getICFormula(
    const std::string &fn_type, const std::vector<double> &params,
    const util::Point3 &x, const size_t &dof, const size_t &dim) {

  if (fn_type == "gaussian" && dim == 2)
    return util::function::gaussian2d(x, dof, params);
  else if (fn_type == "double_gaussian" && dim == 2)
    return util::function::doubleGaussian2d(x, dof, params);
  else {

    std::cerr << "Error: Check initial condition flag = "
              << fn_type << ". Currently only guassian and double_guassian "
                            "type functions are supported in 2-d.\n";
    exit(1);
  }
}
