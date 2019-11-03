////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "dampingGeom.h"
#include "fe/mesh.h"
#include "inp/decks/absborbingCondDeck.h"
#include "util/compare.h"
#include "util/utilGeom.h"
#include <iostream>

//
// DampingGeom
//
geometry::DampingGeom::DampingGeom(inp::AbsorbingCondDeck *deck,
const fe::Mesh *mesh)
    : d_dim(mesh->getDimension()), d_absorbingDeck_p(deck) {

  if (!d_absorbingDeck_p->d_dampingActive)
    return;

  d_coefficients.resize(mesh->getNumNodes());
  computeDampingCoefficient(mesh);
}

void geometry::DampingGeom::computeDampingCoefficient(
    const fe::Mesh *mesh) {

  d_coefficients = std::vector<double>(mesh->getNumNodes());

  for (size_t i=0; i<mesh->getNumNodes(); i++) {

    auto x = mesh->getNode(i);

    // find which damping region this node belongs to
    int loc_r = -1;
    for (size_t r = 0; r < d_absorbingDeck_p->d_dampingGeoms.size(); r++)
      if (util::geometry::isPointInsideCuboid(
              d_dim, x, d_absorbingDeck_p->d_dampingGeoms[r].d_p1,
              d_absorbingDeck_p->d_dampingGeoms[r].d_p2))
        loc_r = r;

    if (loc_r == -1)
      continue;

    auto dg = d_absorbingDeck_p->d_dampingGeoms[loc_r];
    double coeff = 0.;
    double exponent = d_absorbingDeck_p->d_dampingCoeffParams[0];

    // check if x coordinate is within the region
    if (util::compare::definitelyLessThan(x.d_x, dg.d_p2.d_x) &&
    util::compare::definitelyGreaterThan(x.d_x, dg.d_p1.d_x)) {
      if (dg.d_isLeft)
        coeff = std::pow((dg.d_p2.d_x - x.d_x) / (dg.d_p2.d_x - dg.d_p1.d_x),
                         exponent);
      else
        coeff = std::pow((x.d_x - dg.d_p1.d_x) / (dg.d_p2.d_x - dg.d_p1.d_x),
                         exponent);
    }

    if (d_dim > 1) {

      // check y
      if (util::compare::definitelyLessThan(x.d_y, dg.d_p2.d_y) &&
          util::compare::definitelyGreaterThan(x.d_y, dg.d_p1.d_y)) {
        bool check = false;
        if (d_dim == 2)
          check = dg.d_isBottom;
        else if (d_dim == 3)
          check = dg.d_isBack;

        if (check)
          coeff *= std::pow((dg.d_p2.d_y - x.d_y) / (dg.d_p2.d_y - dg.d_p1.d_y),
                           exponent);
        else
          coeff *= std::pow((x.d_y - dg.d_p1.d_y) / (dg.d_p2.d_y - dg.d_p1.d_y),
                           exponent);
      }
    }

    if (d_dim > 2) {

      // check z
      if (util::compare::definitelyLessThan(x.d_z, dg.d_p2.d_z) &&
          util::compare::definitelyGreaterThan(x.d_z, dg.d_p1.d_z)) {

        if (dg.d_isBottom)
          coeff *= std::pow((dg.d_p2.d_z - x.d_z) / (dg.d_p2.d_z - dg.d_p1.d_z),
                            exponent);
        else
          coeff *= std::pow((x.d_z - dg.d_p1.d_z) / (dg.d_p2.d_z - dg.d_p1.d_z),
                            exponent);
      }
    }

    d_coefficients[i] = coeff;
  }
}

double geometry::DampingGeom::getCoefficient(const size_t &i) { return
d_coefficients[i]; }
double geometry::DampingGeom::getCoefficient(const size_t &i) const { return
      d_coefficients[i]; }

bool geometry::DampingGeom::isDampingActive() { return d_absorbingDeck_p->d_dampingActive; }
bool geometry::DampingGeom::isDampingActive() const { return
d_absorbingDeck_p->d_dampingActive; }

bool geometry::DampingGeom::isViscousDamping() { return
      d_absorbingDeck_p->d_isViscousDamping; }
bool geometry::DampingGeom::isViscousDamping() const { return
d_absorbingDeck_p->d_isViscousDamping; }
