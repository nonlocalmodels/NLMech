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

static double damping_geom_tol = 1.0e-8;

//
// DampingGeom
//
geometry::DampingGeom::DampingGeom(inp::AbsorbingCondDeck *deck,
                                   const fe::Mesh *mesh)
    : d_dim(mesh->getDimension()), d_absorbingDeck_p(deck) {
  if (!d_absorbingDeck_p->d_dampingActive) return;

  d_coefficients.resize(mesh->getNumNodes());
  computeDampingCoefficient(mesh);
}

void geometry::DampingGeom::computeDampingCoefficient(const fe::Mesh *mesh) {
  d_coefficients = std::vector<double>(mesh->getNumNodes(), 0.);

  for (size_t i = 0; i < mesh->getNumNodes(); i++) {
    auto x = mesh->getNode(i);

    // find which damping region this node belongs to
    int loc_r = -1;
    for (size_t r = 0; r < d_absorbingDeck_p->d_dampingGeoms.size(); r++)
      if (util::geometry::isPointInsideCuboid(
              d_dim, x, d_absorbingDeck_p->d_dampingGeoms[r].d_p1,
              d_absorbingDeck_p->d_dampingGeoms[r].d_p2))
        loc_r = r;

    if (loc_r == -1) continue;

    auto dg = d_absorbingDeck_p->d_dampingGeoms[loc_r];
    double coeff = 1.;
    double exponent = d_absorbingDeck_p->d_dampingCoeffParams[0];

    double check = 0.;
    double lower = 0.;
    double upper = 0.;
    double thickness = 0.;

    if (d_dim > 0 && dg.d_checkX) {
      thickness = dg.d_layerThicknessX;
      check = x.d_x;

      // check if x coordinate is within the region
      if (dg.d_relativeLoc == "left") {
        // special care when point is exactly on the upper line
        lower = dg.d_p2.d_x - thickness;
        upper = dg.d_p2.d_x;
        if (util::compare::definitelyLessThan(check, upper - damping_geom_tol))
          coeff *= std::pow((upper - check) / thickness, exponent);
        else
          coeff *= 0.;

        // std::cout << "Left layer, coeff = " << coeff << ", y = " << check <<
        // ", lower = " << lower << ", upper = " << upper << std::endl;
      } else if (dg.d_relativeLoc == "right") {
        lower = dg.d_p1.d_x;
        upper = dg.d_p1.d_x + thickness;
        if (util::compare::definitelyGreaterThan(check,
                                                 lower + damping_geom_tol))
          coeff *= std::pow((check - lower) / thickness, exponent);
        else
          coeff *= 0.;

        // std::cout << "Right layer, coeff = " << coeff << ", y = " << check <<
        // ", lower = " << lower << ", upper = " << upper << std::endl;
      } else {
        lower = dg.d_p1.d_x;
        upper = dg.d_p1.d_x + thickness;
        if (util::compare::definitelyLessThan(check, upper - damping_geom_tol))
          coeff *= std::pow((upper - check) / thickness, exponent);
        else
          coeff *= 0.;

        lower = dg.d_p2.d_x - thickness;
        upper = dg.d_p2.d_x;
        if (util::compare::definitelyGreaterThan(check,
                                                 lower + damping_geom_tol))
          coeff *= std::pow((check - lower) / thickness, exponent);
        else
          coeff *= 0.;
      }
    }

    if (d_dim > 1 && dg.d_checkY) {
      thickness = dg.d_layerThicknessY;
      check = x.d_y;

      // check if x coordinate is within the region
      if (dg.d_relativeLoc == "bottom") {
        // special care when point is exactly on the upper line
        lower = dg.d_p2.d_y - thickness;
        upper = dg.d_p2.d_y;
        if (util::compare::definitelyLessThan(check, upper - damping_geom_tol))
          coeff *= std::pow((upper - check) / thickness, exponent);
        else
          coeff *= 0.;

        // std::cout << "Bottom layer, coeff = " << coeff << ", y = " << check
        // << ", lower = " << lower << ", upper = " << upper << std::endl;
      } else if (dg.d_relativeLoc == "top") {
        lower = dg.d_p1.d_y;
        upper = dg.d_p1.d_y + thickness;
        if (util::compare::definitelyGreaterThan(check,
                                                 lower + damping_geom_tol))
          coeff *= std::pow((check - lower) / thickness, exponent);
        else
          coeff *= 0.;

        // std::cout << "Top layer, coeff = " << coeff << ", y = " << check <<
        // ", lower = " << lower << ", upper = " << upper << std::endl;
      } else {
        lower = dg.d_p1.d_y;
        upper = dg.d_p1.d_y + thickness;
        if (util::compare::definitelyLessThan(check, upper - damping_geom_tol))
          coeff *= std::pow((upper - check) / thickness, exponent);
        else
          coeff *= 0.;

        lower = dg.d_p2.d_y - thickness;
        upper = dg.d_p2.d_y;
        if (util::compare::definitelyGreaterThan(check,
                                                 lower + damping_geom_tol))
          coeff *= std::pow((check - lower) / thickness, exponent);
        else
          coeff *= 0.;
      }
    }

    d_coefficients[i] = coeff;
  }
}

double geometry::DampingGeom::getCoefficient(const size_t &i) {
  return d_coefficients[i];
}
double geometry::DampingGeom::getCoefficient(const size_t &i) const {
  return d_coefficients[i];
}

const std::vector<double> *geometry::DampingGeom::getCoefficientDataP() {
  return &d_coefficients;
}
const std::vector<double> *geometry::DampingGeom::getCoefficientDataP() const {
  return &d_coefficients;
}

bool geometry::DampingGeom::isDampingActive() {
  return d_absorbingDeck_p->d_dampingActive;
}
bool geometry::DampingGeom::isDampingActive() const {
  return d_absorbingDeck_p->d_dampingActive;
}

bool geometry::DampingGeom::isViscousDamping() {
  return d_absorbingDeck_p->d_isViscousDamping;
}
bool geometry::DampingGeom::isViscousDamping() const {
  return d_absorbingDeck_p->d_isViscousDamping;
}
