// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "triElem.h"
#include "../util/feElementDefs.h" // global definition of elements

fe::TriElem::TriElem(const size_t &order)
    : fe::BaseElem(order, util::vtk_type_triangle) {

  // compute quad data
  this->init();
}

std::vector<fe::QuadData>
fe::TriElem::getQuadPoints(const std::vector<util::Point3> &nodes) {
  //
  // Map vertices of given triangle to reference triangle. Map first vertex
  // in nodes to the first vertex (0,0) of triangle, similarly, map second
  // and third to the second and third vertices (1,0) and (0,1) of reference
  // triangle.
  //
  // Caller needs to ensure that order does not go higher than 5.
  std::vector<fe::QuadData> qds = d_quads;

  // Since mapping will leave values of shape function unchanged, we only
  // need to modify the positions of quad points in qds and map it to the
  // given triangle, and we also need to modify the weights.
  for (auto qd : qds) {

    qd.d_w *= mapRefElemToElem(qd.d_p, qd.d_shapes, qd.d_derShapes, nodes);
  }

  return qds;
}

std::vector<double> fe::TriElem::getShapes(const util::Point3 &p) {

  // N1 = 1 - xi - eta, N2 = xi, N3 = eta
  return std::vector<double>{1. - p.d_x - p.d_y, p.d_x, p.d_y};
}

std::vector<std::vector<double>>
fe::TriElem::getDerShapes(const util::Point3 &p) {

  // d N1/d xi = -1, d N1/d eta = -1, d N2/ d xi = 1, d N2/d eta = 0,
  // d N3/ d xi = 0, d N3/d eta = 1
  std::vector<std::vector<double>> r;
  r.push_back(std::vector<double>{-1., -1.});
  r.push_back(std::vector<double>{1., 0.});
  r.push_back(std::vector<double>{0., 1.});

  return r;
}

double fe::TriElem::mapRefElemToElem(
    util::Point3 &p, const std::vector<double> &shapes,
    const std::vector<std::vector<double>> &der_shapes,
    const std::vector<util::Point3> &nodes) {

  //
  // see function descriptor for details
  //
  p.d_x = shapes[0] * nodes[0].d_x + shapes[1] * nodes[1].d_x +
          shapes[2] * nodes[2].d_x;
  p.d_y = shapes[0] * nodes[0].d_y + shapes[1] * nodes[1].d_y +
          shapes[2] * nodes[2].d_y;

  return (nodes[1].d_x - nodes[0].d_x) * (nodes[2].d_y - nodes[0].d_y) -
         (nodes[2].d_x - nodes[1].d_x) * (nodes[1].d_y - nodes[0].d_y);
}

void fe::TriElem::init() {

  //
  // compute quad data for reference triangle with vertex at
  // (0,0), (1,0), (0,1)
  //

  if (!d_quads.empty())
    return;

  // no point in zeroth order
  if (d_quadOrder == 0) {
    d_quads.resize(0);
  }

  //
  // first order quad points for triangle
  //
  if (d_quadOrder == 1) {
    d_quads.clear();
    fe::QuadData qd;
    qd.d_w = 0.5;
    qd.d_p = util::Point3(1. / 3., 1. / 3., 0.);
    // N1 = 1 - xi - eta, N2 = xi, N3 = eta
    qd.d_shapes = getShapes(qd.d_p);
    // d N1/d xi = -1, d N1/d eta = -1, d N2/ d xi = 1, d N2/d eta = 0,
    // d N3/ d xi = 0, d N3/d eta = 1
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);
  }

  //
  // second order quad points for triangle
  //
  if (d_quadOrder == 2) {
    d_quads.clear();
    fe::QuadData qd;
    // point 1
    qd.d_w = 1. / 3.;
    qd.d_p = util::Point3(1. / 6., 1. / 6., 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);
    // point 2
    qd.d_w = 1. / 3.;
    qd.d_p = util::Point3(2. / 3., 1. / 6., 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);
    // point 3
    qd.d_w = 1. / 3.;
    qd.d_p = util::Point3(1. / 6., 2. / 3., 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);
  }

  //
  // third order quad points for triangle
  //
  if (d_quadOrder == 3) {
    d_quads.clear();
    fe::QuadData qd;
    // point 1
    qd.d_w = -27. / 48.;
    qd.d_p = util::Point3(1. / 3., 1. / 3., 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);
    // point 2
    qd.d_w = 25. / 48.;
    qd.d_p = util::Point3(1. / 5., 3. / 5., 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);
    // point 3
    qd.d_w = 25. / 48.;
    qd.d_p = util::Point3(1. / 5., 1. / 5., 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);
    // point 4
    qd.d_w = 25. / 48.;
    qd.d_p = util::Point3(3. / 5., 1. / 5., 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);
  }

  //
  // fourth order quad points for triangle
  //
  if (d_quadOrder == 4) {
    d_quads.clear();
    fe::QuadData qd;
    // point 1
    qd.d_w = 0.22338158967801;
    qd.d_p = util::Point3(0.44594849091597, 0.44594849091597, 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);
    // point 2
    qd.d_w = 0.22338158967801;
    qd.d_p = util::Point3(0.44594849091597, 0.10810301816807, 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);
    // point 3
    qd.d_w = 0.22338158967801;
    qd.d_p = util::Point3(0.10810301816807, 0.44594849091597, 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);
    // point 4
    qd.d_w = 0.10995174365532;
    qd.d_p = util::Point3(0.09157621350977, 0.09157621350977, 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);
    // point 5
    qd.d_w = 0.10995174365532;
    qd.d_p = util::Point3(0.09157621350977, 0.81684757298046, 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);
    // point 6
    qd.d_w = 0.10995174365532;
    qd.d_p = util::Point3(0.81684757298046, 0.09157621350977, 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);
  }

  //
  // fifth order quad points for triangle
  //
  if (d_quadOrder == 5) {
    d_quads.clear();
    fe::QuadData qd;
    // point 1
    qd.d_w = 0.22500000000000;
    qd.d_p = util::Point3(0.33333333333333, 0.33333333333333, 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);
    // point 2
    qd.d_w = 0.13239415278851;
    qd.d_p = util::Point3(0.47014206410511, 0.47014206410511, 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);
    // point 3
    qd.d_w = 0.13239415278851;
    qd.d_p = util::Point3(0.47014206410511, 0.05971587178977, 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);
    // point 4
    qd.d_w = 0.13239415278851;
    qd.d_p = util::Point3(0.05971587178977, 0.47014206410511, 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);
    // point 5
    qd.d_w = 0.12593918054483;
    qd.d_p = util::Point3(0.10128650732346, 0.10128650732346, 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);
    // point 6
    qd.d_w = 0.12593918054483;
    qd.d_p = util::Point3(0.10128650732346, 0.79742698535309, 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);
    // point 7
    qd.d_w = 0.12593918054483;
    qd.d_p = util::Point3(0.79742698535309, 0.10128650732346, 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);
  }
}