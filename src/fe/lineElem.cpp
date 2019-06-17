// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "lineElem.h"
#include "util/feElementDefs.h"     // global definition of elements

fe::LineElem::LineElem(size_t order)
    : fe::BaseElem(order, util::vtk_type_line) {

  // compute quad data
  this->init();
}

std::vector<fe::QuadData>
fe::LineElem::getQuadPoints(const std::vector<util::Point3> &nodes) {
  //
  // Map vertices of given line to reference line.
  //
  // Caller needs to ensure that order does not go higher than 5.
  std::vector<fe::QuadData> qds = d_quads;

  // Since mapping will leave values of shape function unchanged, we only
  // need to modify the positions of quad points in qds and map it to the
  // given line element, and we also need to modify the weights.
  for (auto &i : qds) {

    fe::QuadData *qd = &i;
    qd->d_w = qd->d_w *
              mapRefElemToElem(qd->d_p, qd->d_shapes, qd->d_derShapes, nodes);
  }

  return qds;
}

std::vector<double> fe::LineElem::getShapes(const util::Point3 &p) {

  // N1 = (1 - xi)/2
  // N2 = (1 + xi)/2
  return std::vector<double>{0.5 * (1. - p.d_x), 0.5 * (1. + p.d_x)};
}

std::vector<std::vector<double>>
fe::LineElem::getDerShapes(const util::Point3 &p) {

  // N1 = (1 - xi)/2
  // --> d N1/d xi = -1/2
  //
  // N2 = (1 + xi)/2
  // --> d N2/d xi = 1/2
  std::vector<std::vector<double>> r;
  r.push_back(std::vector<double>{-0.5});
  r.push_back(std::vector<double>{0.5});
  return r;
}

double fe::LineElem::mapRefElemToElem(
    util::Point3 &p, const std::vector<double> &shapes,
    const std::vector<std::vector<double>> &der_shapes,
    const std::vector<util::Point3> &nodes) {

  //
  // see function descriptor for details
  //
  p.d_x = shapes[0] * nodes[0].d_x + shapes[1] * nodes[1].d_x;

  return der_shapes[0][0] * nodes[0].d_x + der_shapes[1][0] * nodes[1].d_x;
}

void fe::LineElem::init() {

  //
  // compute quad data for reference line element with vertex at p1 = -1 and
  // p2 = 1
  //
  // Shape functions are
  // N1 = (1 - xi)/2
  // N2 = (1 + xi)/2

  if (!d_quads.empty())
    return;

  // no point in zeroth order
  if (d_quadOrder == 0)
    d_quads.resize(0);

  //
  // first order quad points
  //
  if (d_quadOrder == 1) {
    d_quads.clear();
    // 1-d points are: {0} and weights are: {2}
    fe::QuadData qd;
    qd.d_w = 2.;
    qd.d_p = util::Point3();
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);
  }

  //
  // second order quad points
  //
  if (d_quadOrder == 2) {
    d_quads.clear();
    // 1-d points are: {-1/sqrt{3], 1/sqrt{3}} and weights are: {1,1}
    fe::QuadData qd;
    qd.d_w = 1.;
    qd.d_p = util::Point3(-1. / std::sqrt(3.), 0., 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);

    qd.d_w = 1.;
    qd.d_p = util::Point3(1. / std::sqrt(3.), 0., 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);
  }

  //
  // third order quad points for triangle
  //
  if (d_quadOrder == 3) {
    d_quads.clear();
    // 1-d points are: {-sqrt{3}/sqrt{5}, 0, sqrt{3}/sqrt{5}}
    // weights are: {5/9, 8/9, 5/9}
    fe::QuadData qd;
    qd.d_w = 5. / 9.;
    qd.d_p = util::Point3(-std::sqrt(3.) / std::sqrt(5.), 0., 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);

    qd.d_w = 8. / 9.;
    qd.d_p = util::Point3();
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);

    qd.d_w = 5. / 9.;
    qd.d_p = util::Point3(std::sqrt(3.) / std::sqrt(5.), 0., 0.);
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
    qd.d_w = 0.6521451548625461;
    qd.d_p = util::Point3(-0.3399810435848563, 0., 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);

    qd.d_w = 0.6521451548625461;
    qd.d_p = util::Point3(0.3399810435848563, 0., 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);

    qd.d_w = 0.3478548451374538;
    qd.d_p = util::Point3(-0.8611363115940526, 0., 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);

    qd.d_w = 0.3478548451374538;
    qd.d_p = util::Point3(0.8611363115940526, 0., 0.);
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
    qd.d_w = 0.5688888888888889;
    qd.d_p = util::Point3();
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);

    qd.d_w = 0.4786286704993665;
    qd.d_p = util::Point3(-0.5384693101056831, 0., 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);

    qd.d_w = 0.4786286704993665;
    qd.d_p = util::Point3(0.5384693101056831, 0., 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);

    qd.d_w = 0.2369268850561891;
    qd.d_p = util::Point3(-0.9061798459386640, 0., 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);

    qd.d_w = 0.2369268850561891;
    qd.d_p = util::Point3(0.9061798459386640, 0., 0.);
    qd.d_shapes = getShapes(qd.d_p);
    qd.d_derShapes = getDerShapes(qd.d_p);
    d_quads.push_back(qd);
  }
}