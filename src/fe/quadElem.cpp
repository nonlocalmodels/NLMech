// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "quadElem.h"
#include "../util/feElementDefs.h" // global definition of elements

fe::QuadElem::QuadElem(size_t order)
    : fe::BaseElem(order, util::vtk_type_quad) {

  // compute quad data
  this->init();
}

std::vector<fe::QuadData>
fe::QuadElem::getQuadPoints(const std::vector<util::Point3> &nodes) {
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

std::vector<double> fe::QuadElem::getShapes(const util::Point3 &p) {

  // N1 = (1 - xi)(1 - eta)/4
  // N2 = (1 + xi)(1 - eta)/4
  // N3 = (1 + xi)(1 + eta)/4
  // N4 = (1 - xi)(1 + eta)/4
  return std::vector<double>{
      0.25 * (1. - p.d_x) * (1. - p.d_y), 0.25 * (1. + p.d_x) * (1. - p.d_y),
      0.25 * (1. + p.d_x) * (1. + p.d_y), 0.25 * (1. - p.d_x) * (1. + p.d_y)};
}

std::vector<std::vector<double>>
fe::QuadElem::getDerShapes(const util::Point3 &p) {

  // N1 = (1 - xi)(1 - eta)/4
  // --> d N1/d xi = -(1 - eta)/4, d N1/d eta = -(1 - xi)/4
  //
  // N2 = (1 + xi)(1 - eta)/4
  // --> d N2/d xi = (1 - eta)/4, d N2/d eta = -(1 + xi)/4
  //
  // N3 = (1 + xi)(1 + eta)/4
  // --> d N3/d xi = (1 + eta)/4, d N3/d eta = (1 + xi)/4
  //
  // N4 = (1 - xi)(1 + eta)/4
  // --> d N4/d xi = -(1 + eta)/4, d N4/d eta = (1 - xi)/4
  std::vector<std::vector<double>> r;
  r.push_back(std::vector<double>{-0.25 * (1. - p.d_y), -0.25 * (1. - p.d_x)});
  r.push_back(std::vector<double>{0.25 * (1. - p.d_y), -0.25 * (1. + p.d_x)});
  r.push_back(std::vector<double>{0.25 * (1. + p.d_y), 0.25 * (1. + p.d_x)});
  r.push_back(std::vector<double>{-0.25 * (1. + p.d_y), 0.25 * (1. - p.d_x)});

  return r;
}

// TODO
double fe::QuadElem::mapRefElemToElem(
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

void fe::QuadElem::init() {

  //
  // compute quad data for reference quadrangle with vertex at
  // p1 = (-1,-1), p2 = (1,-1), p3 = (1,1), p4 = (-1,1)
  //
  // Shape functions are
  // N1 = (1 - xi)(1 - eta)/4
  // N2 = (1 + xi)(1 - eta)/4
  // N3 = (1 + xi)(1 + eta)/4
  // N4 = (1 - xi)(1 + eta)/4
  //
  //
  //  Let [-1,1] is the 1-d reference element and {x1, x2, x3,.., xN} are N
  //  quad points for 1-d domain and {w1, w2, w3,..., wN} are respective
  //  weights. Then, the Nth order quad points in Quadrangle [-1,1]x[-1,1] is
  //  simply given by N^2 points and
  //
  //  (i,j) point is (xi, xj) and weight is wi \times wj
  //

  if (!d_quads.empty())
    return;

  // no point in zeroth order
  if (d_quadOrder == 0)
    d_quads.resize(0);

  //
  // first order quad points for quadrangle
  //
  if (d_quadOrder == 1) {
    d_quads.clear();
    // 1-d points are: {0} and weights are: {2}
    int npts = 1;
    std::vector<double> x = std::vector<double>(1, 0.);
    std::vector<double> w = std::vector<double>(1, 2.);
    for (size_t i = 0; i < npts; i++)
      for (size_t j = 0; j < npts; j++) {

        fe::QuadData qd;
        qd.d_w = w[i] * w[j];
        qd.d_p = util::Point3(x[i], x[j], 0.);
        qd.d_shapes = getShapes(qd.d_p);
        qd.d_derShapes = getDerShapes(qd.d_p);
        d_quads.push_back(qd);
      }
  }

  //
  // second order quad points for triangle
  //
  if (d_quadOrder == 2) {
    d_quads.clear();
    // 1-d points are: {-1/sqrt{3], 1/sqrt{3}} and weights are: {1,1}
    int npts = 2;
    std::vector<double> x =
        std::vector<double>{-1. / std::sqrt(3.), 1. / std::sqrt(3.)};
    std::vector<double> w = std::vector<double>{1., 1.};
    for (size_t i = 0; i < npts; i++)
      for (size_t j = 0; j < npts; j++) {

        fe::QuadData qd;
        qd.d_w = w[i] * w[j];
        qd.d_p = util::Point3(x[i], x[j], 0.);
        qd.d_shapes = getShapes(qd.d_p);
        qd.d_derShapes = getDerShapes(qd.d_p);
        d_quads.push_back(qd);
      }
  }

  //
  // third order quad points for triangle
  //
  if (d_quadOrder == 3) {
    d_quads.clear();
    // 1-d points are: {-sqrt{3}/sqrt{5}, 0, sqrt{3}/sqrt{5}}
    // weights are: {5/9, 8/9, 5/9}
    int npts = 3;
    std::vector<double> x = std::vector<double>{
        -std::sqrt(3.) / std::sqrt(5.), 0., std::sqrt(3.) / std::sqrt(5.)};
    std::vector<double> w = std::vector<double>{5. / 9., 8. / 9., 5. / 9.};
    for (size_t i = 0; i < npts; i++)
      for (size_t j = 0; j < npts; j++) {

        fe::QuadData qd;
        qd.d_w = w[i] * w[j];
        qd.d_p = util::Point3(x[i], x[j], 0.);
        qd.d_shapes = getShapes(qd.d_p);
        qd.d_derShapes = getDerShapes(qd.d_p);
        d_quads.push_back(qd);
      }
  }

  //
  // fourth order quad points for triangle
  //
  if (d_quadOrder == 4) {
    d_quads.clear();
    int npts = 4;
    std::vector<double> x =
        std::vector<double>{-0.3399810435848563, 0.3399810435848563,
                            -0.8611363115940526, 0.8611363115940526};
    std::vector<double> w =
        std::vector<double>{0.6521451548625461, 0.6521451548625461,
                            0.3478548451374538, 0.3478548451374538};
    for (size_t i = 0; i < npts; i++)
      for (size_t j = 0; j < npts; j++) {

        fe::QuadData qd;
        qd.d_w = w[i] * w[j];
        qd.d_p = util::Point3(x[i], x[j], 0.);
        qd.d_shapes = getShapes(qd.d_p);
        qd.d_derShapes = getDerShapes(qd.d_p);
        d_quads.push_back(qd);
      }
  }

  //
  // fifth order quad points for triangle
  //
  if (d_quadOrder == 5) {
    d_quads.clear();
    int npts = 5;
    std::vector<double> x =
        std::vector<double>{0., -0.5384693101056831, 0.5384693101056831,
                            -0.9061798459386640, 0.9061798459386640};
    std::vector<double> w = std::vector<double>{
        0.5688888888888889, 0.4786286704993665, 0.4786286704993665,
        0.2369268850561891, 0.2369268850561891};
    for (size_t i = 0; i < npts; i++)
      for (size_t j = 0; j < npts; j++) {

        fe::QuadData qd;
        qd.d_w = w[i] * w[j];
        qd.d_p = util::Point3(x[i], x[j], 0.);
        qd.d_shapes = getShapes(qd.d_p);
        qd.d_derShapes = getDerShapes(qd.d_p);
        d_quads.push_back(qd);
      }
  }
}