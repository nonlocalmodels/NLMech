// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "quadrature.h"
#include "../inp/decks/quadratureDeck.h"
#include <cmath>
#include <util/feElementDefs.h>

int fe::Quadrature::d_debug = 0;

fe::Quadrature::Quadrature(inp::QuadratureDeck *deck, size_t element_type)
    : d_quadOrder(0), d_quadOrderM(0), d_numQuadPts(0), d_numQuadPtsM(0),
      d_eType(element_type) {

  d_quadOrder = deck->d_quadOrder;
  d_quadOrderM = deck->d_quadOrderM;

  if (d_quadOrder == 0)
    d_quadOrder = 1;
  if (d_quadOrderM == 0)
    d_quadOrderM = d_quadOrder;

  // populate quad points
  if (d_eType == util::vtk_type_triangle) {
    initTri(d_quads, d_quadOrder);
    initTri(d_quadsM, d_quadOrderM);
  } else if (d_eType == util::vtk_type_quad) {
    initQuadrangle(d_quads, d_quadOrder);
    initQuadrangle(d_quadsM, d_quadOrderM);
  }
}

void fe::Quadrature::turnDebugOff() { d_debug = 0; }
void fe::Quadrature::turnDebugOn() { d_debug = 1; }

std::vector<fe::QuadData>
fe::Quadrature::getQuadPoints(const std::vector<util::Point3> &nodes,
                              const int &order) {

  if (d_eType == util::vtk_type_triangle)
    getQuadPointsTriangle(nodes, order);
  else if (d_eType == util::vtk_type_quad)
    getQuadPointsQuadrangle(nodes, order);
}

std::vector<fe::QuadData>
fe::Quadrature::getQuadPointsTriangle(const std::vector<util::Point3> &nodes,
                                      const int &order) {
  //
  // Map vertices of given triangle to reference triangle. Map first vertex
  // in nodes to the first vertex (0,0) of triangle, similarly, map second
  // and third to the second and third vertices (1,0) and (0,1) of reference
  // triangle.
  //
  // Caller needs to ensure that order does not go higher than 5.
  std::vector<fe::QuadData> qds = d_quads[order];

  // Since mapping will leave values of shape function unchanged, we only
  // need to modify the positions of quad points in qds and map it to the
  // given triangle, and we also need to modify the weights.
  for (auto qd : qds) {

    qd.d_w *= mapRefTriToTri(qd.d_p, qd.d_shapes, qd.d_derShapes, nodes);
  }

  return qds;
}

std::vector<fe::QuadData>
fe::Quadrature::getQuadPointsQuadrangle(const std::vector<util::Point3> &nodes,
                                        const int &order) {
  return std::vector<fe::QuadData>();
}

std::vector<double> fe::Quadrature::getTriShapes(const util::Point3 &p) {

  // N1 = 1 - xi - eta, N2 = xi, N3 = eta
  return std::vector<double>{1. - p.d_x - p.d_y, p.d_x, p.d_y};
}

std::vector<std::vector<double>> fe::Quadrature::getTriDerShapes() {

  // d N1/d xi = -1, d N1/d eta = -1, d N2/ d xi = 1, d N2/d eta = 0,
  // d N3/ d xi = 0, d N3/d eta = 1
  std::vector<std::vector<double>> r;
  r.push_back(std::vector<double>{-1., -1.});
  r.push_back(std::vector<double>{1., 0.});
  r.push_back(std::vector<double>{0., 1.});

  return r;
}

double fe::Quadrature::mapRefTriToTri(util::Point3 &p, const std::vector<double> &shapes,
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

std::vector<double> fe::Quadrature::getQuadShapes(const util::Point3 &p) {

  // N1 = (1 - xi)(1 - eta)/4
  // N2 = (1 + xi)(1 - eta)/4
  // N3 = (1 + xi)(1 + eta)/4
  // N4 = (1 - xi)(1 + eta)/4
  return std::vector<double>{
      0.25 * (1. - p.d_x) * (1. - p.d_y), 0.25 * (1. + p.d_x) * (1. - p.d_y),
      0.25 * (1. + p.d_x) * (1. + p.d_y), 0.25 * (1. - p.d_x) * (1. + p.d_y)};
}

std::vector<std::vector<double>> fe::Quadrature::getQuadDerShapes(const util::Point3 &p) {

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

void fe::Quadrature::initTri(std::vector<fe::QuadData> &quads,
                             const int &order) {

  //
  // compute quad data for reference triangle with vertex at
  // (0,0), (1,0), (0,1)
  //

  if (!quads.empty())
    return;

  // no point in zeroth order
  if (order == 0) {
    quads.resize(0);
  }

  //
  // first order quad points for triangle
  //
  if (order == 1) {
    quads.clear();
    fe::QuadData qd;
    qd.d_w = 0.5;
    qd.d_p = util::Point3(1. / 3., 1. / 3., 0.);
    // N1 = 1 - xi - eta, N2 = xi, N3 = eta
    qd.d_shapes = getTriShapes(qd.d_p);
    // d N1/d xi = -1, d N1/d eta = -1, d N2/ d xi = 1, d N2/d eta = 0,
    // d N3/ d xi = 0, d N3/d eta = 1
    qd.d_derShapes = getTriDerShapes();
    quads.push_back(qd);
  }

  //
  // second order quad points for triangle
  //
  if (order == 2) {
    quads.clear();
    fe::QuadData qd;
    // point 1
    qd.d_w = 1. / 3.;
    qd.d_p = util::Point3(1. / 6., 1. / 6., 0.);
    qd.d_shapes = getTriShapes(qd.d_p);
    qd.d_derShapes = getTriDerShapes();
    quads.push_back(qd);
    // point 2
    qd.d_w = 1. / 3.;
    qd.d_p = util::Point3(2. / 3., 1. / 6., 0.);
    qd.d_shapes = getTriShapes(qd.d_p);
    qd.d_derShapes = getTriDerShapes();
    quads.push_back(qd);
    // point 3
    qd.d_w = 1. / 3.;
    qd.d_p = util::Point3(1. / 6., 2. / 3., 0.);
    qd.d_shapes = getTriShapes(qd.d_p);
    qd.d_derShapes = getTriDerShapes();
    quads.push_back(qd);
  }

  //
  // third order quad points for triangle
  //
  if (order == 3) {
    quads.clear();
    fe::QuadData qd;
    // point 1
    qd.d_w = -27. / 48.;
    qd.d_p = util::Point3(1. / 3., 1. / 3., 0.);
    qd.d_shapes = getTriShapes(qd.d_p);
    qd.d_derShapes = getTriDerShapes();
    quads.push_back(qd);
    // point 2
    qd.d_w = 25. / 48.;
    qd.d_p = util::Point3(1. / 5., 3. / 5., 0.);
    qd.d_shapes = getTriShapes(qd.d_p);
    qd.d_derShapes = getTriDerShapes();
    quads.push_back(qd);
    // point 3
    qd.d_w = 25. / 48.;
    qd.d_p = util::Point3(1. / 5., 1. / 5., 0.);
    qd.d_shapes = getTriShapes(qd.d_p);
    qd.d_derShapes = getTriDerShapes();
    quads.push_back(qd);
    // point 4
    qd.d_w = 25. / 48.;
    qd.d_p = util::Point3(3. / 5., 1. / 5., 0.);
    qd.d_shapes = getTriShapes(qd.d_p);
    qd.d_derShapes = getTriDerShapes();
    quads.push_back(qd);
  }

  //
  // fourth order quad points for triangle
  //
  if (order == 4) {
    quads.clear();
    fe::QuadData qd;
    // point 1
    qd.d_w = 0.22338158967801;
    qd.d_p = util::Point3(0.44594849091597, 0.44594849091597, 0.);
    qd.d_shapes = getTriShapes(qd.d_p);
    qd.d_derShapes = getTriDerShapes();
    quads.push_back(qd);
    // point 2
    qd.d_w = 0.22338158967801;
    qd.d_p = util::Point3(0.44594849091597, 0.10810301816807, 0.);
    qd.d_shapes = getTriShapes(qd.d_p);
    qd.d_derShapes = getTriDerShapes();
    quads.push_back(qd);
    // point 3
    qd.d_w = 0.22338158967801;
    qd.d_p = util::Point3(0.10810301816807, 0.44594849091597, 0.);
    qd.d_shapes = getTriShapes(qd.d_p);
    qd.d_derShapes = getTriDerShapes();
    quads.push_back(qd);
    // point 4
    qd.d_w = 0.10995174365532;
    qd.d_p = util::Point3(0.09157621350977, 0.09157621350977, 0.);
    qd.d_shapes = getTriShapes(qd.d_p);
    qd.d_derShapes = getTriDerShapes();
    quads.push_back(qd);
    // point 5
    qd.d_w = 0.10995174365532;
    qd.d_p = util::Point3(0.09157621350977, 0.81684757298046, 0.);
    qd.d_shapes = getTriShapes(qd.d_p);
    qd.d_derShapes = getTriDerShapes();
    quads.push_back(qd);
    // point 6
    qd.d_w = 0.10995174365532;
    qd.d_p = util::Point3(0.81684757298046, 0.09157621350977, 0.);
    qd.d_shapes = getTriShapes(qd.d_p);
    qd.d_derShapes = getTriDerShapes();
    quads.push_back(qd);
  }

  //
  // fifth order quad points for triangle
  //
  if (order == 5) {
    quads.clear();
    fe::QuadData qd;
    // point 1
    qd.d_w = 0.22500000000000;
    qd.d_p = util::Point3(0.33333333333333, 0.33333333333333, 0.);
    qd.d_shapes = getTriShapes(qd.d_p);
    qd.d_derShapes = getTriDerShapes();
    quads.push_back(qd);
    // point 2
    qd.d_w = 0.13239415278851;
    qd.d_p = util::Point3(0.47014206410511, 0.47014206410511, 0.);
    qd.d_shapes = getTriShapes(qd.d_p);
    qd.d_derShapes = getTriDerShapes();
    quads.push_back(qd);
    // point 3
    qd.d_w = 0.13239415278851;
    qd.d_p = util::Point3(0.47014206410511, 0.05971587178977, 0.);
    qd.d_shapes = getTriShapes(qd.d_p);
    qd.d_derShapes = getTriDerShapes();
    quads.push_back(qd);
    // point 4
    qd.d_w = 0.13239415278851;
    qd.d_p = util::Point3(0.05971587178977, 0.47014206410511, 0.);
    qd.d_shapes = getTriShapes(qd.d_p);
    qd.d_derShapes = getTriDerShapes();
    quads.push_back(qd);
    // point 5
    qd.d_w = 0.12593918054483;
    qd.d_p = util::Point3(0.10128650732346, 0.10128650732346, 0.);
    qd.d_shapes = getTriShapes(qd.d_p);
    qd.d_derShapes = getTriDerShapes();
    quads.push_back(qd);
    // point 6
    qd.d_w = 0.12593918054483;
    qd.d_p = util::Point3(0.10128650732346, 0.79742698535309, 0.);
    qd.d_shapes = getTriShapes(qd.d_p);
    qd.d_derShapes = getTriDerShapes();
    quads.push_back(qd);
    // point 7
    qd.d_w = 0.12593918054483;
    qd.d_p = util::Point3(0.79742698535309, 0.10128650732346, 0.);
    qd.d_shapes = getTriShapes(qd.d_p);
    qd.d_derShapes = getTriDerShapes();
    quads.push_back(qd);
  }
}

void fe::Quadrature::initQuadrangle(std::vector<fe::QuadData> &quads,
                                    const int &order) {

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

  if (!quads.empty())
    return;

  // no point in zeroth order
  if (order == 0)
    quads.resize(0);

  //
  // first order quad points for quadrangle
  //
  if (order == 1) {
    quads.clear();
    // 1-d points are: {0} and weights are: {2}
    int npts = 1;
    std::vector<double> x = std::vector<double>(1, 0.);
    std::vector<double> w = std::vector<double>(1, 2.);
    for (size_t i = 0; i < npts; i++)
      for (size_t j = 0; j < npts; j++) {

        fe::QuadData qd;
        qd.d_w = w[i] * w[j];
        qd.d_p = util::Point3(x[i], x[j], 0.);
        qd.d_shapes = getQuadShapes(qd.d_p);
        qd.d_derShapes = getQuadDerShapes(qd.d_p);
        quads.push_back(qd);
      }
  }

  //
  // second order quad points for triangle
  //
  if (order == 2) {
    quads.clear();
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
        qd.d_shapes = getQuadShapes(qd.d_p);
        qd.d_derShapes = getQuadDerShapes(qd.d_p);
        quads.push_back(qd);
      }
  }

  //
  // third order quad points for triangle
  //
  if (order == 3) {
    quads.clear();
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
        qd.d_shapes = getQuadShapes(qd.d_p);
        qd.d_derShapes = getQuadDerShapes(qd.d_p);
        quads.push_back(qd);
      }
  }

  //
  // fourth order quad points for triangle
  //
  if (order == 4) {
    quads.clear();
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
        qd.d_shapes = getQuadShapes(qd.d_p);
        qd.d_derShapes = getQuadDerShapes(qd.d_p);
        quads.push_back(qd);
      }
  }

  //
  // fifth order quad points for triangle
  //
  if (order == 5) {
    quads.clear();
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
        qd.d_shapes = getQuadShapes(qd.d_p);
        qd.d_derShapes = getQuadDerShapes(qd.d_p);
        quads.push_back(qd);
      }
  }
}
