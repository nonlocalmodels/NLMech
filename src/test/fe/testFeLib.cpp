// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "testFeLib.h"

#include "fe/massMatrix.h"
#include "fe/mesh.h"
#include "fe/quadElem.h"
#include "fe/quadrature.h"
#include "fe/triElem.h"
#include "util/point.h"

#include <fstream>
#include <string>

static const double tol = 1.0E-12;

double getNChooseR(size_t n, size_t r) {

  if (r == 0)
    return 1.;

  double a = 1.;
  for (size_t i = 1; i <= r; i++)
    a *= double(n - i + 1) / double(i);

  return a;
}

double getExactIntegrationRefTri(size_t alpha, size_t beta) {

  // compute exact integration of s^\alpha t^\beta
  double I = 0.;
  for (size_t k = 0; k <= beta + 1; k++) {
    if (k % 2 == 0)
      I += getNChooseR(beta + 1, k) / double((alpha + 1 + k) * (beta + 1));
    else
      I -= getNChooseR(beta + 1, k) / double((alpha + 1 + k) * (beta + 1));
  }

  return I;
}

void checkRefIntegration(const size_t &n, const size_t &i, const size_t &j,
                            const std::vector<fe::QuadData> &qds,
                            double &I_exact) {

  double I_approx = 0.;
  for (auto qd : qds)
    I_approx += qd.d_w * std::pow(qd.d_p.d_x, i) * std::pow(qd.d_p.d_y, j);

  if (std::abs(I_exact - I_approx) > tol)
    std::cout << "Error in order = " << n << ". Exact integration = " << I_exact
              << " and approximate integration = " << I_approx
              << " of polynomial of order (i = " << i << " + j = " << j
              << ") = " << i + j << " is not matching with approximate "
              << "integration = " << I_approx << " using quadrature points.\n";
}

void testTriRef(size_t n) {

  //
  // Test1: We test accuracy of integrals of polynomials over reference
  // triangle. Reference triangle {(0,0), (1,0), (0,1)}.
  //
  // Test2: We consider simple mesh in meshFeTest.txt over square domain
  // [0,1]^2 and test the accuracy of polynomials over square domain.
  //

  // get Quadrature
  fe::Quadrature<fe::TriElem> quad(n);

  //
  // Test 1
  //
  {
    // T1 (reference triangle)
    // get quad points at reference triangle
    std::vector<util::Point3> nodes = {util::Point3(), util::Point3(1., 0., 0.),
                                       util::Point3(0., 1., 0.)};
    std::vector<fe::QuadData> qds = quad.getQuadPoints(nodes);

    std::string filename = "tri_quads_" + std::to_string(n) + ".txt";
    std::ofstream myfile(filename.c_str());
    myfile.precision(10);
    size_t counter = 0;
    double sum = 0.;
    for (auto qd : qds) {
      myfile << ++counter << " " << qd.d_p.d_x << " " << qd.d_p.d_y << " "
             << qd.d_w << " " << qd.d_shapes[0] << " " << qd.d_shapes[1] << " "
             << qd.d_shapes[2] << " " << qd.d_derShapes[0][0] << " "
             << qd.d_derShapes[0][1] << " " << qd.d_derShapes[1][0] << " "
             << qd.d_derShapes[0][1] << " " << qd.d_derShapes[2][0] << " "
             << qd.d_derShapes[2][1] << "\n";
      sum += qd.d_w;
    }

    if (std::abs(sum - 0.5) > tol)
      std::cout << "Error in order = " << n
                << ". Sum of quad weights is not "
                   "equal to area of reference "
                   "triangle.\n";

    //
    // test the exactness of integration for polynomial
    //
    for (size_t i = 0; i <= n; i++)
      for (size_t j = 0; j <= n; j++) {

        if (i + j > n)
          continue;

        //
        // when {(0,0), (1,0), (0,1)}
        //
        nodes = {util::Point3(), util::Point3(1., 0., 0.),
                 util::Point3(0., 1., 0.)};
        qds = quad.getQuadPoints(nodes);
        // test integration of polynomial f(s,t) = s^i t^j
        // get the exact integration
        double I_exact = getExactIntegrationRefTri(i, j);
        checkRefIntegration(n, i, j, qds, I_exact);

        //
        // when vertices are {(1,0), (0,1), (0,0)}
        //
        nodes = {util::Point3(1., 0., 0.), util::Point3(0., 1., 0.),
                 util::Point3()};
        qds = quad.getQuadPoints(nodes);
        //
        // After changing the order of vertices, we have got a new
        // triangle which is in coordinate system (x,y) and we are
        // integrating function f(x,y) = x^i y^j
        //
        // The quad data we have got is such that quad point is in (x,y)
        // coordinate, weight is such that determinant of the Jacobian is
        // included in the weight.
        //
        // Thus the following method for I_approx is correct.
        checkRefIntegration(n, i, j, qds, I_exact);

        //
        // when vertices are {(0,1), (0,0), (1,0)}
        //
        nodes = {util::Point3(0., 1., 0.), util::Point3(),
                 util::Point3(1., 0., 0.)};
        qds = quad.getQuadPoints(nodes);
        checkRefIntegration(n, i, j, qds, I_exact);
      }
  } // Test 1

  //
  // Test 2
  //


}
