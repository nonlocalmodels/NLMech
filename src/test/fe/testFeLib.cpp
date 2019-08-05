////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "testFeLib.h"

#include "../../external/csv.h"
#include "fe/massMatrix.h"
#include "fe/mesh.h"
#include "fe/quadElem.h"
#include "fe/triElem.h"
#include "util/point.h"
#include <util/feElementDefs.h>

#include <algorithm>
#include <fstream>
#include <string>

//
// static methods
//
static const double tol = 1.0E-12;

static void readNodes(const std::string &filename,
                      std::vector<util::Point3> &nodes) {

  // csv reader
  io::CSVReader<3> in(filename);
  double x, y, z;
  while (in.read_row(x, y, z))
    nodes.emplace_back(x, y, z);
}

static size_t readElements(const std::string &filename, const size_t &elem_type,
                           std::vector<size_t> &elements) {

  if (elem_type == util::vtk_type_triangle) {
    io::CSVReader<3> in(filename);
    std::vector<size_t> ids(3, 0);
    while (in.read_row(ids[0], ids[1], ids[2])) {
      for (auto id : ids)
        elements.emplace_back(id);
    }

    size_t num_vertex = util::vtk_map_element_to_num_nodes[elem_type];
    return elements.size() / num_vertex;
  } else if (elem_type == util::vtk_type_quad) {
    io::CSVReader<4> in(filename);
    std::vector<size_t> ids(4, 0);
    while (in.read_row(ids[0], ids[1], ids[2], ids[3])) {
      for (auto id : ids)
        elements.emplace_back(id);
    }

    size_t num_vertex = util::vtk_map_element_to_num_nodes[elem_type];
    return elements.size() / num_vertex;
  }
}

static bool checkRefIntegration(const size_t &n, const size_t &i,
                                const size_t &j,
                                const std::vector<fe::QuadData> &qds,
                                double &I_exact) {

  double I_approx = 0.;
  for (auto qd : qds)
    I_approx += qd.d_w * std::pow(qd.d_p.d_x, i) * std::pow(qd.d_p.d_y, j);

  if (std::abs(I_exact - I_approx) > tol) {
    std::cout << "Error in order = " << n << ". Exact integration = " << I_exact
              << " and approximate integration = " << I_approx
              << " of polynomial of order (i = " << i << " + j = " << j
              << ") = " << i + j << " over reference element "
              << "is not matching using quadrature points.\n";

    return false;
  }

  return true;
}

//
// Interface methods
//

double test::getNChooseR(size_t n, size_t r) {

  if (r == 0)
    return 1.;

  double a = 1.;
  for (size_t i = 1; i <= r; i++)
    a *= double(n - i + 1) / double(i);

  return a;
}

double test::getExactIntegrationRefTri(size_t alpha, size_t beta) {

  // compute exact integration of s^\alpha t^\beta
  double I = 0.;
  for (size_t k = 0; k <= beta + 1; k++) {
    if (k % 2 == 0)
      I += test::getNChooseR(beta + 1, k) / double((alpha + 1 + k) * (beta + 1));
    else
      I -= test::getNChooseR(beta + 1, k) / double((alpha + 1 + k) * (beta + 1));
  }

  return I;
}

double test::getExactIntegrationRefQuad(size_t alpha, size_t beta) {

  // compute exact integration of s^\alpha t^\beta
  if (alpha % 2 == 0 and beta % 2 == 0)
    return 4. / double((alpha + 1) * (beta + 1));
  else
    return 0.;
}

void test::testLineElem(size_t n) {
  return;
}

void test::testTriElem(size_t n) {

  //
  // Test1: We test accuracy of integrals of polynomials over reference
  // triangle. Reference triangle {(0,0), (1,0), (0,1)}.
  //
  // Test2: We consider simple mesh in meshFeTest.txt over square domain
  // [0,1]^2 and test the accuracy of polynomials over square domain.
  //

  // get Quadrature
  auto quad = fe::TriElem(n);

  //
  // Test 1
  //
  size_t error_test_1 = 0;
  {
    // T1 (reference triangle)
    // get quad points at reference triangle
    std::vector<util::Point3> nodes = {util::Point3(), util::Point3(1., 0., 0.),
                                       util::Point3(0., 1., 0.)};
    std::vector<fe::QuadData> qds = quad.getQuadPoints(nodes);
    double sum = 0.;
    for (auto qd : qds)
      sum += qd.d_w;

    if (std::abs(sum - 0.5) > tol) {
      std::cout << "Error in order = " << n
                << ". Sum of quad weights is not "
                   "equal to area of reference "
                   "triangle.\n";
      error_test_1++;
    }

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
        double I_exact = test::getExactIntegrationRefTri(i, j);
        if (!checkRefIntegration(n, i, j, qds, I_exact))
          error_test_1++;

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
        if (!checkRefIntegration(n, i, j, qds, I_exact))
          error_test_1++;

        //
        // when vertices are {(0,1), (0,0), (1,0)}
        //
        nodes = {util::Point3(0., 1., 0.), util::Point3(),
                 util::Point3(1., 0., 0.)};
        qds = quad.getQuadPoints(nodes);
        if (!checkRefIntegration(n, i, j, qds, I_exact))
          error_test_1++;
      }
  } // Test 1

  //
  // Test 2
  //
  size_t error_test_2 = 0;
  {
    static std::vector<util::Point3> nodes;
    static std::vector<size_t> elements;
    static size_t num_vertex = 3;
    static size_t elem_type = util::vtk_type_triangle;
    static size_t num_elems = 0;
    if (num_elems == 0) {
      readNodes("triMesh_nodes.csv", nodes);
      num_elems = readElements("triMesh_elements.csv", elem_type, elements);
    }

    // loop over polynomials
    for (size_t i = 0; i <= n; i++)
      for (size_t j = 0; j <= n; j++) {

        if (i + j > n)
          continue;

        double I_exact = 1. / (double(i + 1) * double(j + 1));
        double I_approx = 0.;
        // loop over elements and compute I_approx
        for (size_t e = 0; e < num_elems; e++) {
          std::vector<util::Point3> enodes = {
              nodes[elements[num_vertex * e + 0]],
              nodes[elements[num_vertex * e + 1]],
              nodes[elements[num_vertex * e + 2]]};
          std::vector<fe::QuadData> qds = quad.getQuadPoints(enodes);
          for (auto qd : qds)
            I_approx +=
                qd.d_w * std::pow(qd.d_p.d_x, i) * std::pow(qd.d_p.d_y, j);
        }

        if (std::abs(I_exact - I_approx) > tol) {
          std::cout << "Error in order = " << n
                    << ". Exact integration = " << I_exact
                    << " and approximate integration = " << I_approx
                    << " of polynomial of order (i = " << i << " + j = " << j
                    << ") = " << i + j << " over square domain [0,1]x[0,1] "
                    << "is not matching using quadrature points.\n";

          error_test_2++;
        }
      }
  }

  if (n == 1) {
    std::cout << "**********************************\n";
    std::cout << "Triangle Quadrature Test\n";
    std::cout << "**********************************\n";
  }
  std::cout << "Quad order = " << n << ". ";
  if (error_test_1 == 0)
    std::cout << "TEST 1 : PASS. ";
  else
    std::cout << "TEST 1 : FAIL. ";
  if (error_test_2 == 0)
    std::cout << "TEST 2 : PASS. ";
  else
    std::cout << "TEST 2 : FAIL. ";
  std::cout << "\n";
}

void test::testQuadElem(size_t n) {

  //
  // Test1: We test accuracy of integrals of polynomials over reference
  // quadrangle. Reference triangle {(-1,-1), (1,-1), (1,1), (-1,1)}.
  //
  // Test2: We consider simple mesh in meshFeTest.txt over square domain
  // [0,1]^2 and test the accuracy of polynomials over square domain.
  //

  // get Quadrature
  auto quad = fe::QuadElem(n);

  //
  // Test 1
  //
  size_t error_test_1 = 0;
  {
    // T1 (reference quadrangle)
    // get quad points at reference triangle
    std::vector<util::Point3> nodes = {
        util::Point3(-1., -1., 0.), util::Point3(1., -1., 0.),
        util::Point3(1., 1., 0.), util::Point3(-1., 1., 0.)};
    std::vector<fe::QuadData> qds = quad.getQuadPoints(nodes);
    double sum = 0.;
    for (auto qd : qds)
      sum += qd.d_w;

    if (std::abs(sum - 4.0) > tol) {
      std::cout << "Error in order = " << n
                << ". Sum of quad weights is not "
                   "equal to area of reference "
                   "quadrangle.\n";
      error_test_1++;
    }

    //
    // test the exactness of integration for polynomial
    //
    for (size_t i = 0; i <= 2 * n - 1; i++)
      for (size_t j = 0; j <= 2 * n - 1; j++) {

        //
        // when {(-1,-1), (1,-1), (1,1), (-1,1)}
        //
        nodes = {util::Point3(-1., -1., 0.), util::Point3(1., -1., 0.),
                 util::Point3(1., 1., 0.), util::Point3(-1., 1., 0.)};
        qds = quad.getQuadPoints(nodes);
        // test integration of polynomial f(s,t) = s^i t^j
        // get the exact integration
        double I_exact = test::getExactIntegrationRefQuad(i, j);
        if (!checkRefIntegration(n, i, j, qds, I_exact))
          error_test_1++;

        //
        // when {(-1,1), (-1,-1), (1,-1), (1,1)}
        //
        nodes = {util::Point3(-1., 1., 0.), util::Point3(-1., -1., 0.),
                 util::Point3(1., -1., 0.), util::Point3(1., 1., 0.)};
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
        if (!checkRefIntegration(n, i, j, qds, I_exact))
          error_test_1++;

        //
        // when {(1,1), (-1,1), (-1,-1), (1,-1)}
        //
        nodes = {util::Point3(1., 1., 0.), util::Point3(-1., 1., 0.),
                 util::Point3(-1., -1., 0.), util::Point3(1., -1., 0.)};
        qds = quad.getQuadPoints(nodes);
        if (!checkRefIntegration(n, i, j, qds, I_exact))
          error_test_1++;

        //
        // when {(1,-1), (1,1), (-1,1), (-1,-1)}
        //
        nodes = {util::Point3(1., -1., 0.), util::Point3(1., 1., 0.),
                 util::Point3(-1., 1., 0.), util::Point3(-1., -1., 0.)};
        qds = quad.getQuadPoints(nodes);
        if (!checkRefIntegration(n, i, j, qds, I_exact))
          error_test_1++;
      }
  } // Test 1

  //
  // Test 2
  //
  size_t error_test_2 = 0;
  {
    static std::vector<util::Point3> nodes;
    static std::vector<size_t> elements;
    static size_t num_vertex = 4;
    static size_t elem_type = util::vtk_type_quad;
    static size_t num_elems = 0;
    if (num_elems == 0) {
      readNodes("quadMesh_nodes.csv", nodes);
      num_elems = readElements("quadMesh_elements.csv", elem_type, elements);
    }

    // loop over polynomials
    for (size_t i = 0; i <= 2 * n - 1; i++)
      for (size_t j = 0; j <= 2 * n - 1; j++) {

        double I_exact = 1. / (double(i + 1) * double(j + 1));
        double I_approx = 0.;
        // loop over elements and compute I_approx
        for (size_t e = 0; e < num_elems; e++) {
          std::vector<util::Point3> enodes = {
              nodes[elements[num_vertex * e + 0]],
              nodes[elements[num_vertex * e + 1]],
              nodes[elements[num_vertex * e + 2]],
              nodes[elements[num_vertex * e + 3]]};
          std::vector<fe::QuadData> qds = quad.getQuadPoints(enodes);
          for (auto qd : qds)
            I_approx +=
                qd.d_w * std::pow(qd.d_p.d_x, i) * std::pow(qd.d_p.d_y, j);
        }

        if (std::abs(I_exact - I_approx) > tol) {
          std::cout << "Error in order = " << n
                    << ". Exact integration = " << I_exact
                    << " and approximate integration = " << I_approx
                    << " of polynomial of order (i = " << i << " + j = " << j
                    << ") = " << i + j << " over square domain [0,1]x[0,1] "
                    << "is not matching using quadrature points.\n";

          error_test_2++;
        }
      }
  }

  if (n == 1) {
    std::cout << "**********************************\n";
    std::cout << "Quadrangle Quadrature Test\n";
    std::cout << "**********************************\n";
  }
  std::cout << "Quad order = " << n << ". ";
  std::cout << (error_test_1 == 0 ? "TEST 1 : PASS. " : "TEST 1 : FAIL. ");
  std::cout << (error_test_2 == 0 ? "TEST 2 : PASS. \n" : "TEST 2 : FAIL. \n");
}

void test::testTriElemTime(size_t n, size_t N) {

  // get Quadrature
  auto quad = fe::TriElem(n);

  //
  // Test 1
  //
  std::vector<util::Point3> nodes = {util::Point3(2., 2., 0.),
                                     util::Point3(4., 2., 0.),
                                     util::Point3(2., 4., 0.)};
  size_t num_vertex = 3;
  //  std::vector<size_t> elements;
  //  for (size_t i = 0; i < 3 * N; i++)
  //    elements.emplace_back(i % 2);
  std::vector<std::vector<size_t>> elements;
  for (size_t i = 0; i < N; i++)
    elements.emplace_back(std::vector<size_t>{0, 1, 2});

  std::uint64_t t11 = hpx::util::high_resolution_clock::now();
  // method 1: Compute quad points on the fly
  // loop over elements and compute I_approx
  double sum = 0.;
  for (size_t e = 0; e < N; e++) {
    //    std::vector<util::Point3> enodes = {nodes[elements[num_vertex * e +
    //    0]],
    //                                        nodes[elements[num_vertex * e +
    //                                        1]], nodes[elements[num_vertex * e
    //                                        + 2]]};
    std::vector<util::Point3> enodes = {
        nodes[elements[e][0]], nodes[elements[e][1]], nodes[elements[e][2]]};
    std::vector<fe::QuadData> qds = quad.getQuadPoints(enodes);
    for (auto qd : qds)
      sum += qd.d_w * (qd.d_shapes[0] + qd.d_shapes[1] + qd.d_shapes[2]);
  }
  std::uint64_t t12 = hpx::util::high_resolution_clock::now();

  // method 2: Compute quad points in the beginning and use it when needed
  size_t num_quad_pts = 0;
  std::vector<fe::QuadData> quad_data;
  for (size_t e = 0; e < N; e++) {
    std::vector<fe::QuadData> qds = quad.getQuadPoints(nodes);
    if (e == 0)
      num_quad_pts = qds.size();
    for (auto qd : qds)
      quad_data.emplace_back(qd);
  }

  std::uint64_t t21 = hpx::util::high_resolution_clock::now();
  sum = 0.;
  for (size_t e = 0; e < N; e++) {
    for (size_t q = 0; q < num_quad_pts; q++) {
      fe::QuadData qd = quad_data[e * num_quad_pts + q];
      sum += qd.d_w * (qd.d_shapes[0] + qd.d_shapes[1] + qd.d_shapes[2]);
    }
  }
  std::uint64_t t22 = hpx::util::high_resolution_clock::now();

  if (n == 1 and N == 1000) {
    std::cout << "**********************************\n";
    std::cout << "Quadrature Time Efficiency Test\n";
    std::cout << "**********************************\n";
  }
  std::cout << "Quad order = " << n << ". Num Elements =  " << N << ".\n ";
  double dt_1 = double(t12 - t11) / 1.0e9;
  double dt_2 = double(t22 - t21) / 1.0e9;
  double perc = (dt_1 - dt_2) * 100. / dt_2;
  double qpt_mem = 13 * sizeof(double);
  double mem2 = double(quad_data.capacity() * qpt_mem) / double(1000000);
  std::cout << "  dt1 = " << dt_1 << ", dt2 = " << dt_2 << ", perc = " << perc
            << ". Mem saved = " << mem2 << " MB.\n";
}
