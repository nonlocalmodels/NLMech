////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
// ////////////////////////////////////////////////////////////////////////////////

#include "utilGeom.h"
#include "compare.h"
#include "fe/quadElem.h"
#include "feElementDefs.h"
#include "transfomation.h"
#include "utilFunction.h"
#include <cmath>  // definition of sin, cosine etc
#include <iostream>

std::vector<util::Point3> util::geometry::getCornerPoints(
    size_t dim, const std::pair<util::Point3, util::Point3> &box) {
  if (dim == 1)
    return {box.first, box.second};
  else if (dim == 2)
    return {box.first, util::Point3(box.second.d_x, box.first.d_y, 0.),
            box.second, util::Point3(box.first.d_x, box.second.d_y, 0.)};
  else if (dim == 3) {
    double a = box.second.d_x - box.first.d_x;
    double b = box.second.d_y - box.first.d_y;
    double c = box.second.d_z - box.first.d_z;
    return {box.first,
            box.first + util::Point3(a, 0., 0.),
            box.first + util::Point3(a, b, 0.),
            box.first + util::Point3(0., b, 0.),
            box.first + util::Point3(0., 0., c),
            box.first + util::Point3(a, 0., c),
            box.second,
            box.first + util::Point3(0., b, c)};
  } else {
    std::cerr << "Error: Check dimension = " << dim << ".\n";
    exit(1);
  }
}

bool util::geometry::isPointInsideRectangle(util::Point3 x, double x_min,
                                            double x_max, double y_min,
                                            double y_max) {
  return !(util::compare::definitelyLessThan(x.d_x, x_min - 1.0E-12) or
           util::compare::definitelyLessThan(x.d_y, y_min - 1.0E-12) or
           util::compare::definitelyGreaterThan(x.d_x, x_max + 1.0E-12) or
           util::compare::definitelyGreaterThan(x.d_y, y_max + 1.0E-12));
}

bool util::geometry::isPointInsideRectangle(util::Point3 x, util::Point3 x_lb,
                                            util::Point3 x_rt) {
  return !(util::compare::definitelyLessThan(x.d_x, x_lb.d_x - 1.0E-12) or
           util::compare::definitelyLessThan(x.d_y, x_lb.d_y - 1.0E-12) or
           util::compare::definitelyGreaterThan(x.d_x, x_rt.d_x + 1.0E-12) or
           util::compare::definitelyGreaterThan(x.d_y, x_rt.d_y + 1.0E-12));
}

bool util::geometry::isPointInsideAngledRectangle(util::Point3 x, double x1,
                                                  double x2, double y1,
                                                  double y2, double theta) {
  // we assume that the rectangle has passed the test

  //
  //                             (x2,y2)
  //                            o
  //
  //
  //
  //
  //
  //        o
  //      (x1,y1)

  // get divisors
  util::Point3 lam = util::transformation::rotateCW2D(
      util::Point3(x2 - x1, y2 - y1, 0.0), theta);

  // double lam1 = (x2-x1) * std::cos(theta) + (y2-y1) * std::sin(theta);
  // double lam2 = -(x2-x1) * std::sin(theta) + (y2-y1) * std::cos(theta);

  // get mapped coordinate of x
  util::Point3 xmap = util::transformation::rotateCW2D(
      util::Point3(x[0] - x1, x[1] - y1, 0.0), theta);

  // double xmap = (x[0]-x1) * std::cos(theta) + (x[1]-y1) * std::sin(theta);
  // double ymap = -(x[0]-x1) * std::sin(theta) + (x[1]-y1) * std::cos(theta);

  // check if mapped coordinate are out of range [0, lam1] and [0, lam2]
  return !(util::compare::definitelyLessThan(xmap[0], -1.0E-12) or
           util::compare::definitelyLessThan(xmap[1], -1.0E-12) or
           util::compare::definitelyGreaterThan(xmap[0], lam[0] + 1.0E-12) or
           util::compare::definitelyGreaterThan(xmap[1], lam[1] + 1.0E-12));
}

bool util::geometry::isPointInsideCuboid(size_t dim, util::Point3 x,
                                         util::Point3 x_lbb,
                                         util::Point3 x_rtf) {
  if (dim == 1)
    return !(util::compare::definitelyLessThan(x.d_x, x_lbb.d_x - 1.0E-12) or
             util::compare::definitelyGreaterThan(x.d_x, x_rtf.d_x + 1.0E-12));
  else if (dim == 2)
    return !(util::compare::definitelyLessThan(x.d_x, x_lbb.d_x - 1.0E-12) or
             util::compare::definitelyLessThan(x.d_y, x_lbb.d_y - 1.0E-12) or
             util::compare::definitelyGreaterThan(x.d_x, x_rtf.d_x + 1.0E-12) or
             util::compare::definitelyGreaterThan(x.d_y, x_rtf.d_y + 1.0E-12));
  else
    return !(util::compare::definitelyLessThan(x.d_x, x_lbb.d_x - 1.0E-12) or
             util::compare::definitelyLessThan(x.d_y, x_lbb.d_y - 1.0E-12) or
             util::compare::definitelyLessThan(x.d_z, x_lbb.d_z - 1.0E-12) or
             util::compare::definitelyGreaterThan(x.d_x, x_rtf.d_x + 1.0E-12) or
             util::compare::definitelyGreaterThan(x.d_y, x_rtf.d_y + 1.0E-12) or
             util::compare::definitelyGreaterThan(x.d_z, x_rtf.d_z + 1.0E-12));
}

double util::geometry::angle(util::Point3 vec_1, util::Point3 vec_2, size_t dim,
                             bool anticlock) {
  if (dim == 2) {
    double dot = vec_1.dot(vec_2);
    double cross = vec_1.d_x * vec_2.d_y - vec_1.d_y * vec_2.d_x;

    double angle = std::atan2(cross, dot);

    if (!anticlock) return angle;

    if (angle < 0.)
      return angle + 2. * M_PI;
    else
      return angle;

  } else if (dim == 3) {
    // TODO
    //  Implement angle between vectors in 3-d
    return 0.;
  }
}

double util::geometry::getTriangleArea(const std::vector<util::Point3> &nodes) {
  return 0.5 * ((nodes[1].d_x - nodes[0].d_x) * (nodes[2].d_y - nodes[0].d_y) -
                (nodes[2].d_x - nodes[0].d_x) * (nodes[1].d_y - nodes[0].d_y));
}

double util::geometry::getTetVolume(const std::vector<util::Point3> &nodes) {
  auto a = nodes[1] - nodes[0];
  auto b = nodes[2] - nodes[0];
  auto c = nodes[3] - nodes[0];

  return (1. / 6.) * util::function::getDeterminant({a, b, c});
}

util::Point3 util::geometry::getCenter(const std::vector<util::Point3> &nodes,
                                       const size_t &elem_type) {
  // process different element types

  if (elem_type == util::vtk_type_line) {
    if (nodes.size() != 2) {
      std::cerr << "Error: Number of nodes and the element type does not "
                   "match.\n";
      exit(1);
    }

    return {0.5 * nodes[0].d_x + 0.5 * nodes[1].d_x,
            0.5 * nodes[0].d_y + 0.5 * nodes[1].d_y,
            0.5 * nodes[0].d_z + 0.5 * nodes[1].d_z};

  } else if (elem_type == util::vtk_type_triangle) {
    if (nodes.size() != 3) {
      std::cerr << "Error: Number of nodes and the element type does not "
                   "match.\n";
      exit(1);
    }

    auto center = util::Point3();
    for (const auto &node : nodes)
      center += util::Point3(1. * node.d_x / 3., 1. * node.d_y / 3.,
                             1. * node.d_z / 3.);

    return center;

  } else if (elem_type == util::vtk_type_quad) {
    if (nodes.size() != 4) {
      std::cerr << "Error: Number of nodes and the element type does not "
                   "match.\n";
      exit(1);
    }

    auto quad_elem = fe::QuadElem(1);
    auto quads = quad_elem.getQuadPoints(nodes);

    return quads[0].d_p;

  } else if (elem_type == util::vtk_type_tetra) {
    if (nodes.size() != 4) {
      std::cerr << "Error: Number of nodes and the element type does not "
                   "match.\n";
      exit(1);
    }

    auto center = util::Point3();
    for (const auto &node : nodes)
      center += util::Point3(1. * node.d_x / 4., 1. * node.d_y / 4.,
                             1. * node.d_z / 4.);

    return center;

  } else {
    std::cerr << "Error: Element type " << elem_type
              << " currently not supported in getCenter() method.\n";
    exit(1);
  }
}

std::pair<util::Point3, double> util::geometry::getCenterAndVol(
    const std::vector<util::Point3> &nodes, const size_t &elem_type) {
  // process different element types

  if (elem_type == util::vtk_type_line) {
    if (nodes.size() != 2) {
      std::cerr << "Error: Number of nodes = " << nodes.size()
                << " and the element type = " << elem_type
                << " does not match.\n";
      exit(1);
    }

    return std::make_pair(util::Point3(0.5 * nodes[0].d_x + 0.5 * nodes[1].d_x,
                                       0.5 * nodes[0].d_y + 0.5 * nodes[1].d_y,
                                       0.5 * nodes[0].d_z + 0.5 * nodes[1].d_z),
                          nodes[0].dist(nodes[1]));

  } else if (elem_type == util::vtk_type_triangle) {
    if (nodes.size() != 3) {
      std::cerr << "Error: Number of nodes = " << nodes.size()
                << " and the element type = " << elem_type
                << " does not match.\n";
      exit(1);
    }

    auto center = util::Point3();
    for (const auto &node : nodes)
      center += util::Point3(1. * node.d_x / 3., 1. * node.d_y / 3.,
                             1. * node.d_z / 3.);

    return std::make_pair(center, std::abs(getTriangleArea(nodes)));

  } else if (elem_type == util::vtk_type_quad) {
    if (nodes.size() != 4) {
      std::cerr << "Error: Number of nodes = " << nodes.size()
                << " and the element type = " << elem_type
                << " does not match.\n";
      exit(1);
    }

    auto quad_elem = fe::QuadElem(1);
    auto quads = quad_elem.getQuadPoints(nodes);

    return std::make_pair(quads[0].d_p, quads[0].d_w);

  } else if (elem_type == util::vtk_type_tetra) {
    if (nodes.size() != 4) {
      std::cerr << "Error: Number of nodes = " << nodes.size()
                << " and the element type = " << elem_type
                << " does not match.\n";
      exit(1);
    }

    auto center = util::Point3();
    for (const auto &node : nodes)
      center += util::Point3(1. * node.d_x / 4., 1. * node.d_y / 4.,
                             1. * node.d_z / 4.);

    return std::make_pair(center, std::abs(getTetVolume(nodes)));

  } else {
    std::cerr << "Error: Element type " << elem_type
              << " currently not supported in getCenterAndVol() method.\n";
    exit(1);
  }
}

bool util::geometry::doLinesIntersect(util::Point3 A, util::Point3 B,
                                      util::Point3 C, util::Point3 D) {
  // four direction for two lines and points of other line
  int dir1 = direction(A, B, C);
  int dir2 = direction(A, B, D);
  int dir3 = direction(C, D, A);
  int dir4 = direction(C, D, B);

  if (dir1 != dir2 && dir3 != dir4) return true;  // they are intersecting

  if (dir1 == 0 && onLine(A, B, C))  // when p2 of line2 are on the line1
    return true;

  if (dir2 == 0 && onLine(A, B, D))  // when p1 of line2 are on the line1
    return true;

  if (dir3 == 0 && onLine(C, D, A))  // when p2 of line1 are on the line2
    return true;

  if (dir4 == 0 && onLine(C, D, B))  // when p1 of line1 are on the line2
    return true;

  return false;
}

bool util::geometry::onLine(util::Point3 A, util::Point3 B, util::Point3 C) {
  if (C.d_x <= std::max(A.d_x, B.d_x) && C.d_x <= std::min(A.d_x, B.d_x) &&
      (C.d_y <= std::max(A.d_y, B.d_y) && C.d_y <= std::min(A.d_y, B.d_y)))
    return true;

  return false;
}

int util::geometry::direction(util::Point3 A, util::Point3 B, util::Point3 C) {
  int val =
      (B.d_y - A.d_y) * (C.d_x - B.d_x) - (B.d_x - A.d_x) * (C.d_y - B.d_y);
  if (val == 0)
    return 0;  // colinear
  else if (val < 0)
    return 2;  // anti-clockwise direction
  return 1;    // clockwise direction
}

bool util::geometry::isPointinCircle(util::Point3 A, util::Point3 center,
                                     double radius) {
  std::cout << A.d_x << " " << A.d_y << " " << center.d_x << " " << center.d_y
            << " " << radius << " " << (A - center).length() << std::endl;

  if ((A - center).length() <= radius) return true;

  return false;
}