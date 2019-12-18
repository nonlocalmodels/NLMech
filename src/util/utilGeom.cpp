////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
// ////////////////////////////////////////////////////////////////////////////////

#include "utilGeom.h"
#include "compare.h"
#include "transfomation.h"
#include <cmath>                      // definition of sin, cosine etc
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
  }
  else {
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

bool util::geometry::isPointInsideCuboid(size_t dim, util::Point3 x, util::Point3 x_lbb,
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

double util::geometry::angle(util::Point3 vec_1, util::Point3 vec_2, size_t dim, bool anticlock) {

  if (dim == 2) {

    double dot = vec_1.dot(vec_2);
    double cross = vec_1.d_x * vec_2.d_y - vec_1.d_y * vec_2.d_x;

    double angle =  std::atan2(cross, dot);

    if (!anticlock)
      return angle;

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
