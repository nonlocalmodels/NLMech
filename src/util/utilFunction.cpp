////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "utilFunction.h"
#include "compare.h"              // compare functions
#include <cmath>                  // definition of sin, cosine etc
#include <iostream>               // cerr

double util::function::hatFunction(const double &x, const double &x_min,
                                   const double &x_max) {

  if (util::compare::definitelyGreaterThan(x, x_min - 1.0E-12) and
      util::compare::definitelyLessThan(x, x_max + 1.0E-12)) {

    double x_mid = 0.5 * (x_min + x_max);
    double l = x_mid - x_min;

    // check if this is essentially a point load (dirac)
    if (l < 1.0E-12)
      return 1.0;

    if (util::compare::definitelyLessThan(x, x_mid))
      return (x - x_min) / l;
    else
      return (x_max - x) / l;
  } else
    return 0.0;
}

double util::function::hatFunctionQuick(const double &x, const double &x_min,
                                        const double &x_max) {

  double x_mid = 0.5 * (x_min + x_max);
  double l = x_mid - x_min;

  // check if this is essentially a point load (dirac)
  if (l < 1.0E-12)
    return 1.0;

  if (util::compare::definitelyLessThan(x, x_mid))
    return (x - x_min) / l;
  else
    return (x_max - x) / l;
}

double util::function::linearStepFunc(const double &x, const double &x1,
                                      const double &x2) {

  double period = std::floor(x / (x1 + x2));

  if (util::compare::definitelyLessThan(x, period * (x1 + x2) + x1))
    return x - period * x2;
  else
    return (period + 1.) * x1;
}

double util::function::derLinearStepFunc(const double &x, const double &x1,
                                      const double &x2) {

  double period = std::floor(x / (x1 + x2));

  if (util::compare::definitelyLessThan(x, period * (x1 + x2) + x1))
    return 1.;
  else
    return 0.;
}

double util::function::gaussian(const double &r, const double &a,
                                const double &beta) {
  return a * std::exp(-std::pow(r, 2) / beta);
}

double util::function::gaussian2d(const util::Point3 &x, const size_t &dof,
                                  const std::vector<double> &params) {

  if (params.size() < 6) {
    std::cerr << "Error: Not enough parameters to compute guassian 2-d "
                 "function.\n";
    exit(1);
  }

  return util::function::gaussian(
             x.dist(util::Point3(params[0], params[1], 0.)), params[5],
             params[4]) *
         params[2 + dof];
}

double util::function::doubleGaussian2d(const util::Point3 &x,
                                        const size_t &dof,
                                        const std::vector<double> &params) {

  if (params.size() < 10) {
    std::cerr << "Error: Not enough parameters to compute guassian 2-d "
                 "function.\n";
    exit(1);
  }

  return util::function::gaussian(
             x.dist(util::Point3(params[0], params[1], 0.)), params[9],
             params[8]) *
             params[4 + dof] +
         util::function::gaussian(
             x.dist(util::Point3(params[2], params[3], 0.)), params[9],
             params[8]) *
             params[6 + dof];
}

util::Point3 util::function::signVector(const util::Point3 &v) {

  auto r = util::Point3(1., 1., 1.);
  if (util::compare::definitelyLessThan(v.d_x, 0.))
    r.d_x = -1.;
  if (util::compare::definitelyLessThan(v.d_y, 0.))
    r.d_y = -1.;
  if (util::compare::definitelyLessThan(v.d_z, 0.))
    r.d_z = -1.;

  return r;
}

double util::function::getDeterminant(const std::vector<util::Point3> &rows) {

  double a =
      rows[0].d_x * (rows[1].d_y * rows[2].d_z - rows[1].d_z * rows[2].d_y);
  double b =
      rows[0].d_y * (rows[1].d_x * rows[2].d_z - rows[1].d_z * rows[2].d_x);
  double c =
      rows[0].d_z * (rows[1].d_x * rows[2].d_y - rows[1].d_y * rows[2].d_x);

  return a - b + c;
}
