// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "transfomation.h"
#include <cmath> // definition of sin, cosine etc

std::vector<double>
util::transformation::rotateCW2D(const std::vector<double> x,
                                 const double theta) {

  return std::vector<double>{x[0] * std::cos(theta) + x[1] * std::sin(theta),
                             -x[0] * std::sin(theta) + x[1] * std::cos(theta),
                             0.0};
}

util::Point3
util::transformation::rotateCW2D(const util::Point3 x,
                                 const double theta) {

  return {x.d_x * std::cos(theta) + x.d_y * std::sin(theta),
          -x.d_x * std::sin(theta) + x.d_y * std::cos(theta), 0.0};
}

std::vector<double>
util::transformation::rotateACW2D(const std::vector<double> x,
                                  const double theta) {

  return rotateCW2D(x, -theta);
}

util::Point3
util::transformation::rotateACW2D(const util::Point3 x,
                                  const double theta) {

  return rotateCW2D(x, -theta);
}