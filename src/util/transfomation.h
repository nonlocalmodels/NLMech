////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTIL_TRANSFORMATION_H
#define UTIL_TRANSFORMATION_H

#include "point.h"              // definition of Point3
#include <vector>

namespace util {

/*! @brief Provides geometrical methods such rotation of a vector */
namespace transformation {

/*!
 * @brief Rotates a vector in xy-plane in clockwise direction
 * @param x Point
 * @param theta Angle
 * @return Point after rotation
 */
std::vector<double> rotateCW2D(const std::vector<double> &x,
                               const double &theta);

/*!
 * @brief Rotates a vector in xy-plane in clockwise direction
 * @param x Point
 * @param theta Angle
 * @return Point after rotation
 */
util::Point3 rotateCW2D(const util::Point3 &x, const double &theta);

/*!
 * @brief Rotates a vector in xy-plane in anti-clockwise direction
 * @param x Point
 * @param theta Angle
 * @return Point after rotation
 */
std::vector<double> rotateACW2D(const std::vector<double> &x,
                                const double &theta);

/*!
 * @brief Rotates a vector in xy-plane in anti-clockwise direction
 * @param x Point
 * @param theta Angle
 * @return Point after rotation
 */
util::Point3 rotateACW2D(const util::Point3 &x, const double &theta);

} // namespace transformation

} // namespace util

#endif // UTIL_TRANSFORMATION_H
