// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef TRANSFORMATION_H
#define TRANSFORMATION_H

#include <vector>

namespace util {

/*! @brief Provides comparison of floating point values  */
namespace transformation {

/*!
 * @brief Rotates a vector in xy-plane in clockwise direction
 * @param x Point
 * @param theta Angle
 * @return Point after rotation
 */
std::vector<double> rotateCW2D(const std::vector<double> x, const double theta);

/*!
 * @brief Rotates a vector in xy-plane in anti-clockwise direction
 * @param x Point
 * @param theta Angle
 * @return Point after rotation
 */
std::vector<double> rotateACW2D(const std::vector<double> x,
                                const double theta);

} // namespace transformation

} // namespace util

#endif
