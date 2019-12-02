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

/*!
 * @brief Returns the list of elements in anti-clockwise order starting from
 * given index i
 * @param i Index i which will be the first element of returned list
 * @param n Number of elements
 * @return list List of indices starting from i in ACW order
 */
std::vector<size_t> cyclicOrderACW(const size_t &i, const size_t &n);

/*!
 * @brief Returns the list of elements in anti-clockwise order starting from
 * given index i and then index j
 * @param i Index i which will be the first element of returned list
 * @param j Index j which will be the second element of returned list
 * @param n Number of elements
 * @return list List of indices starting from i and j in ACW order
 */
std::vector<size_t> cyclicOrderACW(const size_t &i, const size_t &j, const size_t &n);

/*!
 * @brief Returns the list of elements in anti-clockwise order with indices
 * i, j, and k as the first three elements
 * @param i Index i which will be the first element of returned list
 * @param j Index j which will be the second element of returned list
 * @param k Index k which will be the third element of returned list
 * @param n Number of elements
 * @return list List of indices starting from i, j and k in ACW order
 */
std::vector<size_t> cyclicOrderACW(const size_t &i, const size_t &j, const size_t &k, const size_t &n);

} // namespace transformation

} // namespace util

#endif // UTIL_TRANSFORMATION_H
