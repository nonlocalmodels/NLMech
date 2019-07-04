////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTIL_GEOMETRY_H
#define UTIL_GEOMETRY_H

#include "point.h"              // definition of Point3

namespace util {

/*! @brief Provides geometrical methods such as point inside rectangle */
namespace geometry {

/*!
 * @brief Checks if point is inside a rectangle
 * @param x Point
 * @param x_min X coordinate of left-bottom corner point
 * @param x_max X coordinate of right-top corner point
 * @param y_min Y coordinate of left-bottom corner point
 * @param y_max Y coordinate of right-top corner point
 * @return bool True if point inside rectangle, else false
 */
bool isPointInsideRectangle(util::Point3 x, double x_min, double x_max,
                            double y_min, double y_max);

/*!
 * @brief Checks if point is inside an angled rectangle
 * @param x Point
 * @param x_min X coordinate of left-bottom corner point
 * @param x_max X coordinate of right-top corner point
 * @param y_min Y coordinate of left-bottom corner point
 * @param y_max Y coordinate of right-top corner point
 * @param theta Angle of orientation of rectangle from x-axis
 * @return bool True if point inside rectangle, else false
 */
bool isPointInsideAngledRectangle(util::Point3 x, double x_min, double x_max,
                                  double y_min, double y_max, double theta);

} // namespace geometry

} // namespace util

#endif // UTIL_GEOMETRY_H
