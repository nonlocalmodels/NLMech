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
 * @brief Checks if point is inside a rectangle
 * @param x Point
 * @param x_lb Coordinate of left-bottom corner point
 * @param x_rt Coordinate of right-top corner point
 * @return bool True if point inside rectangle, else false
 */
bool isPointInsideRectangle(util::Point3 x, util::Point3 x_lb,
        util::Point3 x_rt);

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

/*!
 * @brief Checks if point is inside a Cuboid
 * @param x Point
 * @param x_min X coordinate of left-bottom-back corner point
 * @param x_max X coordinate of right-top-front corner point
 * @param y_min Y coordinate of left-bottom-back corner point
 * @param y_max Y coordinate of right-top-front corner point
 * @param z_min Z coordinate of left-bottom-back corner point
 * @param z_max Z coordinate of right-top-front corner point
 * @return bool True if point inside
 */
bool isPointInsideCuboid(util::Point3 x, double x_min, double x_max,
                            double y_min, double y_max, double z_min, double
                            z_max);

/*!
 * @brief Checks if point is inside a Cuboid
 * @param x Point
 * @param x_lbb Coordinate of left-bottom-back corner point
 * @param x_rtf Coordinate of right-top-front corner point
 * @return bool True if point inside
 */
bool isPointInsideCuboid(util::Point3 x, util::Point3 x_lbb,
                            util::Point3 x_rtf);

/*!
 * @brief Computes the area of triangle
 * @param v1 First vertex
 * @param v2 Second vertex
 * @param v3 Third vertex
 * @return area Triangle area
 */
double triangleArea(const util::Point3 &v1, const util::Point3 &v2,
                    const util::Point3 &v3);

} // namespace geometry

} // namespace util

#endif // UTIL_GEOMETRY_H
