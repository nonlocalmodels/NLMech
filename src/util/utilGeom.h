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
#include <vector>

namespace util {

/*! @brief Provides geometrical methods such as point inside rectangle */
namespace geometry {

/*!
 * @brief Returns all corner points in the box
 * @param dim Dimension of the box
 * @param box Pair of points representing cuboid (rectangle in 2d)
 * @return Vector Vector of corner points
 */
std::vector<util::Point3> getCornerPoints(size_t dim, const
std::pair<util::Point3, util::Point3> &box);

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
 * @brief Checks if point is inside a cuboid (rectangle in 2-d, line in 1-d)
 * @param x Point
 * @param x_lbb Coordinate of left-bottom-back corner point
 * @param x_rtf Coordinate of right-top-front corner point
 * @return bool True if point inside cuboid, else false
 */
bool isPointInsideCuboid(size_t dim, util::Point3 x, util::Point3 x_lbb,
                         util::Point3 x_rtf);

/*!
 * @brief Computes angle between two vectors
 * @param vec_1 Vector 1
 * @param vec_2 Vector 2
 * @param anticlock Angle convention
 * @return angle Angle between vector 1 and 2
 */
double angle(util::Point3 vec_1, util::Point3 vec_2,
size_t dim, bool anticlock = false);

/*!
 * @brief Computes the area of triangle
 * @param nodes Vertices of the triangle
 * @return area Area of triangle
 */
double getTriangleArea(const std::vector<util::Point3> &nodes);

/*!
 * @brief Computes the volume of tetrahedron
 * @param nodes Vertices of the tetrahedron
 * @return area Area of tetrahedron
 */
double getTetVolume(const std::vector<util::Point3> &nodes);

/*!
 * @brief Computes the centroid of element
 * @param nodes Vertices of the element
 * @param elem_type Element type
 * @return Point Centroid
 */
util::Point3 getCenter(const std::vector<util::Point3> &nodes, const size_t
&elem_type);

/*!
 * @brief Computes the centroid and volume of element
 * @param nodes Vertices of the element
 * @param elem_type Element type
 * @return Point-volume Pair of centroid and volume
 */
std::pair<util::Point3, double> getCenterAndVol(const std::vector<util::Point3>
    &nodes, const size_t &elem_type);

} // namespace geometry

} // namespace util

#endif // UTIL_GEOMETRY_H
