// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef UTIL_POINT_H
#define UTIL_POINT_H

#include <cmath>
#include <iomanip>
#include <stdlib.h>
#include <vector>

/*!
 * @brief Collection of methods useful in simulation
 *
 * This namespace provides number of useful functions and struct definition.
 *
 * @sa Point3, Matrix3, SymMatrix3, compare, transformation
 */
namespace util {

/*! @brief A structure to represent 3d vectors */
struct Point3 {

  /*! @brief the x coordinate */
  double d_x;

  /*! @brief the y coordinate */
  double d_y;

  /*! @brief the z coordinate */
  double d_z;

  /*!
   * @brief Constructor
   */
  Point3() : d_x(0.), d_y(0.), d_z(0.){};

  /*!
   *  @brief Constructor
   *  @param x The x coordinate
   *  @param y The y coordinate
   *  @param z The z coordinate
   */
  Point3(double x, double y, double z) : d_x(x), d_y(y), d_z(z){};

  /*!
   * @brief Computes the Euclidean length of the vector
   * @return Length Euclidean length of the vector
   */
  double length() { return std::sqrt(d_x * d_x + d_y * d_y + d_z * d_z); }

  /*!
   * @brief Computes the dot product of this vector with another point
   * @param b Another vector
   * @return Value a dot product
   */
  double dot(Point3 b) { return d_x * b.d_x + d_y * b.d_y + d_z * b.d_z; }

  /**
   * @name Group operators
   */
  /**@{*/

  friend Point3 operator+(Point3 lhs, const Point3 &rhs) {
    lhs += rhs;
    return lhs;
  }

  friend Point3 operator-(Point3 lhs, const Point3 &rhs) {
    lhs -= rhs;
    return lhs;
  }

  friend Point3 operator*(Point3 lhs, const double rhs) {
    lhs *= rhs;
    return lhs;
  }

  friend Point3 operator/(Point3 lhs, const double rhs) {
    lhs /= rhs;
    return lhs;
  }

  Point3 &operator+=(const double b) {

    d_x += b;
    d_y += b;
    d_z += b;
    return *this;
  }

  Point3 &operator-=(const double b) {

    d_x -= b;
    d_y -= b;
    d_z -= b;
    return *this;
  }

  Point3 &operator*=(const double b) {

    d_x *= b;
    d_y *= b;
    d_z *= b;
    return *this;
  }

  Point3 &operator+=(const Point3 &b) {

    d_x += b.d_x;
    d_y += b.d_y;
    d_z += b.d_z;
    return *this;
  }

  Point3 &operator-=(const Point3 &b) {

    d_x -= b.d_x;
    d_y -= b.d_y;
    d_z -= b.d_z;
    return *this;
  }

  Point3 &operator*=(const Point3 &b) {

    d_x *= b.d_x;
    d_y *= b.d_y;
    d_z *= b.d_z;
    return *this;
  }

  Point3 &operator/=(const double b) {

    d_x /= b;
    d_y /= b;
    d_z /= b;
    return *this;
  }

  double &operator[](size_t i) {

    if (i == 0)
      return d_x;
    else if (i == 1)
      return d_y;
    else
      return d_z;
  }

  /** @}*/
};

/*!
 * @brief Adds point b to point a inplace
 * @param a Point a
 * @param b Point b
 */
static inline void addInplace(std::vector<util::Point3> &a,
                              std::vector<util::Point3> b) {

  for (size_t i = 0; i < a.size(); i++) {

    a[i].d_x = a[i].d_x + b[i].d_x;
    a[i].d_y = a[i].d_y + b[i].d_y;
    a[i].d_z = a[i].d_z + b[i].d_z;
  }
}

/*!
 * @brief Subtracts point b to point a inplace
 * @param a Point a
 * @param b Point b
 */
static inline void subInplace(std::vector<util::Point3> &a,
                              std::vector<util::Point3> b) {

  for (size_t i = 0; i < a.size(); i++) {

    a[i].d_x = a[i].d_x - b[i].d_x;
    a[i].d_y = a[i].d_y - b[i].d_y;
    a[i].d_z = a[i].d_z - b[i].d_z;
  }
}

/*!
 * @brief Copys vector's a content to b
 * @param a vector of points
 * @param b vector of points
 */
inline void copy(std::vector<util::Point3> a, std::vector<util::Point3> &b) {

  b.clear();
  std::copy(a.begin(), a.end(), std::back_inserter(b));
}

/*!
 * @brief Adds vector of points b to vector of points a in a element wise
 * manner
 * @param a vector of points
 * @param b vector of points
 * @return Vector vector of point after element wise addition
 */
inline std::vector<util::Point3> add(std::vector<util::Point3> a,
                                     std::vector<util::Point3> b) {

  std::vector<util::Point3> res;
  for (size_t i = 0; i < a.size(); i++)
    res.push_back(a[i] + b[i]);

  return res;
}

/*!
 * @brief Subtracts vector of points b to vector of points a in a
 * element-wise manner
 * @param a Point a
 * @param b Point b
 * @return Vector vector of points after element-wise subtraction
 */
inline std::vector<util::Point3> subb(std::vector<util::Point3> a,
                                      std::vector<util::Point3> b) {

  std::vector<util::Point3> res;
  for (size_t i = 0; i < a.size(); i++)
    res.push_back(a[i] - b[i]);

  return res;
}

} // namespace util

#endif
