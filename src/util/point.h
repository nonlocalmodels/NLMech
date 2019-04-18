// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef UTIL_POINT_H
#define UTIL_POINT_H

#include <cmath>
#include <iomanip>
#include <stdlib.h>

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
   * @brief Computes the Euclidean length of the vector
   * @return Length Euclidean length of the vector
   */
  double length() const { return std::sqrt(d_x * d_x + d_y * d_y + d_z * d_z); }

  /*!
   * @brief Computes the dot product of this vector with another point
   * @param b Another vector
   * @return Value a dot product
   */
  double dot(Point3 b) { return d_x * b.d_x + d_y * b.d_y + d_z * b.d_z; }

  /*!
   * @brief Computes the dot product of this vector with another point
   * @param b Another vector
   * @return Value a dot product
   */
  double dot(Point3 b) const { return d_x * b.d_x + d_y * b.d_y + d_z * b.d_z; }

  /*!
   * @brief Computes the distance between a given point from this point
   * @param b Another point
   * @return Value Distance between the two points
   */
  double dist(Point3 b) {
    return std::sqrt((d_x - b.d_x) * (d_x - b.d_x) +
                     (d_y - b.d_y) * (d_y - b.d_y) +
                     (d_z - b.d_z) * (d_z - b.d_z));
  }

  /*!
   * @brief Computes the distance between a given point from this point
   * @param b Another point
   * @return Value Distance between the two points
   */
  double dist(Point3 b) const {
    return std::sqrt((d_x - b.d_x) * (d_x - b.d_x) +
                     (d_y - b.d_y) * (d_y - b.d_y) +
                     (d_z - b.d_z) * (d_z - b.d_z));
  }

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

} // namespace util

#endif
