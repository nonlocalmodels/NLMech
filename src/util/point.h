////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTIL_POINT_H
#define UTIL_POINT_H

#include "matrixBlaze.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
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
  template <class T>
  Point3(T x, T y, T z) : d_x(x), d_y(y), d_z(z){};

  /*!
   *  @brief Constructor
   *  @param x The coordinate vector
   */
  template <class T>
  explicit Point3(T x[3]) : d_x(x[0]), d_y(x[1]), d_z(x[2]){};

  /*!
   *  @brief Constructor
   *  @param p Point 
   */
  explicit Point3(const std::vector<double> &p) {

    if (p.empty())
      return;
    else if (p.size() == 1)
      d_x = p[0];
    else if (p.size() == 2) {
      d_x = p[0];
      d_y = p[1];
    } else if (p.size() == 3) {
      d_x = p[0];
      d_y = p[1];
      d_z = p[2];
    }
  }

  /*!
   *  @brief Copy constructor
   * @param p Point
   */
  Point3(const Point3 &p) : d_x(p.d_x), d_y(p.d_y), d_z(p.d_z) {};

  /*!
   * @brief Prints the information
   *
   * @param nt Number of tabs to append before printing
   * @param lvl Information level (higher means more information)
   * @return Manipulated string
   */
  std::string printStr(int nt = 0, int lvl = 0) const {

    std::string tabS = "";
    for (int i = 0; i < nt; i++)
      tabS += "\t";

    std::ostringstream oss;
    oss << tabS << "(" << d_x << ", " << d_y << ", " << d_z << ")";

    return oss.str();
  }


  /*!
   * @brief Prints the information
   *
   * @param nt Number of tabs to append before printing
   * @param lvl Information level (higher means more information)
   */
  void print(int nt = 0, int lvl = 0) const { std::cout << printStr(nt, lvl); }

  /*!
   * @brief Computes the Euclidean length of the vector
   * @return Length Euclidean length of the vector
   */
  double length() const { return std::sqrt(d_x * d_x + d_y * d_y + d_z * d_z); }

  /*!
   * @brief Computes the Euclidean length of the vector
   * @return Length Euclidean length of the vector
   */
  double lengthSq() const { return d_x * d_x + d_y * d_y + d_z *
  d_z; }

  /*
   * @brief Returns the unit vector
   * @return Vector Unit vector
   */
  Point3 unit() {
    auto l = this->length();
    if (l < 1.0E-10)
      return {};
    else
      return {d_x/l, d_y/l, d_z/l};
  }

  /*!
   * @brief Returns the unit vector
   * @return Vector Unit vector
   */
  Point3 unit() const {
    auto l = this->length();
    if (l < 1.0E-10)
      return {};
    else
      return {d_x/l, d_y/l, d_z/l};
  }

  /*!
   * @brief Computes the dot product of this vector with another point
   * @param b Another vector
   * @return Value a dot product
   */
  double dot(const Point3 &b) const { return d_x * b.d_x + d_y * b.d_y + d_z * b
  .d_z; }

  /*!
   * @brief Computes the distance between a given point from this point
   * @param b Another point
   * @return Value Distance between the two points
   */
  double dist(const Point3 &b) const {
    return std::sqrt((d_x - b.d_x) * (d_x - b.d_x) +
                     (d_y - b.d_y) * (d_y - b.d_y) +
                     (d_z - b.d_z) * (d_z - b.d_z));
  }

  /*!
   * @brief Computes the cross product between this vector and given vector
   * @param b Another vector
   * @return Vector Cross product
   */
  Point3 cross(const Point3 &b) const {
    return {-d_z * b.d_y + d_y * b.d_z,
            d_z * b.d_x - d_x * b.d_z,
            -d_y * b.d_x + d_x * b.d_y};
  }

  /*!
   * @brief Computes projection of vector on this vector
   * @param b Another vector
   * @param is_unit Is a unit vector
   * @return Vector Projection vector
   */
  Point3 project(const Point3 &b, bool is_unit = false) const {
    auto l_sq = (is_unit ? 1. : this->length() * this->length());
    auto dot = this->dot(b);
    return {dot * d_x / l_sq, dot * d_y / l_sq, dot * d_z / l_sq};
  }

  /*!
   * @brief Computes projection of vector on plane with normal as this vector
   * @param b Another vector
   * @param is_unit Is a unit vector
   * @return Vector Projection vector
   */
  Point3 projectNormal(const Point3 &b, bool is_unit = false) const {
    auto l_sq = (is_unit ? 1. : this->length() * this->length());
    auto dot = this->dot(b);
    return b - Point3(dot * d_x / l_sq, dot * d_y / l_sq, dot * d_z / l_sq);
  }


  /*!
   * @brief Computes the angle between vector given by this and the vector b
   * @param b Another vector
   * @return Value a dot product
   */
  double angle(Point3 b) const {

    auto ahat = this->unit();
    auto bhat = b.unit();
    return std::acos(ahat.dot(bhat));
  }

  /*!
   * @brief Computes the dot product of the vector and its transpose
   * @return dot(x,x.T)
   */
  util::Matrix33 toMatrix() {

    util::Matrix33 tmp = util::Matrix33();

    tmp(0, 0) = this->d_x * this->d_x;
    tmp(0, 1) = this->d_x * this->d_y;
    tmp(0, 2) = this->d_x * this->d_z;

    tmp(1, 0) = this->d_y * this->d_x;
    tmp(1, 1) = this->d_y * this->d_y;
    tmp(1, 2) = this->d_y * this->d_z;

    tmp(2, 0) = this->d_z * this->d_x;
    tmp(2, 1) = this->d_z * this->d_y;
    tmp(2, 2) = this->d_z * this->d_z;

    return tmp;
  }

  /*!
   * @brief Computes the dot product of the vector x and the transpose of y
   * @param b The node b
   * @return dot(x,y.T)
   */
  util::Matrix33 toMatrix(const util::Point3 b) {

    util::Matrix33 tmp = util::Matrix33();

    tmp(0, 0) = this->d_x * b.d_x;
    tmp(0, 1) = this->d_x * b.d_y;
    tmp(0, 2) = this->d_x * b.d_z;

    tmp(1, 0) = this->d_y * b.d_x;
    tmp(1, 1) = this->d_y * b.d_y;
    tmp(1, 2) = this->d_y * b.d_z;

    tmp(2, 0) = this->d_z * b.d_x;
    tmp(2, 1) = this->d_z * b.d_y;
    tmp(2, 2) = this->d_z * b.d_z;

    return tmp;
  }

  /*!
   * @brief Converts the point to an blaze vector
   * @return Point as blaze vector 3D
   */
  util::Vector3 toVector() {

    return util::Vector3{this->d_x, this->d_y, this->d_z};
  }

  /**
   * @name Group operators
   */
  /**@{*/


  /*!
  * @brief Adds two points
  * @param lhs Left-hand side
  * @param rhs Left-hand side
  * @return Sum of the two points
  */
  friend Point3 operator+(Point3 lhs, const Point3 &rhs) {
    lhs += rhs;
    return lhs;
  }

  /*!
  * @brief Substracts two points
  * @param lhs Left-hand side
  * @param rhs Left-hand side
  * @return Difference of the two points
  */
  friend Point3 operator-(Point3 lhs, const Point3 &rhs) {
    lhs -= rhs;
    return lhs;
  }

  /*!
  * @brief Multiplies a constant to the point
  * @param lhs Left-hand side
  * @param rhs Left-hand side
  * @return Scaled point
  */
  friend double operator*(Point3 lhs, const Point3 rhs) {
    return lhs.d_x * rhs.d_x + lhs.d_y * rhs.d_y + lhs.d_z * rhs.d_z;
  }

  /*!
  * @brief Adds two points
  * @param lhs Left-hand side
  * @param rhs Left-hand side
  * @return Sum of the two points
  */
  friend Point3 operator*(Point3 lhs, const double rhs) {
    lhs *= rhs;
    return lhs;
  }

  /*!
  * @brief Adds a scalar to the point
  * @param lhs Left-hand side
  * @param rhs Left-hand side
  * @return Shifted point
  */
  friend Point3 operator+(Point3 lhs, const double rhs) {
    return {lhs.d_x + rhs, lhs.d_y + rhs, lhs.d_z + rhs};
  }

    /*!
  * @brief Adds two points
  * @param lhs Left-hand side
  * @param rhs Left-hand side
  * @return Sum of the two points
  */
  friend Point3 operator+(const double lhs, Point3 rhs) {
    return {lhs + rhs.d_x, lhs + rhs.d_y, lhs + rhs.d_z};
  }

  /*!
  * @brief Substracts a scalar to the point
  * @param lhs Left-hand side
  * @param rhs Left-hand side
  * @return Shifted point
  */
  friend Point3 operator-(Point3 lhs, const double rhs) {
    return {lhs.d_x - rhs, lhs.d_y - rhs, lhs.d_z - rhs};
  }

  /*!
  * @brief Substracts a scalar to the point
  * @param lhs Left-hand side
  * @param rhs Left-hand side
  * @return Shifted point
  */
  friend Point3 operator-(const double lhs, Point3 rhs) {
    return {lhs - rhs.d_x, lhs - rhs.d_y, lhs - rhs.d_z};
  }

  /*!
  * @brief Substracts a scalar to the point
  * @param lhs Left-hand side
  * @param rhs Left-hand side
  * @return Shifted point
  */
  friend Point3 operator*(const double lhs, Point3 rhs) {
    rhs *= lhs;
    return rhs;
  }

  /*!
  * @brief Divedes a scalar to the point
  * @param lhs Left-hand side
  * @param rhs Left-hand side
  * @return Scaled point
  */
  friend Point3 operator/(Point3 lhs, const double rhs) {
    lhs /= rhs;
    return lhs;
  }

  /*!
  * @brief Add a scalar to the point
  * @param b Scalar factor
  * @return Scaled point
  */
  Point3 &operator+=(const double b) {

    d_x += b;
    d_y += b;
    d_z += b;
    return *this;
  }

  /*!
  * @brief Substracts a scalar to the point
  * @param b Scalar factor
  * @return Shifted point
  */
  Point3 &operator-=(const double b) {

    d_x -= b;
    d_y -= b;
    d_z -= b;
    return *this;
  }

  /*!
  * @brief Multiplies a scalar to the point
  * @param b Scalar factor
  * @return Scaled point
  */
  Point3 &operator*=(const double b) {

    d_x *= b;
    d_y *= b;
    d_z *= b;
    return *this;
  }

  /*!
  * @brief Adds a point 
  * @param b Point
  * @return Sum of the two points
  */
  Point3 &operator+=(const Point3 &b) {

    d_x += b.d_x;
    d_y += b.d_y;
    d_z += b.d_z;
    return *this;
  }

  /*!
  * @brief Substracts a point 
  * @param b Point
  * @return Difference of the two points
  */
  Point3 &operator-=(const Point3 &b) {

    d_x -= b.d_x;
    d_y -= b.d_y;
    d_z -= b.d_z;
    return *this;
  }

  /*!
  * @brief Multilpies a point 
  * @param b Point
  * @return Point-wise product of the two points
  */
  Point3 &operator*=(const Point3 &b) {

    d_x *= b.d_x;
    d_y *= b.d_y;
    d_z *= b.d_z;
    return *this;
  }

  /*!
  * @brief Divides a point by a scalr 
  * @param b Scalar factor
  * @return Scaled point
  */
  Point3 &operator/=(const double b) {

    d_x /= b;
    d_y /= b;
    d_z /= b;
    return *this;
  }

  /*!
  * @brief Access the 0th componnent of a point
  * @param i Index
  * @return The i-th index of the point
  */
  double &operator[](size_t i) {

    if (i == 0)
      return d_x;
    else if (i == 1)
      return d_y;
    else
      return d_z;
  }

  /*!
  * @brief Access the 0th componnent of a point
  * @param i Index
  * @return The i-th index of the point
  */
  const double &operator[](size_t i) const {

    if (i == 0)
      return d_x;
    else if (i == 1)
      return d_y;
    else
      return d_z;
  }
  
  /*!
  * @brief Print the point's vlaue to the standard output stream
  * @param os Standard output stream
  * @param p Point
  * @return The standard out stream
  */ 
  friend std::ostream &operator<<(std::ostream &os, const Point3 p);

  /** @}*/
};

inline std::ostream &operator<<(std::ostream &os, const Point3 p) {
  os << p.d_x << " " << p.d_y << " " << p.d_z << std::endl;
  return os;
}

} // namespace util

#endif
