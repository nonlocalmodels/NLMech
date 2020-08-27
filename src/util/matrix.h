////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTIL_MATRIX_H
#define UTIL_MATRIX_H

#include "point.h"

namespace util {

/*! @brief A structure to represent 3d matrices */
struct Matrix3 {

  /**
   * @name Data members
   */
  /**@{*/

  /*! @brief data */
  float d_data[3][3]{};

  /** @}*/

  /*!
   * @brief Constructor
   */
  Matrix3() {

    d_data[0][0] = 0.;
    d_data[0][1] = 0.;
    d_data[0][2] = 0.;

    d_data[1][0] = 0.;
    d_data[1][1] = 0.;
    d_data[1][2] = 0.;

    d_data[2][0] = 0.;
    d_data[2][1] = 0.;
    d_data[2][2] = 0.;
  };

  /*!
   * @brief Constructor
   *
   * @param diagonal Diagonal vector
   */
  Matrix3(const util::Point3& diagonal) {

    d_data[0][0] = diagonal.d_x;
    d_data[0][1] = 0.;
    d_data[0][2] = 0.;

    d_data[1][0] = 0.;
    d_data[1][1] = diagonal.d_y;
    d_data[1][2] = 0.;

    d_data[2][0] = 0.;
    d_data[2][1] = 0.;
    d_data[2][2] = diagonal.d_z;
  };

  /*!
   * @brief Constructor
   *
   * @param a1 first row
   * @param a2 second row
   * @param a3 third row
   */
  Matrix3(const util::Point3& a1, const util::Point3& a2, const util::Point3& a3) {

    d_data[0][0] = a1.d_x;
    d_data[0][1] = a1.d_y;
    d_data[0][2] = a1.d_z,
        d_data[1][0] = a2.d_x;
    d_data[1][1] = a2.d_y;
    d_data[1][2] = a2.d_z;
    d_data[2][0] = a3.d_x;
    d_data[2][1] = a3.d_y;
    d_data[2][2] = a3.d_z;
  };

  /*!
   * @brief Constructor
   *
   * @param m Matrix in vector template
   */
  Matrix3(const std::vector<std::vector<double>> &m) {
    for (size_t i=0; i<3; i++)
      for (size_t j=0; j<3; j++)
        d_data[i][j] = m[i][j];
  }

  /*!
   * @brief Constructor
   *
   * @param m Matrix in vector template
   */
  Matrix3(const Matrix3 &m) {
    for (size_t i=0; i<3; i++)
      for (size_t j=0; j<3; j++)
        d_data[i][j] = m(i,j);
  }

  /*!
   * @brief Prints the information
   *
   * @param nt Number of tabs to append before printing
   * @param lvl Information level (higher means more information)
   * @return The resulting string
   */
  std::string printStr(int nt = 0, int lvl = 0) const {

    std::string tabS = "";
    for (int i = 0; i < nt; i++)
      tabS += "\t";

    std::ostringstream oss;
    for (size_t i=0; i<3; i++)
      oss << tabS << "[" << (*this)(i, 0) << ", " << (*this)(i, 1) << ", "
          << (*this)(i, 2) << "]" << std::endl;
    oss << std::endl;

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
   * @brief Access operator of the i-th element
   * @param i Index i
   * @return The i-th element
   */
  Point3 operator()(size_t i) {
    return Point3(d_data[i]);
  }

  /*!
   * @brief Access operator of the i-th element
   * @param i Index i
   * @return The i-th element
   */
  Point3 operator()(size_t i) const {
    return Point3(d_data[i]);
  }

  /*!
   * @brief Access operator of the Matrix element M(i,j)
   * @param i Index i
   * @param j Index j
   * @return Matrix element M(i,j)
   */
  float &operator()(size_t i, size_t j) { return d_data[i][j]; }

  /*!
   * @brief Access operator of the Matrix element M(i,j)
   * @param i Index i
   * @param j Index j
   * @return Matrix element M(i,j)
   */
  const float &operator()(size_t i, size_t j) const { return d_data[i][j]; }

  /*!
 * @brief Computes the dot product between matrix and vector
 * @param v vector
 * @return vector Dot product
 */
  std::vector<double> dot(const std::vector<double> &v) const {

    auto r = std::vector<double>(3,0.);
    for (size_t i=0; i<3; i++)
      for (size_t j=0; j<3; j++)
        r[i] += (*this)(i,j) * v[j];

    return r;
  }

  /*!
   * @brief Computes the dot product of the matrix and the vector v
   * @param v The vector
   * @return The evaluation of the dot product
   * 
   */
  util::Point3 dot(const util::Point3 &v) {

    return {(*this)(0) * v, (*this)(1) * v, (*this)(2) * v};
  }

  /*!
 * @brief Computes the tranpose of matrix
 * @return Matrix Transpose of m
 */
  Matrix3 transpose() const {

    Matrix3 m = Matrix3(*this);

    m(0,1) = (*this)(1,0);
    m(0,2) = (*this)(2,0);

    m(1,0) = (*this)(0,1);
    m(1,2) = (*this)(2,1);

    m(2,0) = (*this)(0,2);
    m(2,1) = (*this)(1,2);

    return m;
  }

  /*!
 * @brief Computes the determinant of matrix
 *
 * @param m Matrix
 * @return det Determinant
 */
  double det() const {
    return (*this)(0,0) * ((*this)(1,1) * (*this)(2,2) - (*this)(2,1) * (*this)(1,2)) -
           (*this)(0,1) * ((*this)(1,0) * (*this)(2,2) - (*this)(2,0) * (*this)(1,2)) +
           (*this)(0,2) * ((*this)(1,0) * (*this)(2,1) - (*this)(2,0) * (*this)(1,1));
  }

  /*!
 * @brief Computes the determinant of matrix
 *
 * @param m Matrix
 * @return inv Inverse of m
 */
  Matrix3 inv() const {

    Matrix3 m = Matrix3();

    auto det_inv = 1. / this->det();

    m(0,0) = det_inv *
             ((*this)(1,1) * (*this)(2,2) - (*this)(2,1) * (*this)(1,2));
    m(0,1) = -det_inv *
             ((*this)(0,1) * (*this)(2,2) - (*this)(2,1) * (*this)(0,2));
    m(0,2) = det_inv *
             ((*this)(0,1) * (*this)(1,2) - (*this)(1,1) * (*this)(0,2));

    m(1,0) = -det_inv *
             ((*this)(1,0) * (*this)(2,2) - (*this)(2,0) * (*this)(1,2));
    m(1,1) = det_inv *
             ((*this)(0,0) * (*this)(2,2) - (*this)(2,0) * (*this)(0,2));
    m(1,2) = -det_inv *
             ((*this)(0,0) * (*this)(1,2) - (*this)(1,0) * (*this)(0,2));

    m(2,0) = det_inv *
             ((*this)(1,0) * (*this)(2,1) - (*this)(2,0) * (*this)(1,1));
    m(2,1) = -det_inv *
             ((*this)(0,0) * (*this)(2,1) - (*this)(2,0) * (*this)(0,1));
    m(2,2) = det_inv *
             ((*this)(0,0) * (*this)(1,1) - (*this)(1,0) * (*this)(0,1));

    return m;
  }
};

/*! @brief A structure to represent 3d matrices */
struct SymMatrix3 {

  /**
   * @name Data members
   *
   * 0 - xx component
   * 1 - yy component
   * 2 - zz component
   * 3 - yz component
   * 4 - xz component
   * 5 - xy component
   */
  /**@{*/

  /*! @brief data */
  float d_data[6]{};

  /** @}*/

  /*!
   * @brief Constructor
   */
  SymMatrix3() {

    d_data[0] = 0.;
    d_data[1] = 0.;
    d_data[2] = 0.;
    d_data[3] = 0.;
    d_data[4] = 0.;
    d_data[5] = 0.;
  };

  /*!
   * @brief Constructor
   *
   * @param diagonal Diagonal vector
   */
  SymMatrix3(const util::Point3& diagonal) {

    d_data[0] = diagonal.d_x;
    d_data[1] = diagonal.d_y;
    d_data[2] = diagonal.d_z;

    d_data[3] = 0.;
    d_data[4] = 0.;
    d_data[5] = 0.;
  };

  /*!
   * @brief Constructor
   *
   * @param m Matrix in vector template
   */
  SymMatrix3(const std::vector<std::vector<double>> &m) {

    d_data[0] = m[0][0];
    d_data[1] = m[1][1];
    d_data[2] = m[2][2];
    d_data[3] = 0.5 * (m[1][2] + m[2][1]);
    d_data[4] = 0.5 * (m[0][2] + m[2][0]);
    d_data[5] = 0.5 * (m[0][1] + m[1][0]);
  }

  /*!
   * @brief Constructor
   *
   * @param m Matrix in vector template
   */
  SymMatrix3(const std::vector<double> &m) {

    d_data[0] = m[0];
    d_data[1] = m[1];
    d_data[2] = m[2];
    d_data[3] = m[3];
    d_data[4] = m[4];
    d_data[5] = m[5];
  }

  /*!
   * @brief Constructor
   *
   * @param m Matrix in vector template
   */
  SymMatrix3(const Matrix3 &m) {

    d_data[0] = m(0,0);
    d_data[1] = m(1,1);
    d_data[2] = m(2,2);
    d_data[3] = 0.5 * (m(1,2) + m(2,1));
    d_data[4] = 0.5 * (m(0,2) + m(2,0));
    d_data[5] = 0.5 * (m(0,1) + m(1,0));
  }

  /*!
   * @brief Constructor
   *
   * @param m Matrix in vector template
   */
  SymMatrix3(const SymMatrix3 &m) {

    for (size_t i=0; i<6; i++)
      d_data[i] = m.d_data[i];
  }

  /*!
   * @brief Prints the information
   *
   * @param nt Number of tabs to append before printing
   * @param lvl Information level (higher means more information)
   * @return resulting string
   */
  std::string printStr(int nt = 0, int lvl = 0) const {

    std::string tabS = "";
    for (int i = 0; i < nt; i++)
      tabS += "\t";

    std::ostringstream oss;
    for (size_t i=0; i<3; i++)
      oss << tabS << "[" << (*this)(i, 0) << ", " << (*this)(i, 1) << ", "
          << (*this)(i, 2) << "]" << std::endl;
    oss << std::endl;

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
   * @brief Access operator of the i-th element
   * @param i Index i
   * @return The i-th element
   */
  Point3 operator()(size_t i) {
    return {(*this)(i, 0), (*this)(i, 1), (*this)(i, 2)};
  }

  /*!
   * @brief Access operator of the i-th element
   * @param i Index i
   * @return The i-th element
   */
  Point3 operator()(size_t i) const {
    return {(*this)(i, 0), (*this)(i, 1), (*this)(i, 2)};
  }

  /*!
   * @brief Access operator of the Matrix element M(i,j)
   * @param i Index i
   * @param j Index j
   * @return Matrix element M(i,j)
   */
  float &operator()(size_t i, size_t j) {
    return d_data[i == j ? i : 6 - i - j];
  }


  /*!
   * @brief Access operator of the Matrix element M(i,j)
   * @param i Index i
   * @param j Index j
   * @return Matrix element M(i,j)
   */
  const float &operator()(size_t i, size_t j) const {
    return d_data[i == j ? i : 6 - i - j];
  }

/*! @brief Return the i-th element
 *@param i Index
 *@return the i-th element
*/
  const float &get(size_t i) const {
    return d_data[i];
  }

/*! @brief Return the i-th element
 *@param i Index
 *@return The i-ith element
*/
  float &get(size_t i) {
    return d_data[i];
  }

  /*!warning
 * @brief Copy
 * @param m Matrix
 */
  void copy(double m[6]) const {

    for (size_t i=0; i<6; i++)
      m[i] = d_data[i];
  }

  /*!
   * @brief Computes the dot product of this matrix with another vector
   * @param v A vector
   * @return Vector Resulting vector
   */
  std::vector<double> dot(const std::vector<double> &v) const {

    auto r = std::vector<double>(3,0.);
    for (size_t i=0; i<3; i++)
      for (size_t j=0; j<3; j++)
        r[i] += (*this)(i,j) * v[j];

    return r;
  }

  /*!
   * @brief Computes the dot product of this matrix with another vector
   * @param v A vector
   * @return Vector Resulting vector
   */
  util::Point3 dot(const util::Point3 &v) {

    return {(*this)(0) * v, (*this)(1) * v, (*this)(2) * v};
  }

  /*!
 * @brief Computes the tranpose of matrix
 * @return Matrix Transpose of m
 */
  SymMatrix3 transpose() const {

    return (*this);
  }

  /*!
 * @brief Computes the determinant of matrix
 *
 * @param m Matrix
 * @return det Determinant
 */
  double det() const {
    return (*this)(0,0) * ((*this)(1,1) * (*this)(2,2) - (*this)(2,1) * (*this)(1,2)) -
           (*this)(0,1) * ((*this)(1,0) * (*this)(2,2) - (*this)(2,0) * (*this)(1,2)) +
           (*this)(0,2) * ((*this)(1,0) * (*this)(2,1) - (*this)(2,0) * (*this)(1,1));
  }

  /*!
 * @brief Computes the determinant of matrix
 *
 * @param m Matrix
 * @return inv Inverse of m
 */
  SymMatrix3 inv() const {

    SymMatrix3 m = SymMatrix3();

    auto det_inv = 1. / this->det();

    m(0,0) = det_inv *
             ((*this)(1,1) * (*this)(2,2) - (*this)(2,1) * (*this)(1,2));
    m(0,1) = -det_inv *
             ((*this)(0,1) * (*this)(2,2) - (*this)(2,1) * (*this)(0,2));
    m(0,2) = det_inv *
             ((*this)(0,1) * (*this)(1,2) - (*this)(1,1) * (*this)(0,2));

    m(1,0) = -det_inv *
             ((*this)(1,0) * (*this)(2,2) - (*this)(2,0) * (*this)(1,2));
    m(1,1) = det_inv *
             ((*this)(0,0) * (*this)(2,2) - (*this)(2,0) * (*this)(0,2));
    m(1,2) = -det_inv *
             ((*this)(0,0) * (*this)(1,2) - (*this)(1,0) * (*this)(0,2));

    m(2,0) = det_inv *
             ((*this)(1,0) * (*this)(2,1) - (*this)(2,0) * (*this)(1,1));
    m(2,1) = -det_inv *
             ((*this)(0,0) * (*this)(2,1) - (*this)(2,0) * (*this)(0,1));
    m(2,2) = det_inv *
             ((*this)(0,0) * (*this)(1,1) - (*this)(1,0) * (*this)(0,1));

    return m;
  }
};

/*!
 * @brief Checks matrix
 *
 * @param m Matrix
 * @return true If matrix is okay
 */
bool checkMatrix(const std::vector<std::vector<double>> &m);

/*!
 * @brief Computes the dot product between matrix and vector
 *
 * @param m Matrix
 * @param v vector
 * @return vector Dot produc
 */
std::vector<double> dot(const std::vector<std::vector<double>> &m, const
std::vector<double> &v);

/*!
 * @brief Computes the tranpose of matrix
 *
 * @param m Matrix
 * @return Matrix Transpose of m
 */
std::vector<std::vector<double>> transpose(const
                                           std::vector<std::vector<double>> &m);

/*!
 * @brief Computes the determinant of matrix
 *
 * @param m Matrix
 * @return det Determinant
 */
double det(const std::vector<std::vector<double>> &m);

/*!
 * @brief Computes the determinant of matrix
 *
 * @param m Matrix
 * @return inv Inverse of m
 */
std::vector<std::vector<double>>
inv(const std::vector<std::vector<double>> &m);

} // namespace util

#endif // UTIL_MATRIX_H
