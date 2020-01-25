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

  /*! @brief the xx component */
  float d_xx;

  /*! @brief the xy component */
  float d_xy;

  /*! @brief the xz component */
  float d_xz;

  /*! @brief the yx component */
  float d_yx;

  /**< @brief the yy component */
  float d_yy;

  /*! @brief the yz component */
  float d_yz;

  /*! @brief the zx component */
  float d_zx;

  /*! @brief the zy component */
  float d_zy;

  /*! @brief the zz component */
  float d_zz;

  /** @}*/

  /*!
   * @brief Constructor
   */
  Matrix3()
      : d_xx(0.), d_xy(0.), d_xz(0.), d_yx(0.), d_yy(0.), d_yz(0.), d_zx(0.),
        d_zy(0.), d_zz(0.){};
};

/*! @brief A structure to represent 3d symmetric matrices */
struct SymMatrix3 {

  /**
   * @name Data members
   */
  /**@{*/

  /*! @brief the xx component */
  float d_xx;

  /*! @brief the yy component */
  float d_yy;

  /*! @brief the zz component */
  float d_zz;

  /*! @brief the xy (yx) component */
  float d_xy;

  /*! @brief the xz (zx) component */
  float d_xz;

  /*! @brief the yz (zy) component */
  float d_yz;

  /** @}*/

  /*!
   * @brief Constructor
   */
  SymMatrix3() : d_xx(0.), d_yy(0.), d_zz(0.), d_xy(0.), d_xz(0.), d_yz(0.){};

  /*!
   * @brief Constructor
   */
  SymMatrix3(const float &xx, const float &yy, const float &zz, const float &xy,
             const float &xz, const float &yz)
      : d_xx(xx), d_yy(yy), d_zz(zz), d_xy(xy), d_xz(xz), d_yz(yz){};
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
