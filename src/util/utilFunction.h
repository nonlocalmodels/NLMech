////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTIL_FUNCTION_H
#define UTIL_FUNCTION_H

#include "point.h"              // definition of Point3
#include <vector>

namespace util {

/*! @brief Provides geometrical methods such as point inside rectangle */
namespace function {

/*!
 * @brief Computes hat function at given point
 *
 * Hat function:
 *      f ^
 *        |
 *        |
 *     1  o
 *        |           /|\
 *        |         /  |  \
 *        |       /    |    \
 *        |     /      |      \
 *        |   /        |        \
 *        | /          |          \
 *        o____________o____________o______\ x
 *                                         /
 *      x_min                      x_max
 *
 * @param x Point in real line
 * @param x_min Left side point in real line
 * @param x_max Right side point in real line
 * @return value Evaluation of hat function at x
 */
double hatFunction(const double &x, const double &x_min, const double &x_max);

/*!
 * @brief Computes hat function at given point
 *
 * This version does not test if point x is in valid interval.
 *
 * Hat function:
 *      f ^
 *        |
 *        |
 *     1  o
 *        |           /|\
 *        |         /  |  \
 *        |       /    |    \
 *        |     /      |      \
 *        |   /        |        \
 *        | /          |          \
 *        o____________o____________o______\ x
 *                                         /
 *      x_min                      x_max
 *
 * @param x Point in real line
 * @param x_min Left side point in real line
 * @param x_max Right side point in real line
 * @return value Evaluation of hat function at x
 */
double hatFunctionQuick(const double &x, const double &x_min,
                        const double &x_max);

/*!
 * @brief Compute linear step function
 *
 * Step function:
 *
 * f ^
 *   |             __________
 *   |            /
 *   |           /
 *   |   _______/
 *   |  /
 *   | /
 *   |/_________________________ t
 *       x1   x1+x2
 *
 *  - Linear (with slope 1) in [0,l1), constant in [l1,l1+l2)
 *  - Periodic with periodicity l1+l2
 *
 * @param x  Point in real line
 * @param x1 Point such that function is linear with slope 1 in [0, x1)
 * @param x2 Point such that function is constant in [x1, x1 + x2)
 * @return value Evaluation of step function at x
 */
double linearStepFunc(const double &x, const double &x1, const double &x2);

/*!
 * @brief Compute derivative of linear step function
 *
 * @param x  Point in real line
 * @param x1 Point such that function is linear with slope 1 in [0, x1)
 * @param x2 Point such that function is constant in [x1, x1 + x2)
 * @return value Derivative of step function at x
 */
double derLinearStepFunc(const double &x, const double &x1, const double &x2);

/*!
 * @brief Compute Gaussian function in 1-d
 *
 * Gaussian (1-d) function: \f$ f(r) = a \exp(-\frac{r^2}{\beta}). \f$
 *
 * Here \f$ a\f$ is the amplitude and \f$ \beta \f$ is the exponential factor.
 *
 * @param r Distance from origin
 * @param a Amplitude
 * @param beta Factor in exponential function
 * @return value Component of Gaussian 1-d function
 */
double gaussian(const double &r, const double &a, const double &beta);

/*!
 * @brief Compute Gaussian function in 2-d
 *
 * Gaussian (2-d) function:
 * \f[ f(x,y) = (f_1(x,y), f_2(x,y)), \f]
 * where
 * \f[ f_1(x,y) =  a \exp(-\frac{(x-x_c)^2 + (y-y_c)^2}{\beta}) d_1,
 * \quad f_1(x,y) =  a \exp(-\frac{(x-x_c)^2 + (y-y_c)^2}{\beta}) d_2.
 * \f]
 * Here \f$ (x_c,y_c) \f$ is the center of the pulse, \f$ a\f$ is the
 * amplitude, \f$ \beta \f$ is the exponential factor, and \f$ (d_1,d_2)\f$
 * is the direction of the pulse.
 *
 * @param x  Coordinates of point
 * @param params List of parameters
 * @param dof Component of Gaussian function
 * @return value Component of Gaussian 2-d vector function along dof
 */
double gaussian2d(const util::Point3 &x, const size_t &dof,
                  const std::vector<double> &params);

/*!
 * @brief Compute sum of two Gaussian function in 2-d
 *
 * Double Gaussian (2-d) function:
 * \f[ f(x,y) = (f_1(x,y), f_2(x,y)) + (g_1(x,y), g_2(x,y)), \f]
 * where \f$ (f_1,f_2)\f$ and \f$(g_1, g_2)\f$ are two Gaussian 2-d function
 * as described in Gaussian2d() with different values of \f$ (x_c, y_c), a,
 * (d_1, d_2)\f$.
 *
 * @param x  Coordinates of point
 * @param params List of parameters
 * @param dof Component of Gaussian function
 * @return value Component of Gaussian 2-d vector function along dof
 */
double doubleGaussian2d(const util::Point3 &x, const size_t &dof,
                        const std::vector<double> &params);

/*!
 * @brief Get vector with plus or minus 1 depending on the sign of component
 * of another vector
 *
 * @param v Vector from which a new vector is created
 * @return vec Vector
 */
util::Point3 signVector(const util::Point3 &v);

/*!
 * @brief Get determinant of 3-d matrix
 *
 * @param rows Three rows of matrix
 * @return det Determinant
 */
double getDeterminant(const std::vector<util::Point3> &rows);

} // namespace function

} // namespace util

#endif // UTIL_FUNCTION_H
