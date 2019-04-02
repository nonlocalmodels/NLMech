// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef COMPARE_H
#define COMPARE_H

// Assign tolerance for comparison
#define COMPARE_EPS 1e-5

namespace util {

/*! @brief Provides comparison of floating point values  */
namespace compare {

/*!
 * @brief Compares if a is approximately equal to b
 * @param a Value a
 * @param b Value b
 * @return Result true if approximately equal else false
 */

bool approximatelyEqual(const double &a, const double &b);

/*!
 * @brief Compares if a is essentially equal to b
 * @param a Value a
 * @param b Value b
 * @return Result true if essentially equal else false
 */

bool essentiallyEqual(const double &a, const double &b);

/*!
 * @brief Compares if a > to b
 * @param a Value a
 * @param b Value b
 * @return Result true if a is definitely greater than b
 */

bool definitelyGreaterThan(const double &a, const double &b);

/*!
 * @brief Compares if a is < to b
 * @param a Value a
 * @param b Value b
 * @return Result true if a is definitely less than b
 */

bool definitelyLessThan(const double &a, const double &b);

} // namespace compare

} // namespace util

#endif
