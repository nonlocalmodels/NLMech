////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTIL_COMPARE_H
#define UTIL_COMPARE_H

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

#endif // UTIL_COMPARE_H
