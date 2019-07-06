// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef UTIL_FAST_METHODS_H
#define UTIL_FAST_METHODS_H

#include "point.h"           // definition of Point3
#include <vector>

namespace util {

/*!
 * @brief Provides fast methods to add/subtract list of data, to find
 * maximum/minimum from list of data
 */
namespace methods {

/*!
 * @brief Returns the sum of data
 * @param data List of real numbers
 * @return sum Sum of the numbers
 */
double add(const std::vector<double> &data);

/*!
 * @brief Returns the maximum from list of data
 * @param data List of real numbers
 * @param i Pointer to store the id where maximum occurs
 * @return max Maximum value
 */
double max(const std::vector<double> &data, size_t *i = nullptr);

/*!
 * @brief Returns the minimum from list of data
 * @param data List of real numbers
 * @param i Pointer to store the id where minimum occurs
 * @return min Minimum value
 */
double min(const std::vector<double> &data, size_t *i = nullptr);

/*!
 * @brief Returns the sum of data
 * @param data List of real numbers
 * @return sum Sum of the numbers
 */
float add(const std::vector<float> &data);

/*!
 * @brief Returns the maximum from list of data
 * @param data List of real numbers
 * @param i Pointer to store the id where maximum occurs
 * @return max Maximum value
 */
float max(const std::vector<float> &data, size_t *i = nullptr);

/*!
 * @brief Returns the minimum from list of data
 * @param data List of real numbers
 * @param i Pointer to store the id where minimum occurs
 * @return min Minimum value
 */
float min(const std::vector<float> &data, size_t *i = nullptr);

/*!
 * @brief Returns the maximum length of point from list of points
 * @param data List of points
 * @return max Maximum length of point
 */
util::Point3 maxLength(const std::vector<util::Point3> &data);

} // namespace methods

} // namespace util

#endif // UTIL_FAST_METHODS_H
