////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILS_UTIL_IO_H
#define UTILS_UTIL_IO_H

#include "point.h"
#include <vector>
#include <iostream>

namespace util {

/*! @brief Provides geometrical methods such as point inside rectangle */
namespace io {


/*! @brief Generate a string contaning nt tabs
* @param nt Number of tabs
* @return String containing nt tabs
*/
inline std::string getTabS(int nt) {
  std::string tabS = "";
  for (int i = 0; i < nt; i++)
    tabS += "\t";

  return tabS;
};

/*!
 * @brief Concatenates the elements of the vector separated with a comma
 * @param list The vector with the elements
 * @param nt Number of tabs prepend to the concatenated vector
 * @return Resulting string
 */
template <class T> inline std::string printStr(const std::vector<T> &list,
    int nt =
    0) {

  auto tabS = getTabS(nt);
  std::ostringstream oss;
  oss << tabS;
  size_t i = 0;
  for (auto l : list) {
    oss << l;
    i++;
    if (i != list.size())
      oss << ", ";
  }

  return oss.str();
};

/*!
 * @brief Concatenates the elements of the vector separated with a comma
 * @param list The vector with the elements
 * @param nt NUmber of tabs prepend to the concatenated vector
 * @return Resulting string
 */
template <> inline std::string printStr(const std::vector<util::Point3> &list,
                                               int nt) {

  auto tabS = getTabS(nt);
  std::ostringstream oss;
  oss << tabS;
  size_t i = 0;
  for (auto l : list) {
    oss << "(" << l[0] << ", " << l[1] << ", " << l[2] << ")";
    i++;
    if (i != list.size())
      oss << ", ";
  }

  return oss.str();
};

/*!
 * @brief Concatenates the elements of the vector separated with a comma and
 * prints them to the standard output stream
 * @param list The vector with the elements
 * @param nt Number of tabs prepend to the concatenated vector
 */
template <class T> inline void print(const std::vector<T> &list, int nt = 0) {

  std::cout << printStr(list, nt);
};

/*!
 * @brief Concatenates the elements of the vector separated with a comma
 * @param list The vector with the elements
 * @param nt Number of tabs prepend to the concatenated vector
 * @return Resulting string
 */
template <class T>
inline std::string printStr(const std::vector<std::vector<T>> &list,
                            int nt = 0) {

  auto tabS = getTabS(nt);
  std::ostringstream oss;
  oss << tabS;
  size_t i = 0;
  for (auto l : list) {
    oss << "(";
    for (size_t k = 0; k < l.size(); k++) {
      oss << l[k];
      if (k < l.size() - 1)
        oss << ", ";
    }
    oss << ")";

    i++;
    if (i != list.size())
      oss << ", ";
  }

  return oss.str();
}

/*!
 * @brief Concatenates the elements of the vector separated with a comma
 * @param list The vector with the elements
 * @param nt Number of tabs prepend to the concatenated vector
 */
template <class T> inline void print(const std::vector<std::vector<T>> &list, int nt = 0) {

  std::cout << printStr(list, nt);
};

/*!
 * @brief Prints the corner points of the to the output stream
 * @param box Corner point of the box
 * @param nt Number of tabs at the beginning of the string
 * @return Resulting string
 */
inline std::string printBoxStr(const std::pair<util::Point3, util::Point3>
    &box, int nt =
0) {
  auto tabS = getTabS(nt);
  std::ostringstream oss;
  oss << tabS << "Corner point 1 = " << box.first.printStr(nt, 0) << std::endl;
  oss << tabS << "Corner point 2 = " << box.second.printStr(nt, 0) << std::endl;

  return oss.str();
};

/*!
 * @brief Prints the corner points of the to the output stream
 * @param box Corner point of the box
 * @param nt Number of tabs at the beginning of the string
 */
inline void printBox(const std::pair<util::Point3, util::Point3> &box, int nt
= 0) {
  std::cout << printBoxStr(box, nt);
};

/*!
 * @brief Prints the corner points of the output stream
 * @param box Corner point of the box
 * @param nt Number of tabs at the beginning of the string
 * @return Resulting string
 */
inline std::string printBoxStr(const std::pair<std::vector<double>, std::vector<double>>
                               &box, int nt =
0) {
  auto tabS = getTabS(nt);
  std::ostringstream oss;
  oss << tabS << "Corner point 1 = (" << printStr<double>(box.first, 0)
              << ")" << std::endl;
  oss << tabS << "Corner point 2 = (" << printStr<double>(box.second, 0)
              << ")" << std::endl;

  return oss.str();
};

/*!
 * @brief Prints the corner points of the to the output stream
 * @param box Corner point of the box
 * @param nt Number of tabs at the beginning of the string
 */
inline void printBox(const std::pair<std::vector<double>, std::vector<double>> &box, int nt
= 0) {
  std::cout << printBoxStr(box, nt);
};
} // namespace io

} // namespace util

#endif // UTILS_UTIL_IO_H
