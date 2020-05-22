////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
// ////////////////////////////////////////////////////////////////////////////////

#include "utilIO.h"
//
// std::string util::io::getTabS(int nt) {
//  std::string tabS = "";
//  for (int i = 0; i < nt; i++)
//    tabS += "\t";
//
//  return tabS;
//};
//
// template <class T> std::string util::io::printStr(const std::vector<T> &list,
// int nt = 0) {
//
//  auto tabS = getTabS(nt);
//  std::ostringstream oss;
//  oss << tabS;
//  size_t i = 0;
//  for (auto l : list) {
//    oss << l;
//    i++;
//    if (i != list.size())
//      oss << ", ";
//  }
//
//  return oss.str();
//};
//
// template <class T> void util::io::print(const std::vector<T> &list, int nt =
// 0) {
//
//  std::cout << printStr(list, nt);
//};
//
///*!
// * @brief Prints box to std::cout
// */
//
// std::string util::io::printBoxStr(const std::pair<util::Point3, util::Point3>
// &box, int nt = 0) {
//  auto tabS = getTabS(nt);
//  std::ostringstream oss;
//  oss << tabS << "Corner point 1 = " << box.first.printStr(nt, 0) <<
//  std::endl; oss << tabS << "Corner point 2 = " << box.second.printStr(nt, 0)
//  << std::endl;
//
//  return oss.str();
//};
//
// void util::io::printBox(const std::pair<util::Point3, util::Point3> &box, int
// nt =
//    0) {
//  std::cout << printBoxStr(box, nt);
//}