////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef LOADING_ULOADING_H
#define LOADING_ULOADING_H

#include "loading.h"        // base class Loading
#include "util/point.h"     // definition of Point3

// forward declaration
namespace fe {
class Mesh;
}

namespace loading {

/*!
 * @brief A class to apply displacement boundary condition
 */
class ULoading : public Loading {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   * @param mesh Mesh object
   */
  ULoading(inp::LoadingDeck *deck, fe::Mesh *mesh);

  /*!
   * @brief Applies displacement boundary condition
   * @param time Current time
   * @param u Vector nodal displacements
   * @param v Vector nodal velocities
   * @param mesh Mesh object
   */
  void apply(const double &time, std::vector<util::Point3> *u,
             std::vector<util::Point3> *v, fe::Mesh *mesh);

  /*!
   * @brief Returns the string containing information about the instance of
   * the object
   *
   * @param nt Number of tabs to append before each line of string
   * @param lvl Level of information sought (higher level means more
   * information)
   * @return string String containing information about this object
   * */
  std::string printStr(int nt = 0, int lvl = 0) const;

  /*!
   * @brief Prints the information about the instance of the object
   *
   * @param nt Number of tabs to append before each line of string
   * @param lvl Level of information sought (higher level means more
   * information)
   * */
  void print(int nt = 0, int lvl = 0) const { std::cout << printStr(nt, lvl); };
};

} // namespace loading

#endif // LOADING_ULOADING_H
