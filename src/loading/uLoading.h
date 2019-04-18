// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

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
};

} // namespace loading

#endif // LOADING_ULOADING_H
