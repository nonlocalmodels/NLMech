////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef GEOM_DAMPINGGEOM_H
#define GEOM_DAMPINGGEOM_H

#include "util/point.h"           // definition of Point3
#include <string>
#include <vector>

// forward declaration of interior flags deck
namespace inp {
struct AbsorbingCondDeck;
}

namespace fe {
class Mesh;
}

namespace geometry {

/*!
 * @brief An abstraction class to process geometry for damping force
 * calculation
 */
class DampingGeom {
public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   * @param mesh Fe mesh object
   */
  DampingGeom(inp::AbsorbingCondDeck *deck,
               const fe::Mesh *mesh);

  bool isDampingActive();
  bool isDampingActive() const;

  bool isViscousDamping();
  bool isViscousDamping() const;

  /*!
   * @brief Get nodal volume of node i
   * @param i Id of the node
   * @return vol Volume
   */
  double getCoefficient(const size_t &i);
  double getCoefficient(const size_t &i) const ;

  const std::vector<double> *getCoefficientDataP() const ;
  const std::vector<double> *getCoefficientDataP();


protected:
  /*!
   * @brief Compute damping coefficients at nodal coordinates
   *
   * @param coefficients Coefficients at nodal coordinates
   */
  void computeDampingCoefficient(const fe::Mesh *mesh);

  size_t d_dim;

  /*! @brief input deck for absorbing condition */
  inp::AbsorbingCondDeck *d_absorbingDeck_p;

  std::vector<double> d_coefficients;
};

} // namespace geometry

#endif // GEOM_DAMPINGGEOM_H
