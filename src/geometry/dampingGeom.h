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


  /*! Checks if damping is active
  * @return If damping is active
  */
  bool isDampingActive();

    /*! Checks if damping is active
  * @return If damping is active
  */
  bool isDampingActive() const;

    /*! Checks if viscous damping is active
  * @return If viscous damping is active
  */
  bool isViscousDamping();

     /*! Checks if viscous damping is active
  * @return If viscous damping is active
  */
  bool isViscousDamping() const;

  /*!
   * @brief Get nodal volume of node i
   * @param i Id of the node
   * @return vol Volume
   */
  double getCoefficient(const size_t &i);

    /*!
   * @brief Get nodal volume of node i
   * @param i Id of the node
   * @return vol Volume
   */
  double getCoefficient(const size_t &i) const ;


  /*!
   * @brief Get pointer to the coefficient data
   * @return The coefficient data
   */
  const std::vector<double> *getCoefficientDataP() const ;

    /*!
   * @brief Get pointer to the coefficient data
   * @return The coefficient data
   */
  const std::vector<double> *getCoefficientDataP();


protected:
  /*!
   * @brief Compute damping coefficients at nodal coordinates
   *
   * @param mesh The mesh
   */
  void computeDampingCoefficient(const fe::Mesh *mesh);

  /*! @brief Dimension */
  size_t d_dim;

  /*! @brief input deck for absorbing condition */
  inp::AbsorbingCondDeck *d_absorbingDeck_p;

  /*! @brief Coefficients */
  std::vector<double> d_coefficients;
};

} // namespace geometry

#endif // GEOM_DAMPINGGEOM_H
