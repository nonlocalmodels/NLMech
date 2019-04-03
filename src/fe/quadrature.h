// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef FE_QUADRATURE_H
#define FE_QUADRATURE_H

#include <vector>

#include "../util/point.h" // definition of point
#include "quadData.h"      // definition of QuadData

namespace fe {

/*!
 * @brief A class for quadrature related operations
 *
 * In this class we store the global quadrature approximation information
 * such as order of quadrature approximation, number of quadrature point per
 * element, quadrature data for reference element.
 *
 * This class also provides method for obtaining quadrature points and
 * quadrature weights.
 *
 * @note Currently, only triangle and quadrangle elements are supported.
 */
template <class T> class Quadrature {

public:
  /*!
   * @brief Constructor
   * @param order Order of quadrature approximation
   */
  explicit Quadrature(size_t order);

  /*!
   * @brief Get element type
   * @return type Type of element
   */
  size_t getElemType();

  /*!
   * @brief Get order of quadrature approximation
   * @return order Order of approximation
   */
  size_t getQuadOrder();

  /*!
   * @brief Get number of quadrature points in the data
   * @return N Number of quadrature points
   */
  size_t getNumQuadPoints();

  /*!
   * @brief Returns the values of shape function at point p
   *
   * @param p Location of point
   * @return Vector of shape functions at point p
   */
  std::vector<double> getShapes(const util::Point3 &p);

  /*!
   * @brief Returns the values of derivative of shape function
   *
   * @param p Location of point
   * @return Vector of derivative of shape functions
   */
  std::vector<std::vector<double>>
  getDerShapes(const util::Point3 &p);

  /*!
   * @brief Maps the point on reference element to given element and
   * returns the determinant of Jacobian
   *
   * @param p Given point in reference element which is to be mapped
   * @param shapes Vector shape functions evaluated at the point p
   * @param der_shapes Vector of derivative of shape functions at point p
   * @param nodes Coordinates of vertices of a given element
   * @return det(J) Determinant of the Jacobian
   */
  double
  mapRefElemToElem(util::Point3 &p, const std::vector<double> &shapes,
                   const std::vector<std::vector<double>> &der_shapes,
                   const std::vector<util::Point3> &nodes);

  /*!
   * @brief Get vector of quadrature data
   *
   * Given vertices of an element, where element type is same as what is
   * specified in the constructor of Quadrature class, it returns the vector
   * of quadrature data.
   *
   * @param nodes Vector of vertices of an element
   * @return Vector of QuadData
   *
   * @note This depends on the element type that is already set in the
   * constructor of Quadrature class.
   *
   * @note Caller needs to ensure that order does not go higher than 5 as at
   * present only upto fifth order quadrature points are implemented.
   */
  std::vector<fe::QuadData>
  getQuadPoints(const std::vector<util::Point3> &nodes);

private:
  /*! @brief Element which computes quadrature point and shape functions */
  T *d_element_p;
};

} // namespace fe

#endif // FE_QUADRATURE_H
