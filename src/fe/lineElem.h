// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef FE_LINEELEM_H
#define FE_LINEELEM_H

#include "baseElem.h"     // base class BaseElem
#include "quadData.h"     // definition of QuadData

namespace fe {

/*!
 * @brief A class for quadrature related operations for linear 2-point line
 * element
 *
 * This class provides methods such as quadrature points for integration,
 * shape functions at quadrature points, and derivative of shape functions.
 * They are specific to linear quadrangle element.
 *
 * The reference line element is made of two vertex at points \f$ -1, 1 \f$.
 *
 */
class LineElem : public BaseElem {

public:
  /*!
   * @brief Constructor for line element
   * @param order Order of quadrature point approximation
   */
  explicit LineElem(size_t order);

  /*!
   * @brief Get vector of quadrature data
   *
   * Given vertices of an element, where element type is same as what is
   * specified in the constructor of Quadrature class, it returns the vector
   * of quadrature data.
   *
   * @param nodes Vector of vertices of an element
   * @return Vector of QuadData
   */
  std::vector<fe::QuadData>
  getQuadPoints(const std::vector<util::Point3> &nodes) override;

  /*!
   * @brief Returns the values of shape function at point p for line element
   *
   * Reference line element is given by vertices at v1 = -1 and v2 =1. The
   * shape function for the reference element is given by
   * \f[N_1(\xi) = \frac{(1- \xi)}{2}, \, N_2(\xi) = \frac{(1+ \xi)}{2}. \f]
   *
   * @param p Location of point
   * @return Vector of shape functions at point p
   */
  std::vector<double> getShapes(const util::Point3 &p) override;

  /*!
   * @brief Returns the values of derivative of shape function for line
   * element at point p
   *
   * For line element, derivative of shape functions are constant and are as
   * follows (for reference line element)
   *
   * \f[\frac{d N_1(\xi)}{d\xi} = \frac{-1}{2}, \, \frac{d N_2
   * (\xi)}{d\xi} = \frac{1}{2}. \f]
   *
   * @param p Point at which derivatives are to be evaluated
   * @return Vector of derivative of shape functions
   */
  std::vector<std::vector<double>> getDerShapes(const util::Point3 &p) override;

  /*!
   * @brief Maps the point on reference line element to given line element and
   * returns Jacobian
   *
   * For line element, the map from \f$ \xi \f$ to \f$ x \f$ coordinate is
   * given by
   * \f[ x = N_1(\xi) x_1 + N_2(\xi) x_2,\f]
   * where \f$ N_1, N_2\f$ are described in LineElem::getShapes() and
   * \f$ x_1, x_2\f$ are vertices of given line element.
   *
   * Jacobian of this map is a matrix
   * \f[ J = \frac{dx}{d\xi} = \frac{dN_1(\xi)}{d\xi} x_1 + \frac{dN_2(\xi)
   * }{d\xi} x_2 = \frac{(x_2 - x_1)}{2}. \f]
   *
   * @param p Given point in reference triangle which is to be mapped
   * @param shapes Vector shape functions evaluated at the point p
   * @param der_shapes Vector of derivative of shape functions at point p
   * @param nodes Coordinates of vertices of a given element
   * @return det(J) Determinant of the Jacobian
   */
  double mapRefElemToElem(util::Point3 &p, const std::vector<double> &shapes,
                          const std::vector<std::vector<double>> &der_shapes,
                          const std::vector<util::Point3> &nodes) override;

private:
  /*!
   * @name Internal methods for triangle element
   */

  /*!
   * @brief Compute the quadrature points for triangle element
   */
  void init() override;

  /** @}*/
};

} // namespace fe

#endif // FE_LINEELEM_H
