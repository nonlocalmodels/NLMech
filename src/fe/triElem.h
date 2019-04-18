// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef FE_TRIELEM_H
#define FE_TRIELEM_H

#include "baseElem.h"       // base class BaseElem
#include "quadData.h"       // definition of QuadData

namespace fe {

/*!
 * @brief A class for quadrature related operations for linear triangle element
 *
 * This class provides methods such as quadrature points for integration,
 * shape functions at quadrature points, and derivative of shape functions.
 * They are specific to linear triangle element.
 *
 * The reference triangle element is made of three vertex at point \f$ (0,0),
 * \, (1,0), \, (0,1) \f$.
 *
 */
class TriElem : public BaseElem {

public:
  /*!
   * @brief Constructor for triangle element
   * @param order Order of quadrature point approximation
   */
  explicit TriElem(size_t order);

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
   * @brief Returns the values of shape function at point p for triangle element
   *
   * Triangle assumed is a reference triangle with vertices v1 = (0,0), v2 =
   * (1,0), v3 = (0,1). The shape function for the reference triangle is
   * given by
   * \f[N_1(\xi, \eta) = 1- \xi - \eta, \, N_2(\xi, \eta) = \xi, \, N_3(\xi,
   * \eta) = \eta. \f]
   *
   * @param p Location of point
   * @return Vector of shape functions at point p
   */
  std::vector<double> getShapes(const util::Point3 &p) override;

  /*!
   * @brief Returns the values of derivative of shape function for
   * triangle element
   *
   * For linear triangle element, derivative of shape functions are constant
   * and are as follows (for reference triangle)
   *
   * \f[\frac{d N_1(\xi, \eta)}{d\xi} = -1, \, \frac{d N_1(\xi, \eta)}{d\eta}
   * = -1, \f]
   * \f[\frac{d N_2(\xi, \eta)}{d\xi} = 1, \, \frac{d N_2(\xi, \eta)}{d\eta}
   * = 0, \f]
   * \f[\frac{d N_3(\xi, \eta)}{d\xi} = 0, \, \frac{d N_3(\xi, \eta)}{d\eta}
   * = 1. \f]
   *
   * @param p Location of point
   * @return Vector of derivative of shape functions
   */
  std::vector<std::vector<double>> getDerShapes(const util::Point3 &p) override;

  /*!
   * @brief Maps the point on reference triangle to given triangle and
   * returns determinant of Jacobian
   *
   * For linear triangle element, the map from \f$ (\xi, \eta) \f$ to \f$ (x,y)
   * \f$ coordinate is given by
   * \f[ x = N_1(\xi, \eta) x_1 + N_2(\xi, \eta) x_2 + N_3(\xi, \eta) x_3,\f]
   * \f[ y = N_1(\xi, \eta) y_1 + N_2(\xi, \eta) y_2 + N_3(\xi, \eta) y_3,\f]
   * where \f$ N_1, N_2, N_3 \f$ are described in TriElem::getShapes() and \f$
   * (x_1,
   * y_1), (x_2,y_2), (x_3,y_3) \f$ are vertices of given triangle.
   *
   * Jacobian of this map is a matrix
   * \f[ J = [\frac{dx}{d\xi},\frac{dy}{d\xi}; \frac{dx}{d\eta},
   * \frac{dy}{d\eta}] \f]
   * and determinant of Jacobian is
   * \f[ det(J) = \frac{dx}{d\xi} \times \frac{dy}{d\eta} -
   * \frac{dy}{d\xi}\times \frac{dx}{d\eta}. \f]
   *
   * It can be easily checked than for linear triangle elements, we simply have
   * \f[ \frac{dx}{d\xi} = x_2 - x_1, \frac{dx}{d\eta} = x_3 - x_1, \f]
   * \f[ \frac{dy}{d\xi} = y_2 - y_1, \frac{dy}{d\eta} = y_3 - y_1. \f]
   * Therefore,
   * \f[ det(J) = (x_2 - x_1)(y_3 - y_1) - (y_2 - y_1)(x_3 - x_1).\f]
   *
   * @param p Given point in reference triangle which is to be mapped
   * @param shapes Vector shape functions evaluated at the point p
   * @param der_shapes Vector of derivative of shape functions at point p
   * @param nodes Coordinates of vertices of a given element   *
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

#endif // FE_TRIELEM_H
