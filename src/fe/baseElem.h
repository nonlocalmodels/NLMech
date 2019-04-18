// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef FE_BASEELEM_H
#define FE_BASEELEM_H

#include "util/point.h" // definition of Point3
#include "quadData.h"      // definition of QuadData

namespace fe {

/*!
 * @brief A base class to compute and store quadrature data
 *
 * This class provides methods such as quadrature points for integration,
 * shape functions at quadrature points, and derivative of shape functions.
 *
 * @sa TriElem, QuadElem
 */
class BaseElem {

public:
  /*!
   * @brief Constructor
   *
   * @param order Order of quadrature point approximation
   * @param element_type Type of element in the mesh
   */
  BaseElem(size_t order, size_t element_type);

  /*!
   * @brief Get element type
   * @return type Type of element
   */
  size_t getElemType() { return d_elemType; }

  /*!
   * @brief Get order of quadrature approximation
   * @return order Order of approximation
   */
  size_t getQuadOrder() { return d_quadOrder; }

  /*!
   * @brief Get number of quadrature points in the data
   * @return N Number of quadrature points
   */
  size_t getNumQuadPoints() { return d_numQuadPts; }

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
  virtual std::vector<fe::QuadData>
  getQuadPoints(const std::vector<util::Point3> &nodes) = 0;

  /*!
   * @brief Returns the values of shape function at point p
   *
   * @param p Location of point
   * @return Vector of shape functions at point p
   */
  virtual std::vector<double> getShapes(const util::Point3 &p) = 0;

  /*!
   * @brief Returns the values of derivative of shape function
   *
   * @param p Location of point
   * @return Vector of derivative of shape functions
   */
  virtual std::vector<std::vector<double>>
  getDerShapes(const util::Point3 &p) = 0;

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
  virtual double
  mapRefElemToElem(util::Point3 &p, const std::vector<double> &shapes,
                   const std::vector<std::vector<double>> &der_shapes,
                   const std::vector<util::Point3> &nodes) = 0;

protected:
  /*!
   * @brief Compute the quadrature points
   *
   * This must be implemented by inheriting classes.
   */
  virtual void init() = 0;

  /*! @brief Order of quadrature point integration approximation */
  size_t d_quadOrder;

  /*! @brief Number of quadrature points for order d_quadOrder */
  size_t d_numQuadPts;

  /*! @brief Element type */
  size_t d_elemType;

  /*! @brief Quadrature data collection */
  std::vector<fe::QuadData> d_quads;
};

} // namespace fe

#endif // FE_BASEELEM_H
