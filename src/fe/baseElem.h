////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef FE_BASEELEM_H
#define FE_BASEELEM_H

#include "quadData.h"   // definition of QuadData
#include "util/point.h" // definition of Point3

namespace fe {

/*!
 * @brief A base class to compute and store quadrature data
 *
 * This class provides methods such as quadrature points for integration,
 * shape functions at quadrature points, and derivative of shape functions.
 *
 * Intrinsic each type of element, there are data types. These are as follows
 *
 * 1. A simple element is considered as reference element \f$ T^0 \f$. For
 * example, in triangle element a triangle formed by points \f$(0,0), (1,0),
 * (0,1) \f$ is considered as reference element.
 *
 * 2. Points in reference element \f$ T^0 \f$ is described by \f$ (\xi, \eta,
 * \zeta)\f$ (in 3-d), \f$ (\xi, \eta) \f$ (in 2-d), and \f$ \xi\f$ (in 1-d).
 * Points in any other element \f$ T \f$ is described by \f$ (x,y,z)\f$ (in
 * 3-d), \f$ (x,y) \f$ (in 2-d), and \f$ x\f$ (in 1-d). In what follows, we
 * will present theory for general 3-d element. For 2-d, reader should ignore
 * coordinates \f$\zeta \f$ and \f$ z\f$, and similarly for 1-d reader should
 * ignore \f$ \eta, \zeta \f$ and \f$ y,z\f$.
 *
 * 3. Associated to vertices of element \f$ T^0 \f$ we have shape function
 * \f$ N^0_1, N^0_2, ..., N^0_n \f$, where \f$ n \f$ is the number of
 * vertices in the element. Shape functions \f$ N^0_i \f$ are functions of
 * point \f$ (\xi, \eta, \zeta) \in T^0 \f$. For each type of fem element, \f$ n
 * \f$ is fixed. E.g., for TriElem, \f$ n = 3\f$, QuadElem \f$ n= 4\f$, LineElem
 * \f$ n = 2 \f$.
 *
 * 4. For any element \f$ T \f$ shape functions are denoted as \f$ N_1, N_2,
 * ..., N_n \f$ and are functions of point \f$ (x,y,z) \in T \f$.
 *
 * 5. Element \f$ T \f$ is described by vertices \f$ v^1, v^2, ..., v^n \f$.
 *
 * 6. A map \f$\Phi : T^0 \to T \f$ , where \f$ T \f$ is a given
 * element formed by vertices \f$ v^1, v^2, ..., v^n \f$, and \f$ T^0 \f$ is a
 * reference element, is defined as follows:
 * \f[x(\xi, \eta, \zeta) = \sum_{i=1}^n N^0_i(\xi, \eta, \zeta) v^i_x, \quad
 * y(\xi, \eta, \zeta) = \sum_{i=1}^n N^0_i(\xi, \eta, \zeta) v^i_y, \quad z
 * (\xi, \eta, \zeta) = \sum_{i=1}^n N^0_i(\xi, \eta, \zeta) v^i_z \f]
 * where \f$ v^i_x, v^i_y, v^i_z \f$ are the x, y, and z component of point
 * \f$ v^i\f$.
 *
 * 7. Jacobian of map \f$ \Phi : T^0 \to T \f$ is given by
 * \f[ J = \left[ {
 * \begin{array}{ccc}
 * \frac{dx}{d\xi} &\frac{dy}{d\xi} & \frac{dz}{d\xi} \\
 * \frac{dx}{d\eta} & \frac{dy}{d\eta} & \frac{dz}{d\eta} \\
 * \frac{dx}{d\zeta} & \frac{dy}{d\zeta} & \frac{dz}{d\zeta} \\
 * \end{array}
 * } \right]. \f]
 * For 2-d element it is
 * \f[ J = \left[ {
 * \begin{array}{cc}
 * \frac{dx}{d\xi} &\frac{dy}{d\xi} \\
 * \frac{dx}{d\eta} & \frac{dy}{d\eta} \\
 * \end{array}
 * } \right]. \f]
 * For 1-d element it is
 * \f[ J = \frac{dx}{d\xi}. \f]
 * For 2-d and 3-d element, determinant of Jacobian is also used in many
 * calculation.
 *
 * @sa LineElem, TriElem, QuadElem
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
   * @brief Returns the size of element (length in 1-d, area in 2-d, volume
   * in 3-d element)
   *
   * @param nodes Vertices of element
   * @return Vector of shape functions at point p
   */
  virtual double elemSize(const std::vector<util::Point3> &nodes) = 0;

  /*!
   * @brief Returns the values of shape function at point p
   *
   * @param p Location of point
   * @param nodes Vertices of element
   * @return Vector of shape functions at point p
   */
  virtual std::vector<double>
  getShapes(const util::Point3 &p, const std::vector<util::Point3> &nodes);

  /*!
   * @brief Returns the values of derivative of shape function at point p
   *
   * @param p Location of point
   * @param nodes Vertices of element
   * @return Vector of derivative of shape functions
   */
  virtual std::vector<std::vector<double>>
  getDerShapes(const util::Point3 &p,
               const std::vector<util::Point3> &nodes);

  /*!
   * @brief Get vector of quadrature data
   *
   * Given element vertices, this method returns the quadrature point, where
   * order of quadrature approximation and element type is set in the
   * constructor. For each quadrature point, data includes
   * - quad point
   * - quad weight
   * - shape function evaluated at quad point
   * - derivative of shape function evaluated at quad point
   * - Jacobian matrix
   * - determinant of the Jacobian
   *
   * @param nodes Vector of vertices of an element
   * @return Vector of QuadData
   */
  virtual std::vector<fe::QuadData>
  getQuadDatas(const std::vector<util::Point3> &nodes) = 0;

  /*!
   * @brief Get vector of quadrature data
   *
   * Given element vertices, this method returns the quadrature point, where
   * order of quadrature approximation and element type is set in the
   * constructor. For each quadrature point, data includes
   * - quad point
   * - quad weight
   * - shape function evaluated at quad point
   *
   * This function is lite version of BaseElem::getQuadDatas.
   *
   * @param nodes Vector of vertices of an element
   * @return Vector of QuadData
   */
  virtual std::vector<fe::QuadData>
  getQuadPoints(const std::vector<util::Point3> &nodes) = 0;

protected:
  /*!
   * @brief Returns the values of shape function at point p (reference element)
   *
   * @param p Location of point
   * @return Vector of shape functions at point p
   */
  virtual std::vector<double> getShapes(const util::Point3 &p) = 0;

  /*!
   * @brief Returns the values of derivative of shape function at point p
   * (reference element)
   *
   * @param p Location of point
   * @return Vector of derivative of shape functions
   */
  virtual std::vector<std::vector<double>>
  getDerShapes(const util::Point3 &p) = 0;

  /*!
   * @brief Maps point p in given element to the reference element and
   * returns the mapped point
   *
   * @param p Location of point
   * @param nodes Vertices of element
   * @return Vector of shape functions at point p
   */
  virtual util::Point3
  mapPointToRefElem(const util::Point3 &p,
                    const std::vector<util::Point3> &nodes);

  /*!
   * @brief Computes Jacobian of map from reference element \f$ T^0 \f$ to
   * given element \f$ T \f$
   *
   * @param p Location of point in reference element
   * @param nodes Vertices of element
   * @param J Matrix to store the Jacobian
   * @return det(J) Determinant of the Jacobain
   */
  virtual double getJacobian(const util::Point3 &p,
                             const std::vector<util::Point3> &nodes,
                             std::vector<std::vector<double>> *J) = 0;

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
