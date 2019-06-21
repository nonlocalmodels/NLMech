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
 * They are specific to 2-point line element.
 *
 * The reference line element \f$T^0 \f$ is made of two vertex at points \f$
 * v^1 = -1, v^2 = 1 \f$.
 *
 * 1. The shape functions at point \f$ \xi \in T^0 \f$ are
 * \f[N^0_1(\xi) = \frac{1 - \xi}{2}, \quad N^0_2(\xi) = \frac{1 + \xi}{2}. \f]
 *
 * 2. Derivative of shape functions are constant and are as follows
 * \f[\frac{d N^0_1(\xi)}{d\xi} = \frac{-1}{2}, \, \frac{d N^0_2(\xi)}{d\xi}
 * = \frac{1}{2}. \f]
 *
 * 3. Map \f$ \xi \in T^0 \to x \in T \f$ is given by
 * \f[ x(\xi) = \sum_{i=1}^2 N^0_i(\xi) v^i_x, \f]
 * where \f$ v^1, v^2\f$ are vertices of element \f$ T \f$. For 1-d
 * points, we simply have \f$ v^i_x = v^i \f$ as \f$ v^i \f$ is a vector with
 * only 1 element.
 *
 * 4. Jacobian of the map \f$ \xi \in T^0 \to x \in T \f$ is given by
 * \f[ J = \frac{dx}{d\xi}. \f]
 * Since this 1-d case, Jacobian and its determinant are same. For line element
 * Jacobian is and is given by
 * \f[ J = \frac{dx}{d\xi} = \frac{v^2_x - v^1_x}{2} = \frac{length(T)
 * }{length(T^0)}. \f]
 *
 * 5. Inverse map from \f$ x \in T \f$ to \f$ \xi \in T^0\f$ is given by
 * \f[ \xi(x) = \frac{2}{v^2_x - v^1_x} (x - \frac{v^2_x + v^1_x}{2}) =
 * \frac{1}{J}(x - \frac{v^2_x + v^1_x}{2}). \f]
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
   * @brief Returns the length of element
   *
   * If line \f$ T \f$ is given by points \f$ v^1, v^2\f$ then the length is
   * simply
   * \f[ length(T) = v^2_x - v^1_x. \f]
   *
   * @param nodes Vertices of element
   * @return Vector of shape functions at point p
   */
  double elemSize(const std::vector<util::Point3> &nodes) override;

  /*!
   * @brief Returns the values of shape function at point p
   *
   * Line \f$ T \f$ is given by points \f$ v^1, v^2\f$. We first map
   * the point p in \f$ T \f$ to reference line \f$ T^0 \f$ using
   * LineElem::mapPointToRefElem and then compute shape functions at the
   * mapped point using LineElem::getShapes(const util::Point3 &).
   *
   * @param p Location of point
   * @param nodes Vertices of element
   * @return Vector of shape functions at point p
   */
  std::vector<double>
  getShapes(const util::Point3 &p,
            const std::vector<util::Point3> &nodes) override;

  /*!
   * @brief Returns the values of derivative of shape function at point p
   *
   * Let \f$ x\f$ is the point on line \f$ T \f$ and
   * let \f$ \xi \f$ is the point on reference line \f$ T^0 \f$.
   * We let shape function on \f$ T\f$ denote as \f$ N_1, N_2 \f$ and
   * shape functions on \f$ T^0 \f$ denote as \f$ N^0_1, N^0_2 \f$.
   *
   * We are interested in \f$ \frac{\partial N_i(x_p)}{\partial x}\f$. By
   * using the map \f$ \xi \to x \f$ (see LineElem) we have
   * \f[ N^0_i(\xi) = N_i(x(\xi)) \f]
   * and therefore we can write
   * \f[ \frac{\partial N^0_i(\xi)}{\partial \xi} = \frac{\partial
   * N_i}{\partial x} \frac{\partial x}{\partial \xi}. \f]
   * Since \f$ \frac{\partial x}{\partial \xi} = J \f$ is Jacobian of map, we
   * can compute it easily given \f$ v^1, v^2 \f$. We then have
   * \f[ \frac{\partial N_i(\xi)}{\partial x} = \frac{1}{J} \frac{\partial
   * N^0_i(\xi)}{\partial \xi}. \f]
   *
   * @param p Location of point
   * @param nodes Vertices of element
   * @return Vector of derivative of shape functions
   */
  std::vector<std::vector<double>>
  getDerShapes(const util::Point3 &p,
               const std::vector<util::Point3> &nodes) override;

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
   * Let \f$ T \f$ is the given line with vertices \f$ v^1, v^2 \f$ and let
   * \f$ T^0 \f$ is the reference line.
   *
   * 1. To compute quadrature point, we first compute the quadrature points
   * on reference line \f$ T^0 \f$, and then map the point on reference
   * line to the current line \f$ T \f$. The map from \f$ \xi\in T^0 \f$ to
   * \f$ x \in T \f$ is given by
   * \f[ x = \sum_{i=1}^2 N^0_i(\xi) v^i_x,\f]
   * where \f$ N^0_1, N^0_2\f$ are shape function associated to the
   * reference line \f$ T^0 \f$ and described in LineElem.
   *
   * 2. To compute the quadrature weight, we compute the quadrature weight
   * associated to the quadrature point in reference line \f$ T^0 \f$.
   * Suppose \f$ w^0_q \f$ is the quadrature weight associated to quadrature
   * point \f$ \xi_q \in T^0 \f$, then the quadrature point \f$ w_q \f$
   * associated to the mapped point \f$ x(\xi_q) \in T \f$ is given by
   * \f[ w_q = w^0_q * J \f]
   * where \f$ J \f$ is the Jacobian of map from \f$ \xi \in T^0 \f$ to \f$ x
   * \in T \f$.
   *
   * 3. We compute shape functions \f$ N_1, N_2\f$ associated to \f$ T \f$ at
   * quadrature point \f$ x(\xi_q) \f$ using formula
   * \f[ N_i(x(\xi_q)) = N^0_i(\xi_q). \f]
   *
   * 5. To compute derivative of shape functions \f$ \frac{\partial
   * N_1}{\partial x}, \frac{\partial N_2}{\partial x}\f$ associated to \f$ T
   * \f$, we use the relation between derivatives of shape function in \f$ T
   * \f$ and \f$ T^0 \f$ described in LineElem::getDerShapes.
   *
   * @param nodes Vector of vertices of an element
   * @return Vector of QuadData
   */
  std::vector<fe::QuadData>
  getQuadDatas(const std::vector<util::Point3> &nodes) override;

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
   * This function is lite version of TriElem::getQuadDatas.
   *
   * @param nodes Vector of vertices of an element
   * @return Vector of QuadData
   */
  std::vector<fe::QuadData>
  getQuadPoints(const std::vector<util::Point3> &nodes) override;

private:
  /*!
   * @brief Returns the values of shape function at point p on reference element
   *
   * @param p Location of point
   * @return Vector of shape functions at point p
   */
  std::vector<double> getShapes(const util::Point3 &p) override;

  /*!
   * @brief Returns the values of derivative of shape function at point p on
   * reference element
   *
   * @param p Location of point
   * @return Vector of derivative of shape functions
   */
  std::vector<std::vector<double>> getDerShapes(const util::Point3 &p) override;

  /*!
   * @brief Maps point p in a given element to the reference element
   *
   * Let \f$ v^1, v^2\f$ are vertices of element \f$ T\f$ and let
   * \f$T^0 \f$ is the reference element. Map \f$ x \in T \f$ to \f$
   * \xi \in T^0 \f$ is given by
   * \f[ \xi(x) = \frac{2}{v^2_x - v^1_x} (x - \frac{v^2_x + v^1_x}{2}) =
   * \frac{1}{J}(x - \frac{v^2_x + v^1_x}{2}). \f]
   *
   * If mapped point \f$ \xi\f$ does not satisfy condition
   * - \f[ -1 \leq \xi \leq 1 \f]
   * then the point \f$ \xi \f$ does not belong to reference line \f$ T^0\f$
   * or equivalently point \f$x  \f$ does not belong to line \f$ T \f$
   * and the method issues error. Otherwise the method returns point \f$ \xi\f$.
   *
   * @param p Location of point
   * @param nodes Vertices of element
   * @return Vector of shape functions at point p
   */
  util::Point3
  mapPointToRefElem(const util::Point3 &p,
                    const std::vector<util::Point3> &nodes) override;

  /*!
   * @brief Computes Jacobian of map from reference element \f$ T^0 \f$ to
   * given element \f$ T \f$
   *
   * @param p Location of point in reference element
   * @param nodes Vertices of element
   * @param J Matrix to store the Jacobian (if not nullptr)
   * @return det(J) Determinant of the Jacobain (same as Jacobian in 1-d)
   */
  double getJacobian(const util::Point3 &p,
                     const std::vector<util::Point3> &nodes,
                     std::vector<std::vector<double>> *J) override;

  /*!
   * @brief Compute the quadrature points for triangle element
   */
  void init() override;

};

} // namespace fe

#endif // FE_LINEELEM_H
