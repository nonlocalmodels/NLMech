// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef FE_QUADRATURE_H
#define FE_QUADRATURE_H

#include <cstdio>
#include <vector>

#include "../util/point.h" // definition of point

// forward declaration
namespace inp {
struct QuadratureDeck;
}

namespace fe {

/*!
 * @brief A struct to store the quadrature data
 */
struct QuadData {

  /*! @brief Quadrature weight */
  double d_w;

  /*! @brief Quadrature point in 2-d or 3-d */
  util::Point3 d_p;

  /*!
   * @brief Value of shape functions at quad point p.
   *
   * Size will be the number of vertices the element has. E.g. for triangle
   * element shapes will have three entries.
   */
  std::vector<double> d_shapes;

  /*!
   * @brief Value of derivative of shape functions at quad point p.
   *
   * Size will be the number of vertices the element has. E.g. for triangle
   * element shapes will have three entries.
   *
   * x-derivative of ith shape function is d_derShapes[i][0]
   *
   * y-derivative of ith shape function is d_derShapes[i][1]
   *
   */
  std::vector<std::vector<double>> d_derShapes;

  /*!
   * @brief Constructor
   */
  QuadData() : d_w(0.), d_p(util::Point3()){};
};

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
class Quadrature {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck
   * @param element_type Type of element in the mesh
   */
  Quadrature(inp::QuadratureDeck *deck, size_t element_type);

  /*! @brief Turn debug off */
  void turnDebugOff();

  /*! @brief Turn debug on */
  void turnDebugOn();

  /*!
   * @brief Get vector of quadrature data
   *
   * Given vertices of an element, where element type is same as what is
   * specified in the constructor of Quadrature class, it returns the vector
   * of quadrature data.
   *
   * @param nodes Vector of vertices of an element
   * @param order Order of quadrature approximation
   * @return Vector of QuadData
   *
   * @note This depends on the element type that is already set in the
   * constructor of Quadrature class.
   *
   * @note Caller needs to ensure that order does not go higher than 5 as at
   * present only upto fifth order quadrature points are implemented.
   */
  std::vector<fe::QuadData>
  getQuadPoints(const std::vector<util::Point3> &nodes, const int &order);

private:
  /*!
   * @name Methods for triangle element
   */

  /*!
   * @brief Compute the quadrature points for triangle element
   * @param quads Vector quad data to which this function will write to
   * @param order Order of quad point approximation
   */
  void initTri(std::vector<fe::QuadData> &quads, const int &order);

  /*!
   * @brief Returns quadrature data on any given triangle
   * @param nodes Vector of vertices of triangle
   * @param order Order of quadrature approximation
   * @return Vector of quad data
   */
  std::vector<fe::QuadData>
  getQuadPointsTriangle(const std::vector<util::Point3> &nodes,
                        const int &order);

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
  std::vector<double> getTriShapes(const util::Point3 &p);

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
   * @return Vector of derivative of shape functions
   */
  std::vector<std::vector<double>> getTriDerShapes();

  /*!
   * @brief Maps the point on reference triangle to given triangle and
   * returns determinant of Jacobian
   *
   * For linear triangle element, the map from \f$ (\xi, \eta) \f$ to \f$ (x,y)
   * \f$ coordinate is given by
   * \f[ x = N_1(\xi, \eta) x_1 + N_2(\xi, \eta) x_2 + N_3(\xi, \eta) x_3,\f]
   * \f[ y = N_1(\xi, \eta) y_1 + N_2(\xi, \eta) y_2 + N_3(\xi, \eta) y_3,\f]
   * where \f$ N_1, N_2, N_3 \f$ are described in getTriShapes() and \f$ (x_1,
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
  double mapRefTriToTri(util::Point3 &p, const std::vector<double> &shapes,
                        const std::vector<std::vector<double>> &der_shapes,
                        const std::vector<util::Point3> &nodes);

  /** @}*/

  /*!
   * @name Methods for quadrangle element
   */

  /*!
   * @brief Compute the quadrature points for quadrangle element
   * @param quads Vector quad data to which this function will write to
   * @param order Order of quad point approximation
   */
  void initQuadrangle(std::vector<fe::QuadData> &quads, const int &order);

  /*!
   * @brief Returns quadrature data on quadrangle
   * @param nodes Vector of vertices of quadrangle
   * @param order Order of quadrature approximation
   * @return Vector of quad data
   */
  std::vector<fe::QuadData>
  getQuadPointsQuadrangle(const std::vector<util::Point3> &nodes,
                          const int &order);

  /*!
   * @brief Returns the values of shape function at point p for quadrangle
   * element
   *
   * Quadrangle assumed is a reference quadrangle with vertices v1 = (-1,-1),
   * v2 =(1,-1), v3 = (1,1), v4 = (-1,1). The shape function for the reference
   * quadrangle is given by
   * \f[N_1(\xi, \eta) = \frac{(1- \xi)(1 - \eta)}{4}, \, N_2(\xi, \eta) =
   * \frac{(1+ \xi)(1 - \eta)}{4}, \f]
   * \f[
   * N_3(\xi, \eta) = \frac{(1+ \xi)(1 + \eta)}{4}, N_4(\xi, \eta) =
   * \frac{(1- \xi)(1 + \eta)}{4}. \f]
   *
   * @param p Location of point
   * @return Vector of shape functions at point p
   */
  std::vector<double> getQuadShapes(const util::Point3 &p);

  /*!
   * @brief Returns the values of derivative of shape function for
   * quadrangle element at point p
   *
   * For bi-linear quadrangle element, derivative of shape functions are
   * linear and are as follows (for reference quadrangle)
   *
   * \f[\frac{d N_1(\xi, \eta)}{d\xi} = \frac{-(1 - \eta)}{4}, \, \frac{d N_1
   * (\xi, \eta)}{d\eta} = \frac{-(1 - \xi)}{4}, \f]
   * \f[\frac{d N_2(\xi, \eta)}{d\xi} = \frac{(1 - \eta)}{4}, \, \frac{d N_2
   * (\xi, \eta)}{d\eta} = \frac{-(1 + \xi)}{4}, \f]
   * \f[\frac{d N_3(\xi, \eta)}{d\xi} = \frac{(1 + \eta)}{4}, \, \frac{d N_3
   * (\xi, \eta)}{d\eta} = \frac{(1 + \xi)}{4}, \f]
   * \f[\frac{d N_4(\xi, \eta)}{d\xi} = \frac{-(1 + \eta)}{4}, \, \frac{d N_4
   * (\xi, \eta)}{d\eta} = \frac{(1 - \xi)}{4}. \f]
   *
   * @param p Point at which derivatives are to be evaluated
   * @return Vector of derivative of shape functions
   */
  std::vector<std::vector<double>> getQuadDerShapes(const util::Point3 &p);

  /** @}*/

  /*! @brief Order of quadrature point integration approximation */
  int d_quadOrder;

  /*!
   * @brief Order of quadrature point integration approximation for
   * mass matrix
   */
  int d_quadOrderM;

  /*! @brief Number of quadrature points for order d_quadOrder */
  int d_numQuadPts;

  /*! @brief Number of quadrature points for order d_quadOrderM */
  int d_numQuadPtsM;

  /*! @brief Type of element */
  size_t d_eType;

  /*! @brief Quadrature data collection */
  static std::vector<fe::QuadData> d_quads;

  /*! @brief Quadrature data collection for mass matrix */
  static std::vector<fe::QuadData> d_quadsM;

  /*! @brief Debug flag (0-off, 1-on) */
  static int d_debug;
};

} // namespace fe

#endif // FE_QUADRATURE_H
