// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef FE_QUADDATA_H
#define FE_QUADDATA_H

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

} // namespace fe

#endif // FE_QUADDATA_H
