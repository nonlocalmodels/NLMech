////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef INP_MESHDECK_H
#define INP_MESHDECK_H

#include <string>
#include "util/utilIO.h"

namespace inp {

/**
 * \ingroup Input
 */
/**@{*/

/*! @brief Structure to read and store mesh related input data */
struct MeshDeck {

  /*! @brief Dimension */
  size_t d_dim;

  /*!
   * @brief Tag for spatial discretization
   *
   * List of allowed values are:
   * - \a finite_difference
   * - \a weak_finite_element
   * - \a nodal_finite_element
   * - \a truss_finite_element
   */
  std::string d_spatialDiscretization;

  /*! @brief Filename to read mesh data */
  std::string d_filename;

  /*! @brief Flag which indicates if mesh size is to be computed */
  bool d_computeMeshSize;

  /*! @brief Mesh size */
  double d_h;

  /*!
   * @brief Centroid based discretization
   */
  bool d_isCentroidBasedDiscretization;

  /*!
   * @brief Specify if we keep the element-node connectivity data
   */
  bool d_keepElementConn;

  /*!
   * @brief Specify if coupling data from PUM is loaded
   */
  bool d_loadPUMData;

  /*!
   * @brief Constructor
   */
  MeshDeck()
      : d_dim(0), d_computeMeshSize(false), d_h(0.),
        d_isCentroidBasedDiscretization(false), d_keepElementConn(false){};

  /*!
   * @brief Returns the string containing information about the instance of
   * the object
   *
   * @param nt Number of tabs to append before each line of string
   * @param lvl Level of information sought (higher level means more
   * information)
   * @return string String containing information about this object
   * */
  std::string printStr(int nt = 0, int lvl = 0) const {
    auto tabS = util::io::getTabS(nt);
    std::ostringstream oss;
    oss << tabS << "------- MeshDeck --------" << std::endl << std::endl;
    oss << tabS << "Dimension = " << d_dim << std::endl;
    oss << tabS << "Spatial discretization type = " << d_spatialDiscretization << std::endl;
    oss << tabS << "Mesh filename = " << d_filename << std::endl;
    oss << tabS << "Compute mesh size = " << d_computeMeshSize << std::endl;
    oss << tabS << "Mesh size = " << d_h << std::endl;
    oss << tabS << "Is this centroid-based particle mesh = " << d_isCentroidBasedDiscretization <<
        std::endl;
    oss << tabS << std::endl;

    return oss.str();
  };

  /*!
   * @brief Prints the information about the instance of the object
   *
   * @param nt Number of tabs to append before each line of string
   * @param lvl Level of information sought (higher level means more
   * information)
   * */
  void print(int nt = 0, int lvl = 0) const { std::cout << printStr(nt, lvl); };
};

/** @}*/

} // namespace inp

#endif // INP_MESHDECK_H
