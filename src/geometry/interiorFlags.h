////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef GEOM_INTERIORFLAGS_H
#define GEOM_INTERIORFLAGS_H

#include "util/point.h"           // definition of Point3
#include <string>
#include <vector>

// forward declaration of interior flags deck
namespace inp {
struct InteriorFlagsDeck;
}

namespace geometry {

/*! @brief An abstraction class to store interior/exterior flags of node
 *
 * This is a default class and is used when *no-fail* region is not specified.
 */
class BaseInterior {
public:

  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   * @param bbox Bounding box to determine the location of nodes relative to
   * the boundary
   * @param no_fail_regions No fail regions other than the boundary of domain
   */
  BaseInterior(inp::InteriorFlagsDeck *deck,
               std::pair<std::vector<double>, std::vector<double>> bbox,
               std::vector<std::pair<std::string, std::vector<double>>> no_fail_regions);

  /*!
   * @brief Get interior flag of node
   *
   * Returns true for all the nodes.
   *
   * @param i Nodal id
   * @param x Nodal coordinate
   * @return True Always returns true as this is a default class
   */
  virtual bool getInteriorFlag(const size_t &i, const util::Point3 &x);

protected:
  /*!
   * @brief Interior flags. For given node i the flag is d_intFlags[i%8]. We
   * use 1 bit per node.
   */
  std::vector<uint8_t> d_intFlags;

  /*! @brief Bounding box */
  std::pair<std::vector<double>, std::vector<double>> d_bbox;

  /*! @brief Specify multiple regions in which we set no-fail flag to true */
  std::vector<std::pair<std::string, std::vector<double>>> d_noFailRegions;

  /*! @brief Tolerance to check if the point is in interior/exterior */
  double d_noFailTol;
};

/*!
 * @brief A class to check if the node is in interior or exterior
 *
 * This class checks the relative position of node on the fly instead of
 * dedicating a data for the flags.
 */
class ComputeInterior : public BaseInterior {
public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   * @param bbox Bounding box to determine the location of nodes relative to
   * the boundary
   */
  ComputeInterior(inp::InteriorFlagsDeck *deck,
      std::pair<std::vector<double>, std::vector<double>> bbox,
                  std::vector<std::pair<std::string, std::vector<double>>> no_fail_regions);

  /*!
   * @brief Returns true if node is in interior, false otherwise
   *
   * This function computes and checks the relative position of x with
   * respect to the boundary.
   * @param i Nodal id
   * @param x Nodal coordinate
   * @return True/False If node is in interior or exterior
   */
  bool getInteriorFlag(const size_t &i, const util::Point3 &x) override;
};

/*!
 * @brief A class to check if the node is in interior or exterior
 *
 * This class stores the interior flags of all nodes so that the flags are
 * not computed every time it is needed.
 */
class DataInterior : public BaseInterior {
public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   * @param nodes Pointer to the list of nodal coordinates
   * @param bbox Bounding box to determine the location of nodes relative to
   * the boundary
   */
  DataInterior(inp::InteriorFlagsDeck *deck,
               const std::vector<util::Point3> *nodes,
               std::pair<std::vector<double>, std::vector<double>> bbox,
               std::vector<std::pair<std::string, std::vector<double>>> no_fail_regions);

  /*!
   * @brief Returns true if node is in interior, false otherwise
   *
   * This function simply looks up at the stored bit corresponding to the
   * node and returns true and false depending on whether bit is 0 or 1.
   * @param i Nodal id
   * @param x Nodal coordinate
   * @return True/False If node is in interior or exterior
   */
  bool getInteriorFlag(const size_t &i, const util::Point3 &x) override;
};

/*! @brief A class to store interior/exterior flags of node
 *
 * In this class we store the the flag which indicates if the node is inside
 * the material domain or if it is near the boundary. This required to
 * implement *no-fail region* in Peridynamics.
 */
class InteriorFlags {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   * @param nodes Pointer to the list of nodal coordinates
   * @param bbox Bounding box to determine the location of nodes relative to
   * the boundary
   */
  InteriorFlags(
      inp::InteriorFlagsDeck *deck, const std::vector<util::Point3> *nodes,
      const std::pair<std::vector<double>, std::vector<double>> &bbox);

  /*!
   * @brief Returns true for all the nodes
   * @param i Nodal id
   * @param x Nodal coordinate
   * @return Flag True if node is in interior otherwise false
   */
  bool getInteriorFlag(const size_t &i, const util::Point3 &x);

private:
  /*! @brief Class providing interior flags and method */
  BaseInterior *d_interior_p;
};

} // namespace geometry

#endif // GEOM_INTERIORFLAGS_H
