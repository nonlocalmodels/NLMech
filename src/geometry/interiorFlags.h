// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef GEOM_INTERIORFLAGS_H
#define GEOM_INTERIORFLAGS_H

#include "util/point.h" // definition of point
#include <string>
#include <vector>

// forward declaration of interior flags deck
namespace inp {
struct InteriorFlagsDeck;
}

namespace geometry {

/*! @brief An abstraction class to store interior/exterior flags of node
 *
 * This class is used if no-fail region is not active. It simply returns true
 * flag for all the nodes.
 */
class BaseInterior {
public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   * @param bbox Bounding box to determine the location of nodes relative to
   * the boundary
   */
  BaseInterior(inp::InteriorFlagsDeck *deck,
               std::pair<std::vector<double>, std::vector<double>> bbox);

  /*!
   * @brief Returns true for all the nodes
   * @param i Nodal id
   * @param x Nodal coordinate
   * @return True Always returns true as this is default class
   */
  virtual bool getInteriorFlag(const size_t &i, const util::Point3 &x);

protected:
  /*!
   * @brief Interior flags. For given node i the flag is d_intFlags[i%8]. We
   * store flag in 1 bit per node.
   */
  std::vector<uint8_t> d_intFlags;

  /*! @brief Bounding box */
  std::pair<std::vector<double>, std::vector<double>> d_bbox;

  /*! @brief Tolerance to decide if the point is in interior/exterior */
  double d_noFailTol;
};

/*!
 * @brief A class to check if the node is in interior or exterior
 *
 * This class is activated when we do not want to store the flag data for
 * each node, instead we compute the flags as and when needed. This method is
 * only preferred if the memory resource is less or if problem is extremely
 * large that we need to reduce the memory load.
 */
class ComputeInterior : public BaseInterior {
public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   * @param bbox Bounding box to determine the location of nodes relative to
   * the boundary
   */
  explicit ComputeInterior(
      inp::InteriorFlagsDeck *deck,
      std::pair<std::vector<double>, std::vector<double>> bbox);

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
               std::pair<std::vector<double>, std::vector<double>> bbox);

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
 * the material domain or it is near the boundary. This is useful when we
 * implement \a no-fail \a region in \b Peridynamics.
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
  InteriorFlags(inp::InteriorFlagsDeck *deck,
                const std::vector<util::Point3> *nodes,
                const std::pair<std::vector<double>, std::vector<double>>& bbox);

private:
  /*! @brief Class providing interior flags and method */
  BaseInterior *d_interior_p;
};

} // namespace geometry

#endif // GEOM_INTERIORFLAGS_H
