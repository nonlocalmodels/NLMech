// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <string>
#include "../util/point.h"         // definition of struct Point3
#include "../util/matrixBlaze.h"   // definition of SymMatrixFij

// forward declaration of geometry deck
namespace inp {
struct GeometryDeck;
}

//!Collection of methods and database purely related to geometry
namespace geometry {

/*! @brief Methods and database associated to the mesh */
class Geometry {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   */
  explicit Geometry(inp::GeometryDeck *deck);

  /**
   * \defgroup Accessor methods
   */
  /**@{*/

  /*!
   * @brief Return the dimension of the domain
   * @return N dimension
   */
  size_t getDimension();

  /*!
   * @brief Return the number of nodes
   * @return N number of nodes
   */
  size_t getNumNodes();

  /*!
   * @brief Return the number of dofs
   * @return N number of dofs
   */
  size_t getNumDofs();

  /*!
   * @brief Return coordinates of node i
   * @return Coordinates
   */
  util::Point3 getNode(size_t i);

  /*!
   * @brief Return nodes data
   * @return nodes vector of nodal coordinates
   */
  std::vector<util::Point3> getNodes();

  /*!
   * @brief Return the pointer to nodes data
   * @return Pointer to nodes data
   */
  const std::vector<util::Point3>* getNodesP();

  /** @}*/

private:

  /**
   * \defgroup Utility methods
   */
  /**@{*/

  /*! @brief Reads mesh data from file */
  void readFile(std::string filename);

  /** @}*/

  /**
   * \defgroup Mesh related data
   */
  /**@{*/

  /*! @brief Number of nodes */
  size_t d_numNodes;

  /*! @brief Number of elements */
  size_t d_numElems;

  /*! @brief Element type. We follow VTK convention to identify the elements:
   * Line element = 3,
   * Triangle element = 5,
   * Pixel element = 8,
   * Quadrilateral element = 9,
   * Tetrahedral element = 10
   */
  size_t d_eType;

  /*! @brief Number of vertex per element.
   * Line element d_eNumVertex = 2,
   * Triangle element d_eNumVertex = 3,
   * Quadrilateral element d_eNumVertex = 4,
   * Tetrahedral element d_eNumVertex = 4
   */
  size_t d_eNumVertex;

  /*! @brief Vector of initial (reference) coordinates of nodes */
  std::vector<util::Point3> d_nodes;

  /*! @brief Element-node connectivity data. I.e. list of ids of nodes which
   * form the element.
   * First d_eNumVertex data are first element, and next d_eNumVertex are for
   * second element, and so on
   */
  std::vector<size_t> d_enc;

  /*! @brief Node-element connectivity data. I.e. list of ids of elements
   * which have particular node as its vertex
   */
  std::vector<std::vector<size_t>> d_nec;

  /*! @brief Vector of fixity mask of each node.
   * First bit represents x-dof, second represents y-dof, and third
   * represents z-dof
   */
  std::vector<char> d_fix;

  /*! @brief Vector of volume of each node. For uniform square mesh, the
   * volume is simply h^2 (h^3 in 3-d), however, for general mesh, the volume
   * is computed using linear interpolation
   */
  std::vector<double> d_vol;

  /** @}*/

  /**
   * \defgroup Data which do not belong to mesh directly
   */
  /**@{*/

  /*! @brief Geometry deck */
  inp::GeometryDeck *d_geometryDeck_p;

  /*! @brief Number of dofs = Dimension times number of nodes */
  size_t d_numDofs;

  /*! @brief Inverse of diagonal mass matrix stored in a vector */
  std::vector<double> d_invMDiag;

  /*! @brief Inverse of exact mass matrix stored in a Blaze matrix data type.
   * We store data in float to reduce memory load. However if results appear
   * to be not so accurate, this should be changed to double
   */
  util::SymMatrixFij d_invM;

  /*! @brief Map from global reduced id to default global id. Each free dof
   * has associated global id, which we refer to as "global reduced id",
   * and d_gMap provides a map from global reduced id to default global id.
   * Use: Only used when the discretization is "weak_finite_element". Used in
   * assembling of mass matrix and also in time integration
   */
  std::vector<size_t> d_gMap;

  /*! @brief Map from Default global id to global reduced global id.
   * This is inverse of d_gMap
   */
  std::vector<int> d_gInvMap;

  /*! @brief Bounding box */
  std::pair<std::vector<double>, std::vector<double>> d_bbox;

  /** @}*/

};

} // namespace geometry

#endif // GEOMETRY_H
