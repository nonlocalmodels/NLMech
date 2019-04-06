// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef FE_MESH_H
#define FE_MESH_H

#include <hpx/config.hpp>
#include <algorithm>

#include "util/point.h" // definition of struct Point3
#include <string>
#include <vector>

// forward declaration of geometry deck
namespace inp {
struct MeshDeck;
}

/*!
 * @brief Collection of methods and data related to mesh and fem
 *
 * This namespace provides methods and data members specific to mesh and
 * finite element.
 *
 * @sa Mesh, MassMatrix, Quadrature
 */
namespace fe {

/*! @brief A class for mesh data
 *
 * In this class the mesh data such as nodes, element-node connectivity,
 * node-element connectivity are stored. The class also stores fixity mask of
 * nodes which indicate if x-, y-, or z-dof of the node is fixed or free.
 *
 * We currently only support mesh with only one type of elements, i.e. mesh
 * can not have mix of two types of elements. For example, we can not have
 * mesh with triangle and quadrangle elements together.
 *
 * This class is used in both finite difference implementation and finite
 * element implementation. For finite difference, we only require nodal
 * volume. If the mesh file contains nodal volume, we skip reading
 * element-node and node-element connectivity, however if mesh file does not
 * have nodal volume data, we read connectivity data and compute nodal volume
 * and then delete the connectivity data.
 */
class Mesh {

public:
  /*!
   * @brief Constructor
   *
   * The constructor initializes the data using input deck, performs checks
   * on input data, and reads mesh file and populates the mesh related data.
   * The mesh file of .csv, .vtu and .msh are supported.
   *
   * @param deck Input deck which contains user-specified information
   *
   * @sa readFile()
   */
  explicit Mesh(inp::MeshDeck *deck);

  /**
   * @name Accessor methods
   */
  /**@{*/

  /*!
   * @brief Get the dimension of the domain
   * @return N dimension
   */
  size_t getDimension();

  /*!
   * @brief Get the number of nodes
   * @return N number of nodes
   */
  size_t getNumNodes();

  /*!
   * @brief Get the number of dofs
   * @return N number of dofs
   */
  size_t getNumDofs();

  /*!
   * @brief Get the type of element in mesh
   * @return N element type (using VTK convention)
   */
  size_t getElementType();

  /*!
   * @brief Get coordinates of node i
   * @param i Id of the node
   * @return Coordinates
   */
  util::Point3 getNode(size_t i);

  /*!
   * @brief Get nodes data
   * @return nodes vector of nodal coordinates
   */
  std::vector<util::Point3> getNodes();

  /*!
   * @brief Get the pointer to nodes data
   * @return Pointer to nodes data
   */
  const std::vector<util::Point3> *getNodesP();

  /*!
   * @brief Get the connectivity of element
   *
   * Since we store connectivity in single vector, we use d_eNumVertex to get
   * the connectivity of element. Given element e, the connectivity of e
   * begins from location \f[ e*d\_eNumVertex + 0 \f] upto \f[e*d\_eNumVertex
   * + d\_eNumVertex - 1\f]
   *
   * So connectivity of e is
   *
   * d_enc[e*d_eNumVertex+0], d_enc[e*d_eNumVertex+1], ...,
   * d_end[e*d_eNumVertex+d_eNumVertex-1]
   *
   * @param i Id of an element
   * @return Vec vector of nodal ids
   */
  std::vector<size_t> getElementConnectivity(size_t i);

  /*!
   * @brief Get the vertices of element
   *
   * @param i Id of an element
   * @return Vec vector of vertices
   */
  std::vector<util::Point3> getElementConnectivityNodes(size_t i);

  /*!
   * @brief Get the bounding box of the mesh
   * @return Box Bounding box
   */
  std::pair<std::vector<double>, std::vector<double>> getBoundingBox();

  /** @}*/

private:
  /**
   * @name Utility methods
   */
  /**@{*/

  /*!
   * @brief Reads mesh data from file and populates other data
   *
   * This function calls reader methods in namespace rw to read the mesh file
   * . For finite difference implementation, .csv mesh file with only nodal
   * coordinates and nodal volumes data. However, for finite element
   * implementation, we require either .vtu or .msh file with connectivity
   * data.
   *
   * @sa rw::reader::readCsvFile(), rw::reader::readVtuFile(),
   * rw::reader::readMshFile()
   *
   * @param filename Name of the mesh file
   * */
  void createData(std::string filename);

  /*!
   * @brief Compute the nodal volume
   *
   * This method uses finite element mesh and computes the nodal volume.
   * Formula for volume of a node \f$ i\f$ is given by
   * \f[ V_i = \sum_{e \in N_i} \int_{T_e} N_i(x) dx, \f]
   * where \f$N_i\f$ is the list of elements which have node \f$ i\f$ as its
   * vertex, \f$ T_e\f$ is the element domain, \f$ N_i\f$ is the shape
   * function of the node \f$ i\f$.
   */
  void computeVol();

  /** @}*/

  /**
   * @name Mesh data
   */
  /**@{*/

  /*! @brief Number of nodes */
  size_t d_numNodes;

  /*! @brief Number of elements */
  size_t d_numElems;

  /*! @brief Element type
   *
   * We follow VTK convention to identify the elements:
   * - Line element = 3,
   * - Triangle element = 5,
   * - Pixel element = 8,
   * - Quadrilateral element = 9,
   * - Tetrahedral element = 10
   */
  size_t d_eType;

  /*! @brief Number of vertex per element
   *
   * This information is useful in getting the connectivity for given element
   * . We assume that mesh has only type of elements and based on that
   * assumption we store the element-node connectivity in a single vector.
   *
   * - Line element: 2,
   * - Triangle element: 3,
   * - Quadrilateral element: 4,
   * - Tetrahedral element: 4
   */
  size_t d_eNumVertex;

  /*! @brief Vector of initial (reference) coordinates of nodes
   *
   * We use struct Point3 which consists of three double data and comes with
   * length() and dot() function. It also provides operators.
   *
   * @sa util::Point3
   */
  std::vector<util::Point3> d_nodes;

  /*! @brief Element-node connectivity data
   *
   * First d_eNumVertex data gives the connectivity of first element, and next
   * d_eNumVertex data gives the connectivity of second element and so on and
   * so fourth.
   */
  std::vector<size_t> d_enc;

  /*! @brief Node-element connectivity data */
  std::vector<std::vector<size_t>> d_nec;

  /*! @brief Vector of fixity mask of each node
   *
   * First bit represents x-dof, second represents y-dof, and third
   * represents z-dof. To check if x-dof of \f$ i^{th} \f$ node is fixed, we
   * check d_fix[i] & FIX_X_MASK and if it is true then x-dof is fixed.
   * Similarly we check for y-dof and z-dof using FIX_Y_MASK and FIX_Z_MASK.
   *
   * We store data in uint8_t type which is 1 byte. Although we only need 3 bits.
   */
  std::vector<uint8_t> d_fix;

  /*! @brief Vector of volume of each node
   *
   * For uniform square mesh, the volume is simply \f$ h^2 \f$ in 2-d and \f$
   * h^3\f$ in 3-d, where \f$ h\f$ is the mesh size. For general mesh, the
   * volume is computed using the finite element mesh.
   */
  std::vector<double> d_vol;

  /** @}*/

  //  /*! @brief Mesh deck */
  //  inp::MeshDeck *d_meshDeck_p;

  /*! @brief Dimension */
  size_t d_dim;

  /*! @brief Tag for spatial discretization
   *
   * List of allowed values are:
   * - finite_difference
   * - weak_finite_element
   * - nodal_finite_element
   * - truss_finite_element
   */
  std::string d_spatialDiscretization;

  /*! @brief Filename to read mesh data */
  std::string d_filename;

  /*! @brief Number of dofs = Dimension times number of nodes */
  size_t d_numDofs;

  /*! @brief Map from global reduced id to default global id
   *
   * Each free dof has associated global id, which we refer to as "global
   * reduced id", and d_gMap provides a map from global reduced id to
   * default global id.
   *
   * @note Needed only when the discretization is "weak_finite_element" for
   * assembly of the mass matrix.
   */
  std::vector<size_t> d_gMap;

  /*! @brief Map from Default global id to global reduced global id.
   *
   * This is inverse of d_gMap
   */
  std::vector<int> d_gInvMap;

  /*! @brief Bounding box */
  std::pair<std::vector<double>, std::vector<double>> d_bbox;
};

} // namespace fe

#endif // FE_MESH_H
