////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef RW_VTKREADER_H
#define RW_VTKREADER_H

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>

#include "util/point.h"           // definition of Point3

namespace rw {

namespace reader {

/*!
 * @brief A class to read vtk (vtu) mesh files
 *
 * @note Depends on VTK library.
 */
class VtkReader {

public:
  /*!
   * @brief Constructor
   * @param filename name of mesh file
   */
  explicit VtkReader(const std::string &filename);

  /*!
   * @brief Reads mesh data into node file and element file
   * @param dim Dimension
   * @param nodes vector of nodes data
   * @param element_type type of element
   * @param num_elem number of elements
   * @param enc vector holding element-node connectivity
   * @param nec vector holding node-element connectivity
   * @param volumes vector holding volume of the nodes
   * @param is_fd flag indicating if this mesh is for finite_difference
   * simulation
   */
  void readMesh(size_t dim, std::vector<util::Point3> *nodes,
                size_t &element_type, size_t &num_elem,
                std::vector<size_t> *enc, std::vector<std::vector<size_t>> *nec,
                std::vector<double> *volumes, bool is_fd = false);

  /*!
   * @brief Reads nodal position
   *
   * @param nodes vector of nodes data
   */
  void readNodes(std::vector<util::Point3> *nodes);

  /*!
   * @brief reads point data from .vtu file
   * @param name Name of data
   * @param data Pointer to the vector of data
   * @return status True if data is found otherwise false
   */
  bool readPointData(std::string name, std::vector<util::Point3> *data);

  /*! @brief Close the file */
  void close();

private:
  /**
   * @name Internal data
   */
  /**@{*/

  /*! @brief Counter */
  static size_t d_count;

  /*! @brief XML unstructured grid writer */
  vtkSmartPointer<vtkXMLUnstructuredGridReader> d_reader_p;

  /*! @brief Unstructured grid */
  vtkSmartPointer<vtkUnstructuredGrid> d_grid_p;

  /** @}*/
};

} // namespace reader

} // namespace rw

#endif // RW_VTKREADER_H
