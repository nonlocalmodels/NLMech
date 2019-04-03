// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef RW_VTKREADER_H
#define RW_VTKREADER_H

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>

#include "../util/point.h"

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
