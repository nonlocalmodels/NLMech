////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef RW_VTKREADER_H
#define RW_VTKREADER_H

#include <util/matrixBlaze.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>

#include "util/matrix.h" // definition of matrices
#include "util/matrixBlaze.h"
#include "util/point.h" // definition of Point3

namespace rw {

namespace reader {

/*!
 * @brief A class to read VTK (.vtu) mesh files
 *
 * @note Depends on VTK library.
 */
class VtkReader {

public:
  /*!
   * @brief Constructor
   * @param filename Name of mesh file
   *
   * @note filename is expected to have .vtu extension
   */
  explicit VtkReader(const std::string &filename);

  /*!
   * @brief Checks if file has needed data
   * @param data_tag Tag name of data
   * @return True If file has data otherwise false
   */
  bool vtuHasPointData(const std::string &data_tag);

  /*!
   * @brief Checks if file has needed data
   * @param data_tag Tag name of data
   * @return True If file has data otherwise false
   */
  bool vtuHasCellData(const std::string &data_tag);

  /*!
   * @brief Reads all point data tags
   * @return List List of point data tags
   */
  std::vector<std::string> readVtuFilePointTags();

  /*!
   * @brief Reads all cell data tags
   * @return List List of point data tags
   */
  std::vector<std::string> readVtuFileCellTags();

  /*!
   * @brief Reads mesh data into node file and element file
   * @param dim Dimension
   * @param nodes Vector of nodes data
   * @param element_type Type of element
   * @param num_elem Number of elements
   * @param enc Vector holding element-node connectivity
   * @param nec Vector holding node-element connectivity
   * @param volumes Vector holding volume of the nodes
   * @param is_fd Flag indicating if this mesh is for finite_difference
   * simulation
   */
  void readMesh(size_t dim, std::vector<util::Point3> *nodes,
                size_t &element_type, size_t &num_elem,
                std::vector<size_t> *enc, std::vector<std::vector<size_t>> *nec,
                std::vector<double> *volumes, bool is_fd = false);

  /*!
   * @brief Reads nodal position
   *
   * @param nodes Vector of nodal coordinates
   */
  void readNodes(std::vector<util::Point3> *nodes);

  /*!
   * @brief Reads cell data, i.e. element-node connectivity and node-element
   * connectivity
   * @param dim Dimension
   * @param element_type Type of element
   * @param num_elem Number of elements
   * @param enc Element-node connectivity
   * @param nec Node-element connectivity
   */
  void readCells(size_t dim, size_t &element_type, size_t &num_elem,
                 std::vector<size_t> *enc,
                 std::vector<std::vector<size_t>> *nec);

  /*!
   * @brief reads point data from .vtu file
   * @param name Name of data
   * @param data Pointer to the vector of data
   * @return status True if data is found otherwise false
   */
  bool readPointData(const std::string &name, std::vector<uint8_t> *data);

  /*!
   * @brief reads point data from .vtu file
   * @param name Name of data
   * @param data Pointer to the vector of data
   * @return status True if data is found otherwise false
   */
  bool readPointData(const std::string &name, std::vector<size_t> *data);

  /*!
   * @brief reads point data from .vtu file
   * @param name Name of data
   * @param data Pointer to the vector of data
   * @return status True if data is found otherwise false
   */
  bool readPointData(const std::string &name, std::vector<int> *data);

  /*!
   * @brief reads point data from .vtu file
   * @param name Name of data
   * @param data Pointer to the vector of data
   * @return status True if data is found otherwise false
   */
  bool readPointData(const std::string &name, std::vector<float> *data);

  /*!
   * @brief reads point data from .vtu file
   * @param name Name of data
   * @param data Pointer to the vector of data
   * @return status True if data is found otherwise false
   */
  bool readPointData(const std::string &name, std::vector<double> *data);

  /*!
   * @brief reads point data from .vtu file
   * @param name Name of data
   * @param data Pointer to the vector of data
   * @return status True if data is found otherwise false
   */
  bool readPointData(const std::string &name, std::vector<util::Point3> *data);

  /*!
   * @brief reads point data from .vtu file
   * @param name Name of data
   * @param data Pointer to the vector of data
   * @return status True if data is found otherwise false
   */
  bool readPointData(const std::string &name,
                     std::vector<util::SymMatrix3> *data);

  /*!
   * @brief reads point data from .vtu file
   * @param name Name of data
   * @param data Pointer to the vector of data
   * @return status True if data is found otherwise false
   */
  bool readPointData(const std::string &name,
                     std::vector<util::Matrix33> *data);

  /*!
   * @brief reads cell data from .vtu file
   * @param name Name of data
   * @param data Pointer to the vector of data
   * @return status True if data is found otherwise false
   */
  bool readCellData(const std::string &name, std::vector<float> *data);

  /*!
   * @brief reads cell data from .vtu file
   * @param name Name of data
   * @param data Pointer to the vector of data
   * @return status True if data is found otherwise false
   */
  bool readCellData(const std::string &name, std::vector<double> *data);

  /*!
   * @brief reads cell data from .vtu file
   * @param name Name of data
   * @param data Pointer to the vector of data
   * @return status True if data is found otherwise false
   */
  bool readCellData(const std::string &name, std::vector<util::Point3> *data);

  /*!
   * @brief reads cell data from .vtu file
   * @param name Name of data
   * @param data Pointer to the vector of data
   * @return status True if data is found otherwise false
   */
  bool readCellData(const std::string &name,
                    std::vector<util::SymMatrix3> *data);

  /*!
   * @brief reads cell data from .vtu file
   * @param name Name of data
   * @param data Pointer to the vector of data
   * @return status True if data is found otherwise false
   */
  bool readCellData(const std::string &name, std::vector<util::Matrix33> *data);

  /*! @brief Close the file */
  void close();

private:
  /*! @brief Counter */
  static size_t d_count;

  /*! @brief XML unstructured grid writer */
  vtkSmartPointer<vtkXMLUnstructuredGridReader> d_reader_p;

  /*! @brief Unstructured grid */
  vtkSmartPointer<vtkUnstructuredGrid> d_grid_p;
};

} // namespace reader

} // namespace rw

#endif // RW_VTKREADER_H
