// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef RW_WRITER_H
#define RW_WRITER_H

#include "util/point.h"           // definition of Point3
#include <string>
#include <vector>

// forward declaration
namespace rw {
namespace writer {
class VtkWriter;
}
} // namespace rw

namespace rw {

/*!
 * @brief Collection of methods and database related to writing
 *
 * This namespace provides methods and data members specific to wiriting of
 * the mesh data and simulation data. Currently, .vtu is supported.
 */
namespace writer {

/*!
 * @brief A interface class writing data using vtk writer
 *
 * This interface separates the caller from vtk library.
 */
class VtkWriterInterface {

public:
  /*!
   * @brief Constructor
   *
   * Creates and opens .vtu file of name given by filename. The file remains
   * open till the close() function is invoked.
   *
   * @param filename Name of file which will be created
   */
  explicit VtkWriterInterface(const std::string &filename);

  /*! @brief Destructor */
  ~VtkWriterInterface();

  /**
   * @name Mesh data
   */
  /**@{*/

  /*!
   * @brief Writes the nodes to the file
   * @param nodes Current positions of the nodes
   */
  void appendNodes(const std::vector<util::Point3> *nodes);

  /*!
   * @brief Writes the nodes to the file
   * @param nodes Reference positions of the nodes
   * @param u Nodal displacements
   */
  void appendNodes(const std::vector<util::Point3> *nodes,
                   const std::vector<util::Point3> *u);

  /*!
   * @brief Writes the mesh data to file
   *
   * @param nodes Vector of nodal coordinates
   * @param element_type Type of element
   * @param en_con Vector of element-node connectivity
   * @param u Vector of nodal displacement
   */
  void appendMesh(const std::vector<util::Point3> *nodes,
                  const size_t &element_type,
                  const std::vector<size_t> *en_con,
                  const std::vector<util::Point3> *u);

  /** @}*/

  /**
   * @name Point data
   */
  /**@{*/

  /*!
   * @brief Writes the scalar point data to the file
   * @param name Name of the data
   * @param data The vector containing the data
   */
  void appendPointData(const std::string &name,
                       const std::vector<uint8_t> *data);

  /*!
   * @brief Writes the scalar point data to the file
   * @param name Name of the data
   * @param data The vector containing the data
   */
  void appendPointData(const std::string &name,
                       const std::vector<size_t> *data);

  /*!
   * @brief Writes the scalar point data to the file
   * @param name Name of the data
   * @param data The vector containing the data
   */
  void appendPointData(const std::string &name, const std::vector<int> *data);

  /*!
   * @brief Writes the scalar point data to the file
   * @param name Name of the data
   * @param data The vector containing the data
   */
  void appendPointData(const std::string &name, const std::vector<float> *data);

  /*!
   * @brief Writes the scalar point data to the file
   * @param name Name of the data
   * @param data The vector containing the data
   */
  void appendPointData(const std::string &name,
                       const std::vector<double> *data);

  /*!
   * @brief Writes the vector point data to the file
   * @param name Name of the data
   * @param data The vector containing the data
   */
  void appendPointData(const std::string &name,
                       const std::vector<util::Point3> *data);

  /** @}*/

  /**
   * @name Field data
   */
  /**@{*/

  /*!
   * @brief Writes the scalar field data to the file
   * @param name Name of the data
   * @param data Value
   */
  void appendFieldData(const std::string &name, const double &data);

  /*!
   * @brief Writes the scalar field data to the file
   * @param name Name of the data
   * @param data Value
   */
  void appendFieldData(const std::string &name, const float &data);

  /*!
   * @brief Writes the time step to the file
   * @param timestep Current time step of the simulation
   */
  void addTimeStep(const double &timestep);

  /** @}*/

  /*!
   * @brief Closes the file and store it to the hard disk
   */
  void close();

private:
  /*! @brief Pointer to the vtk writer class */
  rw::writer::VtkWriter *d_vtkWriter_p;
};

} // namespace writer

} // namespace rw

#endif // RW_WRITER_H
