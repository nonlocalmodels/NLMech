////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

//Todo: make member functions private after refactoring

#ifndef RW_WRITER_H
#define RW_WRITER_H

#include "util/point.h"           // definition of Point3
#include "util/matrix.h"          // definition of SymMatrix3
#include <string>
#include <vector>

namespace inp {
class Input;
struct ModelDeck;
struct OutputDeck;
} // namespace inp

namespace geometry {
class Neighbor;
} // namespace inp


// forward declaration
namespace data {
class DataManager;
}

// forward declaration
namespace rw {
namespace writer {
class VtkWriter;
class LegacyVtkWriter;
class MshWriter;
}
} // namespace rw

namespace rw {

/*!
 * @brief Collection of methods and database related to writing
 *
 * This namespace provides methods and data members specific to writing of
 * the mesh data and simulation data. Currently, .vtu and .msh is supported.
 */
namespace writer {

/*!
 * @brief A interface class writing data
 *
 * This interface separates the caller from vtk library.
 */
class Writer {

public:
  /*!
   * @brief Constructor
   */
  Writer();

  /*!
   * @brief Constructor
   *
   * Creates and opens .vtu file of name given by filename. The file remains
   * open till the close() function is invoked or if the instance of this
   * class is destroyed.
   *
   * @param filename Name of file which will be created
   * @param format Format of the output file, e.g. "vtu", "msh"
   * @param compress_type Specify the compression type (optional)
   */
  Writer(const std::string &filename, const std::string &format, const std::string &compress_type);


  /*!
   * @brief Constructor
   *
   * Creates and opens .vtu file of name given by filename. The file remains
   * open till the close() function is invoked or if the instance of this
   * class is destroyed.
   *
   * @param format Format of the output file, e.g. "vtu", "msh"
   * @param compress_type Specify the compression type (optional)
   */
   Writer(const std::string &format, const std::string &compress_type);




  /*! @brief Destructor */
  ~Writer();

  /*!
   * @brief Function to trigger the output
   */
  void write(inp::Input *deck,data::DataManager* dataManager,size_t d_n, double d_time);


  /*!
   * @brief Open a .vtu file
   *
   * @param filename Name of file which will be created
   * @param format Format of the output file, e.g. "vtu", "msh"
   * @param compress_type Compression type (optional)
   */
  void open(const std::string &filename, const std::string &format = "vtu",
            const std::string &compress_type = "");

  /**
   * @name Mesh data
   */
  /**@{*/

  /*!
   * @brief Writes the nodes to the file
   * @param nodes Reference positions of the nodes
   * @param u Nodal displacements
   */
  void appendNodes(const std::vector<util::Point3> *nodes,
                   const std::vector<util::Point3> *u = nullptr);

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
                  const std::vector<util::Point3> *u = nullptr);

  /** @}*/

  /**
   * @name Point data
   */
  /**@{*/

  /*!
   * @brief Writes the scalar point data to the file
   * @param name Name of the data
   * @param data Vector containing the data
   */
  void appendPointData(const std::string &name,
                       const std::vector<uint8_t> *data);

  /*!
   * @brief Writes the scalar point data to the file
   * @param name Name of the data
   * @param data Vector containing the data
   */
  void appendPointData(const std::string &name,
                       const std::vector<size_t> *data);

  /*!
   * @brief Writes the scalar point data to the file
   * @param name Name of the data
   * @param data Vector containing the data
   */
  void appendPointData(const std::string &name, const std::vector<int> *data);

  /*!
   * @brief Writes the scalar point data to the file
   * @param name Name of the data
   * @param data Vector containing the data
   */
  void appendPointData(const std::string &name, const std::vector<float> *data);

  /*!
   * @brief Writes the scalar point data to the file
   * @param name Name of the data
   * @param data Vector containing the data
   */
  void appendPointData(const std::string &name,
                       const std::vector<double> *data);

  /*!
   * @brief Writes the vector point data to the file
   * @param name Name of the data
   * @param data Vector containing the data
   */
  void appendPointData(const std::string &name,
                       const std::vector<util::Point3> *data);

  /*!
 * @brief Writes the symmetric matrix data associated to nodes to the
 * file
 * @param name Name of the data
 * @param data Vector containing the data
 */
  void appendPointData(const std::string &name,
                       const std::vector<util::SymMatrix3> *data);

 /*!
 * @brief Writes the matrix data associated to nodes to the
 * file
 * @param name Name of the data
 * @param data Vector containing the data
 */
  void appendPointData(const std::string &name,
                       const std::vector<util::Matrix33> *data);


  /** @}*/

  /**
   * @name Cell data
   */
  /**@{*/

  /*!
   * @brief Writes the float data associated to cells to the file
   * @param name Name of the data
   * @param data Vector containing the data
   */
  void appendCellData(const std::string &name, const std::vector<float> *data);

  /*!
   * @brief Writes the symmetric matrix data associated to cells to the file
   * @param name Name of the data
   * @param data Vector containing the data
   */
  void appendCellData(const std::string &name,
                      const std::vector<util::SymMatrix3> *data);

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

  /*! @brief Pointer to the vtk writer class */
  rw::writer::LegacyVtkWriter *d_legacyVtkWriter_p;

  /*! @brief Pointer to the vtk writer class */
  rw::writer::MshWriter *d_mshWriter_p;

  /*! @brief Format of output file */
  std::string d_format;

  /*! @brief Format of output file */
  std::string d_compression_type;

}; // class Writer

} // namespace writer

} // namespace rw

#endif // RW_WRITER_H
