////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef RW_MSHREADER_H
#define RW_MSHREADER_H

#include "util/point.h"           // definition of Point3
#include <string>
#include <vector>

namespace rw {

namespace reader {

/*!
 * @brief A class to read Gmsh (msh) mesh files
 */
class MshReader {

public:
  /*!
   * @brief Constructor
   * @param filename name of mesh file
   */
  explicit MshReader(const std::string &filename);

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

private:
  /**
   * @name Internal data
   */
  /**@{*/

  /*! @brief filename */
  const std::string d_filename;

  /** @}*/
};

} // namespace reader

} // namespace rw

#endif // RW_MSHREADER_H
