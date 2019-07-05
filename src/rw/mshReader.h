// Copyright (c)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

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
   * @param filename Name of the mesh file
   */
  explicit MshReader(const std::string &filename);

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

private:
  /*! @brief filename */
  const std::string d_filename;
};

} // namespace reader

} // namespace rw

#endif // RW_MSHREADER_H
