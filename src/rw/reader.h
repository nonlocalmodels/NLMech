// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef RW_READER_H
#define RW_READER_H

#include "../util/point.h"
#include <vector>

//! Collection of methods and database related to reading and writing mesh data
namespace rw {

/*! @brief Methods and database used in reading operation */
namespace reader {

/*!
 * @brief Reads mesh data into node file and element file
 * @param filename name of mesh file
 * @param dim dimension
 * @param nodes vector of nodes data
 * @param volumes vector holding volume of the nodes
 */
void readCsvFile(const std::string &filename, size_t dim,
                 std::vector<util::Point3> *nodes,
                 std::vector<double> *volumes);

/*!
 * @brief Reads mesh data into node file and element file
 * @param filename name of mesh file
 * @param dim dimension
 * @param nodes vector of nodes data
 * @param element_type type of element
 * @param enc vector holding element-node connectivity
 * @param nec vector holding node-element connectivity
 * @param volumes vector holding volume of the nodes
 * @param is_fd flag indicating if this mesh is for finite_difference
 * simulation
 */
void readVtuFile(const std::string &filename, size_t dim,
                 std::vector<util::Point3> *nodes, size_t &element_type,
                 std::vector<size_t> *enc,
                 std::vector<std::vector<size_t>> *nec,
                 std::vector<double> *volumes, bool is_fd = false);

/*!
 * @brief Reads mesh data into node file and element file
 * @param filename name of mesh file
 * @param dim dimension
 * @param nodes vector of nodes data
 * @param element_type type of element
 * @param enc vector holding element-node connectivity
 * @param nec vector holding node-element connectivity
 * @param volumes vector holding volume of the nodes
 * @param is_fd flag indicating if this mesh is for finite_difference
 * simulation
 */
void readMshFile(const std::string &filename, size_t dim,
                 std::vector<util::Point3> *nodes, size_t &element_type,
                 std::vector<size_t> *enc,
                 std::vector<std::vector<size_t>> *nec,
                 std::vector<double> *volumes, bool is_fd = false);

} // namespace reader

} // namespace rw

#endif // READER_H
