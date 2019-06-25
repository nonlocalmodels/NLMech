// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef RW_READER_H
#define RW_READER_H

#include "util/point.h"           // definition of Point3
#include <vector>

/*!
 * @brief Collection of methods and database related to reading and writing
 *
 * This namespace provides methods and data members specific to reading and
 * writing of the mesh data and simulation data.
 *
 * @sa VtkReader, MshReader
 */
namespace rw {

/*!
 * @brief Collection of methods and database related to reading
 *
 * This namespace provides methods and data members specific to reading of
 * the mesh data. Currently, .csv, .vtu and .msh files are supported.
 *
 * @sa VtkReader, MshReader
 */
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
 * @param num_elem number of elements
 * @param enc vector holding element-node connectivity
 * @param nec vector holding node-element connectivity
 * @param volumes vector holding volume of the nodes
 * @param is_fd flag indicating if this mesh is for finite_difference
 * simulation
 */
void readVtuFile(const std::string &filename, size_t dim,
                 std::vector<util::Point3> *nodes, size_t &element_type,
                 size_t &num_elem, std::vector<size_t> *enc,
                 std::vector<std::vector<size_t>> *nec,
                 std::vector<double> *volumes, bool is_fd = false);

/*!
 * @brief Reads mesh data into node file and element file
 * @param filename name of mesh file
 * @param dim dimension
 * @param nodes vector of nodes data
 * @param element_type type of element
 * @param num_elem number of elements
 * @param enc vector holding element-node connectivity
 * @param nec vector holding node-element connectivity
 * @param volumes vector holding volume of the nodes
 * @param is_fd flag indicating if this mesh is for finite_difference
 * simulation
 */
void readMshFile(const std::string &filename, size_t dim,
                 std::vector<util::Point3> *nodes, size_t &element_type,
                 size_t &num_elem, std::vector<size_t> *enc,
                 std::vector<std::vector<size_t>> *nec,
                 std::vector<double> *volumes, bool is_fd = false);

/*!
 * @brief Reads mesh data into node file and element file
 * @param filename name of mesh file
 * @param u Pointer to vector of nodal displacement
 * @param v Pointer to vector of nodal velocity
 * @param X Pointer to vector of nodal reference position (Optional)
 */
void readVtuFileRestart(const std::string &filename,
                        std::vector<util::Point3> *u,
                        std::vector<util::Point3> *v,
                        const std::vector<util::Point3> *X = nullptr);

} // namespace reader

} // namespace rw

#endif // RW_READER_H
