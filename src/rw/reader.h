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
 */
namespace rw {

/*!
 * @brief Collection of methods and database related to reading
 *
 * This namespace provides methods and data members specific to reading of
 * the mesh data. Currently, .csv, .vtu and .msh files are supported.
 */
namespace reader {

/*!
 * @brief Reads mesh data into node file and element file
 * @param filename Name of mesh file
 * @param dim Dimension
 * @param nodes Vector of nodes data
 * @param volumes Vector holding volume of the nodes
 */
void readCsvFile(const std::string &filename, size_t dim,
                 std::vector<util::Point3> *nodes,
                 std::vector<double> *volumes);

/*!
 * @brief Reads mesh data into node file and element file
 * @param filename Name of mesh file
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
void readVtuFile(const std::string &filename, size_t dim,
                 std::vector<util::Point3> *nodes, size_t &element_type,
                 size_t &num_elem, std::vector<size_t> *enc,
                 std::vector<std::vector<size_t>> *nec,
                 std::vector<double> *volumes, bool is_fd = false);

/*!
 * @brief Reads mesh data into node file and element file
 * @param filename Name of mesh file
 * @param u Pointer to vector of nodal displacement
 * @param v Pointer to vector of nodal velocity
 * @param X Pointer to vector of nodal reference position (Optional)
 */
void readVtuFileRestart(const std::string &filename,
                        std::vector<util::Point3> *u,
                        std::vector<util::Point3> *v,
                        const std::vector<util::Point3> *X = nullptr);

/*!
 * @brief Reads data of specified tag from the vtu file
 * @param filename Name of mesh file
 * @param tag Name of point data to be read from .vtu file
 * @param data Pointer to vector of point data
 * @return bool True if found the data in file
 */
bool readVtuFilePointData(const std::string &filename,
                          const std::string &tag,
                        std::vector<double> *data);

/*!
 * @brief Reads mesh data into node file and element file
 * @param filename Name of mesh file
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
void readMshFile(const std::string &filename, size_t dim,
                 std::vector<util::Point3> *nodes, size_t &element_type,
                 size_t &num_elem, std::vector<size_t> *enc,
                 std::vector<std::vector<size_t>> *nec,
                 std::vector<double> *volumes, bool is_fd = false);

/*!
 * @brief Reads mesh data into node file and element file
 * @param filename Name of mesh file
 * @param u Pointer to vector of nodal displacement
 * @param v Pointer to vector of nodal velocity
 * @param X Pointer to vector of nodal reference position (Optional)
 */
void readMshFileRestart(const std::string &filename,
                        std::vector<util::Point3> *u,
                        std::vector<util::Point3> *v,
                        const std::vector<util::Point3> *X = nullptr);

/*!
 * @brief Reads data of specified tag from the vtu file
 * @param filename Name of mesh file
 * @param tag Name of point data to be read from .vtu file
 * @param data Pointer to vector of point data
 * @return bool True if found the data in file
 */
bool readMshFilePointData(const std::string &filename,
                          const std::string &tag,
                          std::vector<double> *data);

} // namespace reader

} // namespace rw

#endif // RW_READER_H
