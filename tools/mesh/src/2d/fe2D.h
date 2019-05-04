// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef TOOLS_MESH_TWOD_H
#define TOOLS_MESH_TWOD_H

#include <string>

namespace tools {
namespace mesh {

/*!
 * @brief Generates uniform mesh in 2-d and writes the data to files
 * specified in input filename
 * @param filename Name of YAML input file
 */
void fe2D(const std::string &filename);

/*!
 * @brief Generates uniform triangular mesh in 2-d and writes the data to files
 * specified in input filename
 * @param filename Name of YAML input file
 */
void uniformTri(const std::string &filename);

/*!
 * @brief Generates uniform square mesh in 2-d and writes the data to files
 * specified in input filename
 * @param filename Name of YAML input file
 */
void uniformSquare(const std::string &filename);

/*!
 * @brief Generates uniform symmetric triangular mesh in 2-d and writes the
 * data to files specified in input filename
 * @param filename Name of YAML input file
 */
void uniformTriSym(const std::string &filename);

} // namespace mesh

} // namespace tools

#endif // TOOLS_MESH_TWOD_H