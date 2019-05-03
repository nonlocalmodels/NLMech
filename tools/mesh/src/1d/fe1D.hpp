// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef TOOLS_MESH_ONED_H
#define TOOLS_MESH_ONED_H

#include <string>

/*!
 * @brief Namespace consisting of libraries which are used either in
 * pre-process or in post-process.
 *
 * This consists of
 *
 * - simple mesh generators
 * - data comparison tool which compares the displacement field of two mesh
 * sizes and returns the \f$ L^2\f$ and \f$ sup \f$ norm of the error
 * - post-processing tool which computes the post-processing fields such as
 * stress, strain, and manipulates the result of the code for better plots
 */
namespace tools {

/*!
 * @brief Namespace for simple mesh generation
 *
 * In this namespace we define the methods for generating simple finite
 * element/finite difference mesh in 1-d and 2-d. For finite element, we can
 * specify three types of triangular mesh which differ in how we create
 * triangle on uniform grid.
 */
namespace mesh {

/*!
 * @brief Generate uniform mesh in 1-d and writes the data to files specified
 * in input filename
 * @param filename Name of YAML input file
 */
void fe1D(const std::string &filename);

} // namespace mesh

} // namespace tools

#endif // TOOLS_MESH_ONED_H