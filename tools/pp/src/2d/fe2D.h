// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef TOOLS_PP_TWOD_H
#define TOOLS_PP_TWOD_H

#include <string>

namespace tools {

/*!
 * @brief Namespace for postprocessing of simulation results
 *
 * In this namespace we define the methods to compute quantities of interest
 * from the simulation results. Examples are: strain and stress, scaling of
 * displacement for fracture plot, symmetrization of plots, etc.
 */
namespace pp {

/*!
 * @brief Processes 2-d simulation results and computes postprocessing
 * quantities
 * @param filename Name of YAML input file
 */
void fe2D(const std::string &filename);

} // namespace pp

} // namespace tools

#endif // TOOLS_PP_TWOD_H