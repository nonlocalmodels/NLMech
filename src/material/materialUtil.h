////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef MATERIAL_UTIL_H
#define MATERIAL_UTIL_H

#include "util/point.h"
#include "pdMaterial.h"
#include "geometry/fracture.h"
#include <limits>
#include <string>
#include <vector>

namespace material {

/*!
   * @brief Computes the moment \f$ m_x \f$ term in state-based peridynamic
   * formulation
   *
   * @param nodes Reference configuration of nodes
   * @param nodal_vol Nodal volumes
   * @param mesh_size Mesh size
   * @param material Material which provides influence function value,
   * horizon, etc
   * @param mx Data which will hold the computed values
   * @param compute_in_parallel Whether to compute in parallel
   */
void computeStateMx(const std::vector<util::Point3> &nodes,
                    const std::vector<double> &nodal_vol,
                    const std::vector<std::vector<size_t>> &neighbor_list,
                    const double &mesh_size,
                    material::pd::Material *material,
                    std::vector<double> &mx, bool compute_in_parallel = false);

/*!
   * @brief Computes the moment \f$ \theta_x \f$ term in state-based
   * peridynamic formulation
   *
   * @param nodes Reference configuration of nodes
   * @param nodes_disp Current displacement of nodes
   * @param nodal_vol Nodal volumes
   * @param material Material which provides influence function value,
   * horizon, etc
   * @param mx Moment mx data
   * @param thetax Data which will hold the computed values
   * @param compute_in_parallel Whether to compute in parallel
   */
void computeStateThetax(const std::vector<util::Point3> &nodes,
                    const std::vector<util::Point3> &nodes_disp,
                    const std::vector<double> &nodal_vol,
                        const std::vector<std::vector<size_t>> &neighbor_list,
                        const double &mesh_size,
                    material::pd::Material *material,
                        const geometry::Fracture *fracture,
                    const std::vector<double> &mx,
                    std::vector<double> &thetax, bool compute_in_parallel = false);

/*!
   * @brief Compute hydrostatic strain for rob's state-based model
   *
   * @param nodes Reference configuration of nodes
   * @param nodes_disp Current displacement of nodes
   * @param nodal_vol Nodal volumes
   * @param material Material which provides influence function value,
   * horizon, etc
   * @param mx Moment mx data
   * @param thetax Data which will hold the computed values
   * @param dim Dimension
   * @param compute_in_parallel Whether to compute in parallel
   */
void computeHydrostaticStrain(const std::vector<util::Point3> &nodes,
                               const std::vector<util::Point3> &nodes_disp,
                               const std::vector<double> &nodal_vol,
                              const std::vector<std::vector<size_t>> &neighbor_list,
                               const double &mesh_size,
                               material::pd::Material *material,
                               const geometry::Fracture *fracture,
                               std::vector<double> &thetax,
                               size_t dim, bool compute_in_parallel = false);

/*!
   * @brief Updates bond's fracture state, i.e. broken or not brokern
   *
   * @param nodes Reference configuration of nodes
   * @param nodes_disp Current displacement of nodes
   * @param material Material which provides influence function value,
   * horizon, etc
   * @param compute_in_parallel Whether to compute in parallel
   */
void updateBondFractureData(const std::vector<util::Point3> &nodes,
                              const std::vector<util::Point3> &nodes_disp,
                            const std::vector<std::vector<size_t>> &neighbor_list,
                              material::pd::Material *material,
                              geometry::Fracture *fracture, bool
                              compute_in_parallel = false);

} // namespace material

#endif // MATERIAL_UTIL_H
