////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "interiorFlags.h"

#include <iostream>

#include "inp/decks/interiorFlagsDeck.h"
#include "util/compare.h"

//
// BaseInterior
//
geometry::BaseInterior::BaseInterior(
    inp::InteriorFlagsDeck *deck, const size_t dim,
    std::pair<std::vector<double>, std::vector<double>> bbox)
    : d_noFailTol(deck->d_noFailTol), d_dim(dim), d_bbox(std::move(bbox)) {}

bool geometry::BaseInterior::getInteriorFlag(const size_t &i,
                                             const util::Point3 &x) {
  return true;
}

//
// ComputeInterior
//
geometry::ComputeInterior::ComputeInterior(
    inp::InteriorFlagsDeck *deck, const size_t dim,
    std::pair<std::vector<double>, std::vector<double>> bbox)
    : geometry::BaseInterior(deck, dim, std::move(bbox)) {}

bool geometry::ComputeInterior::getInteriorFlag(const size_t &i,
                                                const util::Point3 &x) {
  // check if x coordinate is on left side of the left interior boundary
  if (util::compare::definitelyLessThan(x.d_x, d_bbox.first[0] + d_noFailTol))
    return false;

  // check if x coordinate is on right side of the right interior
  // boundary
  if (util::compare::definitelyGreaterThan(x.d_x,
                                           d_bbox.second[0] - d_noFailTol))
    return false;

  if (d_dim > 1) {
    // check if y coordinate is below the bottom interior boundary
    if (util::compare::definitelyLessThan(x.d_y, d_bbox.first[1] + d_noFailTol))
      return false;

    // check if y coordinate is on top of the top interior boundary
    if (util::compare::definitelyGreaterThan(x.d_y,
                                             d_bbox.second[1] - d_noFailTol))
      return false;
  }

  if (d_dim > 2) {
    // check if z coordinate is below the bottom interior boundary
    if (util::compare::definitelyLessThan(x.d_z, d_bbox.first[2] + d_noFailTol))
      return false;

    // check if z coordinate is on top of the top interior boundary
    if (util::compare::definitelyGreaterThan(x.d_z,
                                             d_bbox.second[2] - d_noFailTol))
      return false;
  }

  return true;
}

//
// DataInterior
//
geometry::DataInterior::DataInterior(
    inp::InteriorFlagsDeck *deck, const size_t dim,
    const std::vector<util::Point3> *nodes,
    std::pair<std::vector<double>, std::vector<double>> bbox)
    : geometry::BaseInterior(deck, dim, std::move(bbox)) {
  size_t s = nodes->size() / 8;
  if (s * 8 < nodes->size()) s++;
  d_intFlags = std::vector<uint8_t>(s, uint8_t(0));

  for (size_t i = 0; i < nodes->size(); i++) {
    bool flag = true;
    auto x = (*nodes)[i];

    // check if x coordinate is on left side of the left interior boundary
    flag = !util::compare::definitelyLessThan(x.d_x,
                                              d_bbox.first[0] + d_noFailTol);

    // check if x coordinate is on right side of the right interior
    // boundary
    flag = !util::compare::definitelyGreaterThan(
        x.d_x, d_bbox.second[0] - d_noFailTol);

    if (d_dim > 1) {
      // check if y coordinate is below the bottom interior boundary
      flag = !util::compare::definitelyLessThan(x.d_y,
                                                d_bbox.first[1] + d_noFailTol);

      // check if y coordinate is on top of the top interior boundary
      flag = !util::compare::definitelyGreaterThan(
          x.d_y, d_bbox.second[1] - d_noFailTol);
    }

    if (d_dim > 2) {
      // check if z coordinate is below the bottom interior boundary
      flag = !util::compare::definitelyLessThan(x.d_z,
                                                d_bbox.first[2] + d_noFailTol);

      // check if z coordinate is on top of the top interior boundary
      flag = !util::compare::definitelyGreaterThan(
          x.d_z, d_bbox.second[2] - d_noFailTol);
    }

    // set loc_bit = i%8 of loc_i = i/8 to flag

    // to set i^th bit as true of integer a,
    // a |= 1UL << (i % 8)

    // to set i^th bit as flase of integer a,
    // a &= ~(1UL << (i % 8))
    flag ? (d_intFlags[i / 8] |= 1UL << (i % 8))
         : (d_intFlags[i / 8] &= ~(1UL << (i % 8)));
  }
}

bool geometry::DataInterior::getInteriorFlag(const size_t &i,
                                             const util::Point3 &x) {
  auto bond = d_intFlags[i / 8];
  return bond >> (i % 8) & 1UL;
}

geometry::InteriorFlags::InteriorFlags(
    inp::InteriorFlagsDeck *deck, const size_t dim,
    const std::vector<util::Point3> *nodes,
    const std::pair<std::vector<double>, std::vector<double>> &bbox) {
  if (deck->d_noFailActive) {
    if (deck->d_computeAndNotStoreFlag)
      d_interior_p = new geometry::ComputeInterior(deck, dim, bbox);
    else
      d_interior_p = new geometry::DataInterior(deck, dim, nodes, bbox);
  } else
    d_interior_p = new geometry::BaseInterior(deck, dim, bbox);
}

bool geometry::InteriorFlags::getInteriorFlag(const size_t &i,
                                              const util::Point3 &x) {
  return d_interior_p->getInteriorFlag(i, x);
}
