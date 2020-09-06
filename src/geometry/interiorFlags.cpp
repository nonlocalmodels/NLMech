////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "interiorFlags.h"
#include "inp/decks/interiorFlagsDeck.h"
#include "util/compare.h"
#include "util/utilGeom.h"
#include <iostream>
#include "util/utilIO.h"

//
// BaseInterior
//
geometry::BaseInterior::BaseInterior(
    inp::InteriorFlagsDeck *deck,
    std::pair<std::vector<double>, std::vector<double>> bbox,
    std::vector<std::pair<std::string, std::vector<double>>> no_fail_regions)
    : d_noFailTol(deck->d_noFailTol),
      d_bbox(std::move(bbox)),
      d_noFailRegions(no_fail_regions) {}

bool geometry::BaseInterior::getInteriorFlag(const size_t &i,
                                             const util::Point3 &x) {
  return true;
}

std::string geometry::BaseInterior::printStr(int nt, int lvl) const {
  auto tabS = util::io::getTabS(nt);
  std::ostringstream oss;
  oss << tabS << "------- BaseInterior --------" << std::endl << std::endl;
  oss << tabS << "Number of interior flags = " << d_intFlags.size()
      << std::endl;
  oss << tabS << "Bounding box = " << util::io::printBoxStr(d_bbox)
      << std::endl;
  oss << tabS << "Number of no-fail region = " << d_noFailRegions.size()
      << std::endl;
  oss << tabS << "No-fail region tol = " << d_noFailTol << std::endl;
  oss << tabS << std::endl;

  return oss.str();
}

//
// ComputeInterior
//
geometry::ComputeInterior::ComputeInterior(
    inp::InteriorFlagsDeck *deck,
    std::pair<std::vector<double>, std::vector<double>> bbox,
    std::vector<std::pair<std::string, std::vector<double>>> no_fail_regions)
    : geometry::BaseInterior(deck, std::move(bbox), no_fail_regions) {}

bool geometry::ComputeInterior::getInteriorFlag(const size_t &i,
                                                const util::Point3 &x) {
  // check for boundary of domain
  {
    // check if x coordinate is on left side of the left interior boundary
    if (util::compare::definitelyLessThan(x.d_x, d_bbox.first[0] + d_noFailTol))
      return false;

    // check if y coordinate is below the bottom interior boundary
    if (util::compare::definitelyLessThan(x.d_y, d_bbox.first[1] + d_noFailTol))
      return false;

    // check if x coordinate is on right side of the right interior
    // boundary
    if (util::compare::definitelyGreaterThan(x.d_x,
                                             d_bbox.second[0] - d_noFailTol))
      return false;

    // check if y coordinate is on top of the top interior boundary
    if (util::compare::definitelyGreaterThan(x.d_y,
                                             d_bbox.second[1] - d_noFailTol))
      return false;
  }

  // check for regions provided
  for (const auto &data : d_noFailRegions) {
    if (data.first == "Circle") {
      auto center =
          util::Point3(data.second[0], data.second[1], data.second[2]);
      auto r = data.second[3];

      if (util::compare::definitelyLessThan(center.dist(x), r))
        return false;  // want to not fail point x

    } else if (data.first == "Rectangle") {
      if (util::geometry::isPointInsideRectangle(x, data.second[0],
                                                 data.second[2], data.second[1],
                                                 data.second[3]))
        return false;  // want to not fail point x
    }
  }

  return true;
}

std::string geometry::ComputeInterior::printStr(int nt, int lvl) const {
  auto tabS = util::io::getTabS(nt);
  std::ostringstream oss;
  oss << tabS << "------- ComputeInterior --------" << std::endl << std::endl;
  oss << tabS << "Number of interior flags = " << d_intFlags.size()
      << std::endl;
  oss << tabS << "Bounding box = " << util::io::printBoxStr(d_bbox)
      << std::endl;
  oss << tabS << "Number of no-fail region = " << d_noFailRegions.size()
      << std::endl;
  oss << tabS << "No-fail region tol = " << d_noFailTol << std::endl;
  oss << tabS << std::endl;

  return oss.str();
}

//
// DataInterior
//
geometry::DataInterior::DataInterior(
    inp::InteriorFlagsDeck *deck, const std::vector<util::Point3> *nodes,
    std::pair<std::vector<double>, std::vector<double>> bbox,
    std::vector<std::pair<std::string, std::vector<double>>> no_fail_regions)
    : geometry::BaseInterior(deck, std::move(bbox), no_fail_regions) {
  size_t s = nodes->size() / 8;
  if (s * 8 < nodes->size()) s++;
  d_intFlags = std::vector<uint8_t>(s, uint8_t(0));
  for (size_t i = 0; i < nodes->size(); i++) {
    bool flag = true;
    auto x = (*nodes)[i];

    // check for boundary of domain
    {
      // check if x coordinate is on left side of the left interior boundary
      if (util::compare::definitelyLessThan(x.d_x,
                                            d_bbox.first[0] + d_noFailTol))
        flag = false;

      // check if y coordinate is below the bottom interior boundary
      if (util::compare::definitelyLessThan(x.d_y,
                                            d_bbox.first[1] + d_noFailTol))
        flag = false;

      // check if x coordinate is on right side of the right interior
      // boundary
      if (util::compare::definitelyGreaterThan(x.d_x,
                                               d_bbox.second[0] - d_noFailTol))
        flag = false;

      // check if y coordinate is on top of the top interior boundary
      if (util::compare::definitelyGreaterThan(x.d_y,
                                               d_bbox.second[1] - d_noFailTol))
        flag = false;
    }

    // check for regions provided
    for (const auto &data : d_noFailRegions) {
      if (data.first == "Circle") {
        auto center =
            util::Point3(data.second[0], data.second[1], data.second[2]);
        auto r = data.second[3];

        if (util::compare::definitelyLessThan(center.dist(x), r))
          flag = false;  // want to not fail point x

      } else if (data.first == "Rectangle") {
        if (util::geometry::isPointInsideRectangle(
                x, data.second[0], data.second[2], data.second[1],
                data.second[3]))
          flag = false;  // want to not fail point x
      }
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

std::string geometry::DataInterior::printStr(int nt, int lvl) const {
  auto tabS = util::io::getTabS(nt);
  std::ostringstream oss;
  oss << tabS << "------- DataInterior --------" << std::endl << std::endl;
  oss << tabS << "Number of interior flags = " << d_intFlags.size()
      << std::endl;
  oss << tabS << "Bounding box = " << util::io::printBoxStr(d_bbox)
      << std::endl;
  oss << tabS << "Number of no-fail region = " << d_noFailRegions.size()
      << std::endl;
  oss << tabS << "No-fail region tol = " << d_noFailTol << std::endl;
  oss << tabS << std::endl;

  return oss.str();
}

//
// InteriorFlags
//

geometry::InteriorFlags::InteriorFlags(
    inp::InteriorFlagsDeck *deck, const std::vector<util::Point3> *nodes,
    const std::pair<std::vector<double>, std::vector<double>> &bbox) {
  if (deck->d_noFailActive) {
    if (deck->d_computeAndNotStoreFlag)
      d_interior_p =
          new geometry::ComputeInterior(deck, bbox, deck->d_noFailRegions);
    else
      d_interior_p =
          new geometry::DataInterior(deck, nodes, bbox, deck->d_noFailRegions);
  } else
    d_interior_p =
        new geometry::BaseInterior(deck, bbox, deck->d_noFailRegions);
}

bool geometry::InteriorFlags::getInteriorFlag(const size_t &i,
                                              const util::Point3 &x) {
  return d_interior_p->getInteriorFlag(i, x);
}

std::string geometry::InteriorFlags::printStr(int nt, int lvl) const {
  auto tabS = util::io::getTabS(nt);
  std::ostringstream oss;
  oss << tabS << "------- InteriorFlags --------" << std::endl << std::endl;
  oss << tabS << "Interior data" << std::endl;
  oss << d_interior_p->printStr(nt + 1, lvl);
  oss << tabS << std::endl;

  return oss.str();
}
