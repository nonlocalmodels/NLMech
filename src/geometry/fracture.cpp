////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "fracture.h"
#include "inp/decks/fractureDeck.h"
#include "util/utilGeom.h"
#include <hpx/include/parallel_algorithm.hpp>

geometry::Fracture::Fracture(inp::FractureDeck *deck)
    : d_fractureDeck_p(deck) {}

geometry::Fracture::Fracture(
    inp::FractureDeck *deck, const std::vector<util::Point3> *nodes,
    const std::vector<std::vector<size_t>> *neighbor_list)
    : d_fractureDeck_p(deck) {
  d_fracture.resize(neighbor_list->size());

  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      neighbor_list->size(), [this, nodes, neighbor_list](boost::uint64_t i) {
        auto ns = (*neighbor_list)[i];

        size_t s = ns.size() / 8;
        if (s * 8 < ns.size()) s++;
        d_fracture[i] = std::vector<uint8_t>(s, uint8_t(0));

        for (auto &crack : d_fractureDeck_p->d_cracks)
          if (crack.d_activationTime < 0.) {
            this->computeFracturedBondFd(i, &crack, nodes, &ns);
            crack.d_crackAcrivated = true;
          }
      });  // end of parallel for loop

  f.get();
}

bool geometry::Fracture::addCrack(
    const double &time, const std::vector<util::Point3> *nodes,
    const std::vector<std::vector<size_t>> *neighbor_list) {
  for (auto &crack : d_fractureDeck_p->d_cracks) {
    if (!crack.d_crackAcrivated) {
      if (util::compare::definitelyLessThan(crack.d_activationTime, time)) {
        std::cout << "Fracture: Adding crack to system\n";

        auto f = hpx::parallel::for_loop(
            hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
            neighbor_list->size(),
            [this, nodes, neighbor_list, &crack](boost::uint64_t i) {
              auto ns = (*neighbor_list)[i];

              this->computeFracturedBondFd(i, &crack, nodes, &ns);
            });  // end of parallel for loop

        f.get();

        crack.d_crackAcrivated = true;

        return true;
      }
    }
  }

  return false;
}

void geometry::Fracture::computeFracturedBondFd(
    const size_t &i, inp::EdgeCrack *crack,
    const std::vector<util::Point3> *nodes,
    const std::vector<size_t> *neighbors) {
  //
  //
  // Here [ ] represents a mesh node and o------o represents a crack.
  //
  //
  //                     pt = pr
  //
  //                      o
  //                     /
  //         [ ]-----[ ]/----[ ]
  //          |       |/      |
  //          |       /       |
  //          |      /|       |
  //         [ ]----/[ ]-----[ ]
  //          |    /  |       |
  //          |   /   |       |
  //          |  /    |       |
  //         [ ]/----[ ]-----[ ]
  //           /
  //          o
  //
  //     pb = pl

  //
  //
  // By design, the crack is offset a very small amount (5.0E-8) to bottom
  // and to right side.
  //
  util::Point3 i_node = (*nodes)[i];

  // we assume pb is below pt i.e. pb.y < pt.y
  util::Point3 pb = crack->d_pb;
  util::Point3 pt = crack->d_pt;

  // check if point is outside crack line
  if (crack->ptOutside(i_node, crack->d_o, pb, pt)) return;

  // find if this node is on right side or left side of the crack line
  bool left_side = crack->ptLeftside(i_node, pb, pt);

  //
  if (left_side) {
    // loop over neighboring nodes
    for (size_t j = 0; j < neighbors->size(); j++) {
      size_t id_j = (*neighbors)[j];
      util::Point3 j_node = (*nodes)[id_j];

      // check if j_node lies on right side of crack line
      bool modify = true;

      // check if point lies outside crack line
      if (crack->ptOutside(j_node, crack->d_o, pb, pt)) modify = false;

      // modify only those nodes which lie on opposite side (in this
      // case on right side)
      if (crack->ptRightside(j_node, pb, pt) and modify)
        this->setBondState(i, j, modify);
    }

  }  // left side
  else {
    // loop over neighboring nodes
    for (size_t j = 0; j < neighbors->size(); j++) {
      size_t id_j = (*neighbors)[j];
      util::Point3 j_node = (*nodes)[id_j];

      // check if j_node lies on left side of crack line
      // As pointed out in the beginning, since we are looking for
      // nodes on left side of the crack, we need to be prepared to
      // handle the case when node is very close to crack line
      // i.e. area (pb, pt, j_node) >= 0.

      auto modify = true;

      if (crack->ptOutside(j_node, crack->d_o, pb, pt)) modify = false;

      // modify only those nodes which lie on opposite side (in this
      // case on left side)
      if (crack->ptLeftside(j_node, pb, pt) and modify)
        this->setBondState(i, j, modify);
    }
  }  // right side
}  // computeFracturedBondFd

void geometry::Fracture::setBondState(const size_t &i, const size_t &j,
                                      const bool &state) {
  // to set i^th bit as true of integer a,
  // a |= 1UL << (i % 8)

  // to set i^th bit as false of integer a,
  // a &= ~(1UL << (i % 8))

  state ? (d_fracture[i][j / 8] |= 1UL << (j % 8))
        : (d_fracture[i][j / 8] &= ~(1UL << (j % 8)));
}

bool geometry::Fracture::getBondState(const size_t &i, const size_t &j) const {
  //  if (d_fracture.size() <= i) {
  //    std::cerr << "Error: Size of fracture data is invalid\n";
  //    exit(1);
  //  }
  //
  //  if (d_fracture[i].size() <= j / 8) {
  //    std::cerr << "Error: Size of fracture data at node " << i << " is
  //    invalid\n"; exit(1);
  //  }

  auto bond = d_fracture[i][j / 8];
  return bond >> (j % 8) & 1UL;
}

const std::vector<uint8_t> geometry::Fracture::getBonds(const size_t &i) const {
  return d_fracture[i];
}
