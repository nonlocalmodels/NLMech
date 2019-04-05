// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "fracture.h"
#include "../inp/decks/fractureDeck.h"

#include <hpx/include/parallel_algorithm.hpp>

geometry::Fracture::Fracture(
    inp::FractureDeck *deck, const std::vector<util::Point3> *nodes,
    const std::vector<std::vector<size_t>> *neighbor_list) {

  d_fractureDeck_p = deck;

  d_fracture.resize(neighbor_list->size());

  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      neighbor_list->size(), [this,nodes,neighbor_list](boost::uint64_t i) {
        std::vector<size_t> ineighs = (*neighbor_list)[i];

        size_t s = ineighs.size() / 8;
        if (s * 8 < ineighs.size())
          s++;
        d_fracture[i] = std::vector<uint8_t>(s, uint8_t(0));

        for (auto crack : d_fractureDeck_p->d_cracks) {
          if (crack.d_o == -1)
            this->computeFracturedBondFdAlongY(i, &crack, nodes, &ineighs);
        } // loop over cracks
      }); // end of parallel for loop

  f.get();
}

void geometry::Fracture::computeFracturedBondFdAlongY(const size_t &i,
    inp::EdgeCrack *crack, const std::vector<util::Point3> *nodes,
    const std::vector<size_t> *neighbors) {

  //
  //
  // Here [ ] represents a mesh node and o------o represents a crack.
  //
  //              pt
  //              o
  //              |
  //             [|]---[ ]
  //              |
  //              |
  //             [|]---[ ]
  //              |
  //              |
  //             [|]---[ ]
  //              |
  //              o
  //             pb
  //
  //
  // By design, the crack is offset a very small amount (5.0E-8) to right
  // side. Which means few mesh nodes which are on left side of crack can get
  // very close to the crack line. Whereas the nodes which are on right can
  // not get too close to crack line.
  //
  // Thus, when dealing with left side, we need to be prepared to
  // handle the case when nodes on left side are very close to the line.
  //
  util::Point3 i_node = (*nodes)[i];

  // we assume pb is below pt i.e. pb.y < pt.y
  util::Point3 pb = util::Point3(crack->d_pb[0], crack->d_pb[1], 0.);
  util::Point3 pt = util::Point3(crack->d_pt[0], crack->d_pt[1], 0.);
  int o = -1;

  // if (compare::definitelyLessThan(i_node.y, pb.y) or compare::definitelyGreaterThan(i_node.y, pt.y))
  // 	return;

  // check if point is outside crack line
  if (crack->ptOutside(i_node, o, pb, pt))
    return;

  // find if this node is on right side or left side of the crack line
  bool left_side = crack->ptLeftside(i_node, o, pb, pt);

  //
  if (left_side) {

    // loop over neighboring nodes
    for (size_t j = 0; j < neighbors->size(); j++) {

      size_t id_j = (*neighbors)[j];
      util::Point3 j_node = (*nodes)[id_j];

      // check if j_node lies on right side of crack line
      bool modify = true;

      // check if point lies outside crack line
      if (crack->ptOutside(j_node, o, pb, pt))
        modify = false;

      // modify only those nodes which lie on opposite side (in this
      // case on right side)
      if (crack->ptRightside(j_node, o, pb, pt) and modify)
        this->setBondState(i,j, modify);
    }

  } // left side
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

      if (crack->ptOutside(j_node, o, pb, pt))
        modify = false;

      // modify only those nodes which lie on opposite side (in this
      // case on left side)
      if (crack->ptLeftside(j_node, o, pb, pt) and modify)
        this->setBondState(i,j, modify);
    }
  } // right side

}

void geometry::Fracture::markBondBroken(const size_t &i, const size_t &j) {

  d_fracture[i][j/8] |= 1UL << (j%8);
}

void geometry::Fracture::setBondState(const size_t &i, const size_t &j, bool
state) {

  if (state) {
    d_fracture[i][j / 8] |= 1UL << (j % 8);
  } else {
    d_fracture[i][j / 8] &= ~(1UL << (j % 8));
  }
}

bool geometry::Fracture::getBondState(const size_t &i, const size_t &j) {

  uint8_t bond = d_fracture[i][j/8];
  return (bond >> (j % 8)) & 1UL;
}

const std::vector<uint8_t> geometry::Fracture::getBonds(const size_t &i) {
  return d_fracture[i];
}
