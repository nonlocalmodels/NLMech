// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "fracture.h"
#include "inp/decks/fractureDeck.h"
#include "util/utilGeom.h"
#include <hpx/include/parallel_algorithm.hpp>

geometry::Fracture::Fracture(inp::FractureDeck *deck)
    : d_fractureDeck_p(deck), d_crackOutData(geometry::CrackOutData()) {}

geometry::Fracture::Fracture(
    inp::FractureDeck *deck, const std::vector<util::Point3> *nodes,
    const std::vector<std::vector<size_t>> *neighbor_list)
    : d_fractureDeck_p(deck), d_crackOutData(geometry::CrackOutData()) {

  d_fracture.resize(neighbor_list->size());

  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      neighbor_list->size(), [this, nodes, neighbor_list](boost::uint64_t i) {
        auto ns = (*neighbor_list)[i];

        size_t s = ns.size() / 8;
        if (s * 8 < ns.size())
          s++;
        d_fracture[i] = std::vector<uint8_t>(s, uint8_t(0));

        for (auto crack : d_fractureDeck_p->d_cracks)
          this->computeFracturedBondFd(i, &crack, nodes, &ns);
      }); // end of parallel for loop

  f.get();
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
  if (crack->ptOutside(i_node, crack->d_o, pb, pt))
    return;

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
      if (crack->ptOutside(j_node, crack->d_o, pb, pt))
        modify = false;

      // modify only those nodes which lie on opposite side (in this
      // case on right side)
      if (crack->ptRightside(j_node, pb, pt) and modify)
        this->setBondState(i, j, modify);
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

      if (crack->ptOutside(j_node, crack->d_o, pb, pt))
        modify = false;

      // modify only those nodes which lie on opposite side (in this
      // case on left side)
      if (crack->ptLeftside(j_node, pb, pt) and modify)
        this->setBondState(i, j, modify);
    }
  } // right side
} // computeFracturedBondFd

void geometry::Fracture::setBondState(const size_t &i, const size_t &j,
                                      const bool &state) {

  // to set i^th bit as true of integer a,
  // a |= 1UL << (i % 8)

  // to set i^th bit as false of integer a,
  // a &= ~(1UL << (i % 8))

  state ? (d_fracture[i][j / 8] |= 1UL << (j % 8))
        : (d_fracture[i][j / 8] &= ~(1UL << (j % 8)));
}

void geometry::Fracture::setUpdateCrack(const size_t &dt) {
  d_fractureDeck_p->d_dtCrackOut = dt;
}

bool geometry::Fracture::getBondState(const size_t &i, const size_t &j) {

  auto bond = d_fracture[i][j / 8];
  return bond >> (j % 8) & 1UL;
}

const std::vector<uint8_t> geometry::Fracture::getBonds(const size_t &i) {
  return d_fracture[i];
}

size_t geometry::Fracture::getDtCrackOut() {
  return d_fractureDeck_p->d_dtCrackOut;
}

void geometry::Fracture::updateCrackAndOutput(
    const size_t &n, const double &time, const std::string &output_path,
    const double &horizon, const std::vector<util::Point3> *nodes,
    std::vector<util::Point3> *u, std::vector<float> *Z) {
  // Idea: Suppose we output crack tip and crack velocity data every dt
  // interval of time step and suppose we use dt_v time step interval to
  // compute the velocity. The velocity is simply given by
  //
  //  v_tip(t) = (x_tip(t) - x_tip(t-dt_v)) / dt_v
  //
  // Suppose k = N*dt for some integer N. Then we first need to compute the
  // crack tip at k - dt  and then compute the crack tip at k.
  //
  // Special case: When dt_v == dt then we use the crack tip computed at
  // previous call to this function, i.e. at k - dt = (N-1) * dt as old crack
  // tip and compute new crack tip at current time step k.
  //
  if (d_fractureDeck_p->d_dtCrackOut == 0)
    return;

  if ((d_fractureDeck_p->d_dtCrackVelocity < d_fractureDeck_p->d_dtCrackOut) &&
      ((n + d_fractureDeck_p->d_dtCrackVelocity) %
           d_fractureDeck_p->d_dtCrackOut ==
       0))
    updateCrack(n, time, horizon, nodes, u, Z);

  if (n % d_fractureDeck_p->d_dtCrackOut == 0) {
    updateCrack(n, time, horizon, nodes, u, Z);
    output(n, time, output_path, nodes, u, Z);
  }
}

void geometry::Fracture::updateCrack(const size_t &n, const double &time,
                                     const double &horizon,
                                     const std::vector<util::Point3> *nodes,
                                     std::vector<util::Point3> *u,
                                     std::vector<float> *Z) {

  // loop over crack lines
  size_t count = 0;
  for (auto &crack : d_fractureDeck_p->d_cracks) {

    if (crack.d_o == 0)
      continue;

    auto pb = crack.d_pb;
    auto pt = crack.d_pt;

    // search length
    double search_length = 1000.;

    // create search rectangle containing crack
    double rect_t[4];
    double rect_b[4];
    if (crack.d_o == -1) {
      rect_t[0] = pt.d_x - 2. * horizon;
      rect_t[1] = pt.d_y;
      rect_t[2] = pt.d_x + 2. * horizon;
      rect_t[3] = pt.d_y + search_length;

      rect_b[0] = pb.d_x - 2. * horizon;
      rect_b[1] = pb.d_y - search_length;
      rect_b[2] = pb.d_x + 2. * horizon;
      rect_b[3] = pb.d_y;
    } else if (crack.d_o == 1) {
      rect_t[0] = pt.d_x;
      rect_t[1] = pt.d_y - 2. * horizon;
      rect_t[2] = pt.d_x + search_length;
      rect_t[3] = pt.d_y + 2. * horizon;

      rect_b[0] = pb.d_x - search_length;
      rect_b[1] = pb.d_y - 2. * horizon;
      rect_b[2] = pb.d_x;
      rect_b[3] = pb.d_y + 2. * horizon;
    }

    // find and store id of node at crack tip
    long ib = -1;
    long it = -1;
    float Zb = 1000.;
    float Zt = 1000.;

    for (size_t i = 0; i < nodes->size(); i++) {
      auto xi = (*nodes)[i];
      auto yi = xi + (*u)[i];
      auto damage = (*Z)[i];

      if (util::geometry::isPointInsideRectangle(yi, rect_t[0], rect_t[2],
                                                 rect_t[1], rect_t[3]) &&
          util::compare::definitelyGreaterThan(damage, 1.) &&
          util::compare::definitelyLessThan(damage, Zt)) {
        it = i;
        Zt = damage;
      }

      if (util::geometry::isPointInsideRectangle(yi, rect_b[0], rect_b[2],
                                                 rect_b[1], rect_b[3]) &&
          util::compare::definitelyGreaterThan(damage, 1.) &&
          util::compare::definitelyLessThan(damage, Zb)) {
        ib = i;
        Zb = damage;
      }
    } // loop over nodes

    if (it != -1) {
      crack.d_it = it;
      crack.d_pt = (*nodes)[it] + (*u)[it];
      auto diff = crack.d_pt - pt;
      auto delta_t = time - d_crackOutData.d_timet;
      crack.d_lt += diff.length();
      crack.d_l += diff.length();
      crack.d_vt = util::Point3(diff.d_x / delta_t, diff.d_y / delta_t,
                                diff.d_z / delta_t);
      d_crackOutData.d_timet = time;
    }

    if (ib != -1) {
      crack.d_ib = ib;
      crack.d_pb = (*nodes)[ib] + (*u)[ib];
      auto diff = crack.d_pb - pb;
      auto delta_t = time - d_crackOutData.d_timeb;
      crack.d_lb += diff.length();
      crack.d_l += diff.length();
      crack.d_vb = util::Point3(diff.d_x / delta_t, diff.d_y / delta_t,
                                diff.d_z / delta_t);
      d_crackOutData.d_timeb = time;
    }
  } // loop over cracks
}

void geometry::Fracture::output(const size_t &n, const double &time,
                                const std::string &output_path,
                                const std::vector<util::Point3> *nodes,
                                std::vector<util::Point3> *u,
                                std::vector<float> *Z) {

  // create new file for every 10000 calls to this function
  int up_bound = 10000;
  size_t mod = d_crackOutData.d_updateCount % up_bound;
  if (d_crackOutData.d_updateCount > 1 && mod == 0) {
    d_crackOutData.d_fileOutCount++;
    d_crackOutData.d_needNewFile = true;
  } else
    d_crackOutData.d_needNewFile = false;

  // if the total call reaches 5*10000, close the open file and return
  if (d_crackOutData.d_updateCount >= 5 * up_bound) {
    if (d_crackOutData.d_file) {
      fclose(d_crackOutData.d_file);
      std::cout << "Warning: Number of times crack data output requested "
                   "exceeds the upper limit 10000.\n";
    }
    return;
  }

  // write
  if (d_crackOutData.d_updateCount == 0 || d_crackOutData.d_needNewFile) {
    std::string filename =
        output_path + "/" + d_fractureDeck_p->d_crackOutFilename + "_" +
        std::to_string(d_crackOutData.d_fileOutCount) + ".csv";

    // before opening a new file, close the previous file
    if (d_crackOutData.d_needNewFile && d_crackOutData.d_file)
      fclose(d_crackOutData.d_file);
    d_crackOutData.d_file = fopen(filename.c_str(), "w");
    // write header
    fprintf(d_crackOutData.d_file,
            "crack_id, time, ib, pb.x, pb.y, Zb, it, pt.x,"
            " pt.y, Zt, l, lb,  lt, vb.x, vb.y, vb_mag,"
            " vt.x, vt.y, vt_mag\n");
  }

  for (size_t i = 0; i < d_fractureDeck_p->d_cracks.size(); i++) {
    auto crack = d_fractureDeck_p->d_cracks[i];
    fprintf(d_crackOutData.d_file,
            "%lu, %6.8e, %u, %6.8e, %6.8e, %6.8e, %u, %6.8e, "
            "%6.8e, %6.8e, %6.8e, %6.8e, %6.8e, %6.8e, %6.8e, %6.8e, %6.8e, "
            "%6.8e, %6.8e\n",
            i, time, crack.d_ib, crack.d_pb.d_x, crack.d_pb.d_y,
            (*Z)[crack.d_ib], crack.d_it, crack.d_pt.d_x, crack.d_pt.d_y,
            (*Z)[crack.d_it], crack.d_l, crack.d_lb, crack.d_lt, crack.d_vb.d_x,
            crack.d_vb.d_y, crack.d_vb.length(), crack.d_vt.d_x, crack.d_vt.d_y,
            crack.d_vt.length());
  }

  d_crackOutData.d_updateCount++;
}