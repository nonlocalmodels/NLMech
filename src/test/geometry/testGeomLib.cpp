////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "testGeomLib.h"
#include "../../external/csv.h"
#include "geometry/fracture.h"
#include "inp/decks/fractureDeck.h"
#include "util/point.h"
#include <bitset>
#include <fstream>
#include <iostream>

static void readNodes(const std::string &filename,
                      std::vector<util::Point3> &nodes) {
  // csv reader
  io::CSVReader<3> in(filename);
  double x, y, z;
  while (in.read_row(x, y, z)) nodes.emplace_back(x, y, z);
}

static void printBits(const std::string &filename,
                      const std::vector<util::Point3> &nodes,
                      geometry::Fracture *fracture) {
  std::ofstream myfile(filename.c_str());
  myfile.precision(8);
  for (size_t i = 0; i < nodes.size(); i++) {
    auto list = fracture->getBonds(i);
    for (auto k : list) {
      std::bitset<8> a(k);
      myfile << a;
    }
    myfile << "\n";
  }
  myfile.close();
}

void test::testFracture() {
  auto *deck = new inp::FractureDeck();

  // create dummy list of nodes and neighborlist
  std::vector<util::Point3> nodes;
  readNodes("testFracture_nodes.csv", nodes);

  // put all the nodes in each others neighborlist
  std::vector<std::vector<size_t>> neighbor_list;
  neighbor_list.resize(nodes.size());
  for (size_t i = 0; i < nodes.size(); i++)
    for (size_t j = 0; j < nodes.size(); j++) neighbor_list[i].emplace_back(j);

  auto *fracture = new geometry::Fracture(deck, &nodes, &neighbor_list);

  //  // print bonds as bits
  //  printBits("bonds_1.csv", nodes, fracture);
  size_t error_check = 0;

  // now set the bit of neighbors of even node as fixed
  for (size_t i = 0; i < nodes.size(); i++) {
    for (size_t k = 0; k < neighbor_list[i].size(); k++) {
      if (fracture->getBondState(i, k)) error_check++;

      if (i % 2 == 0) fracture->setBondState(i, k, true);

      if (i % 2 == 0 and !fracture->getBondState(i, k)) error_check++;
    }
  }

  // now set the bit of neighbors of even node as free and odd node as fixed
  for (size_t i = 0; i < nodes.size(); i++) {
    for (size_t k = 0; k < neighbor_list[i].size(); k++) {
      fracture->setBondState(i, k, i % 2 != 0);

      if (i % 2 == 0 and fracture->getBondState(i, k)) error_check++;
      if (i % 2 == 1 and !fracture->getBondState(i, k)) error_check++;
    }
  }

  std::cout << "**********************************\n";
  std::cout << "Fracture Class Test\n";
  std::cout << "**********************************\n";
  std::cout << (error_check == 0 ? "TEST 1 : PASS. \n" : "TEST 1 : FAIL. \n");
}
