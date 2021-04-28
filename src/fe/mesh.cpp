////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "mesh.h"

#include <stdint.h>
#include <util/compare.h>

#include <cstdint>
#include <fstream>
#include <hpx/include/parallel_algorithm.hpp>
#include <iostream>

#include "inp/decks/meshDeck.h"
#include "inp/policy.h"
#include "lineElem.h"
#include "quadElem.h"
#include "rw/reader.h"
#include "triElem.h"
#include "tetElem.h"
#include "util/compare.h"
#include "util/feElementDefs.h"
#include "util/utilGeom.h"
#include "util/utilIO.h"

fe::Mesh::Mesh(size_t dim)
    : d_numNodes(0),
      d_numElems(0),
      d_eType(1),
      d_eNumVertex(0),
      d_numDofs(0),
      d_h(0.),
      d_dim(dim),
      d_keepElementConn(false) {}

fe::Mesh::Mesh(inp::MeshDeck *deck)
    : d_numNodes(0),
      d_numElems(0),
      d_eType(1),
      d_eNumVertex(0),
      d_numDofs(0),
      d_h(deck->d_h),
      d_dim(deck->d_dim),
      d_spatialDiscretization(deck->d_spatialDiscretization),
      d_filename(deck->d_filename),
      d_keepElementConn(deck->d_keepElementConn) {
  // perform check on input data
  if (d_spatialDiscretization != "finite_difference" and
      d_spatialDiscretization != "weak_finite_element" and
      d_spatialDiscretization != "weak_finite_element" and
      d_spatialDiscretization != "nodal_finite_element" and
      d_spatialDiscretization != "truss_finite_element") {
    std::cerr << "Error: Spatial discretization type not known. Check input "
                 "data.\n";
    exit(1);
  }

  if (d_dim < 0 or d_dim > 3) {
    std::cerr << "Error: Check Dimension in input data.\n";
    exit(1);
  }

  if (d_filename.empty()) {
    std::cerr << "Error: Filename for mesh data not specified.\n";
    exit(1);
  }

  // update policy for data population
  if (d_spatialDiscretization == "weak_finite_element")
    inp::Policy::getInstance()->addToTags(1, "Mesh_d_vol");

  // read mesh data from file
  createData(d_filename, false, deck->d_isCentroidBasedDiscretization,
             deck->d_loadPUMData);

  // check if we need to compute mesh size
  if (deck->d_computeMeshSize) computeMeshSize();

  // check nodal volume
  size_t counter = 0;
  for (const auto &v : d_vol) {
    if (v < 0.01 * std::pow(d_h, d_dim)) {
      std::cerr << "Error: Check nodal volume " << v << " is less than "
                << 0.01 * std::pow(d_h, d_dim) << ", Node = " << counter
                << " at position = " << d_nodes[counter].printStr() << "\n";

      exit(1);
    }

    counter++;
  }
}

//
// Utility functions
//
void fe::Mesh::createData(const std::string &filename, bool ref_config,
                          bool is_centroid_based, bool has_coupling_data) {
  int file_type = -1;

  // find the extension of file and call correct reader
  if (filename.substr(filename.find_last_of(".") + 1) == "csv") file_type = 0;
  if (filename.substr(filename.find_last_of(".") + 1) == "msh") file_type = 1;
  if (filename.substr(filename.find_last_of(".") + 1) == "vtu") file_type = 2;

  if (d_spatialDiscretization != "finite_difference" and file_type == 0) {
    std::cerr << "Error: For discretization = " << d_spatialDiscretization
              << " .vtu or .msh mesh file is required.\n";
    exit(1);
  }

  //
  bool is_fd = false;
  if (d_spatialDiscretization == "finite_difference") is_fd = true;

  //
  // read node and elements
  //
  std::cout << "Mesh: Reading mesh.\n";
  if (file_type == 0)
    rw::reader::readCsvFile(filename, d_dim, &d_nodes, &d_vol);
  else if (file_type == 1)
    rw::reader::readMshFile(filename, d_dim, &d_nodes, d_eType, d_numElems,
                            &d_enc, &d_nec, &d_vol, false);
  else if (file_type == 2) {
    //
    // old reading of mesh
    //
    // rw::reader::readVtuFile(filename, d_dim, &d_nodes, d_eType, d_numElems,
    //                         &d_enc, &d_nec, &d_vol, false);

    //
    // new reading of mesh
    // We read the data from file one by one depending on what data we need
    //

    // read node
    rw::reader::readVtuFileNodes(filename, d_dim, &d_nodes, ref_config);

    // read volume if required
    bool found_volume_data = false;
    if (is_fd) {
      found_volume_data =
          rw::reader::readVtuFilePointData(filename, "Node_Volume", &d_vol);

      // try another tag for nodal volume
      if (!found_volume_data)
        found_volume_data =
            rw::reader::readVtuFilePointData(filename, "Volume", &d_vol);
    }

    // read element data (only if this is fe simulation or if we need
    // element-node connectivity data for nodal volume calculation)
    if (!is_fd || !found_volume_data || d_keepElementConn)
      rw::reader::readVtuFileCells(filename, d_dim, d_eType, d_numElems, &d_enc,
                                   &d_nec);

    // check if file has fixity data
    rw::reader::readVtuFilePointData(filename, "Fixity", &d_fix);

    // Read the data for coupling
    if (has_coupling_data) {
      rw::reader::readVtuFilePointData(filename, "Boundary-Layer",
                                       &d_prescribed_nodes);

      rw::reader::readVtuFilePointData(filename, "PUM-Displacement",
                                       &d_prescribed_values);
    }
  }

  // compute data from mesh data
  d_numNodes = d_nodes.size();
  d_eNumVertex = util::vtk_map_element_to_num_nodes[d_eType];
  d_numDofs = d_numNodes * d_dim;

  //
  // compute nodal volume if required
  //
  bool compute_vol = false;
  if (is_fd and d_vol.empty()) compute_vol = true;

  // see if we want to create nodes at the centroid of each element
  if (d_numElems > 0 and is_centroid_based and is_fd) {
    nodesAtCentroid();

    // above also computes the vol so we do need to compute again
    compute_vol = false;

    // element data is now useless
    d_enc.clear();
    d_numElems = 0;
    // check if we should retain element data
    if (!d_keepElementConn) {
      d_enc.clear();
      d_numElems = 0;
    }

    // update
    d_numNodes = d_nodes.size();
    d_eNumVertex = 1;
    d_numDofs = d_numNodes * d_dim;
  }

  //
  // assign default values to fixity
  //
  d_fix = std::vector<uint8_t>(d_nodes.size(), uint8_t(0));

  // if this is weak finite element simulation then check from policy if
  // volume is to be computed
  if (d_spatialDiscretization == "weak_finite_element" and
      !inp::Policy::getInstance()->populateData("Mesh_d_vol"))
    compute_vol = false;

  if (compute_vol) {
    std::cout << "Mesh: Computing nodal volume.\n";
    computeVol();
  }

  //
  // compute bounding box
  //
  computeBBox();
}

void fe::Mesh::computeVol() {
  // initialize quadrature data
  fe::BaseElem *quads;
  if (d_eType == util::vtk_type_triangle)
    quads = new fe::TriElem(2);
  else if (d_eType == util::vtk_type_quad)
    quads = new fe::QuadElem(2);
  else if (d_eType == util::vtk_type_tetra)
    quads = new fe::TetElem(2);

  // check if we have valid element-node connectivity data for nodal volume
  // calculations
  if (d_nec.size() != d_numNodes || d_enc.empty()) {
    std::cerr << "Error: Can not compute nodal volume for given finite "
                 "element mesh as the element-node connectivity data is "
                 "invalid."
              << std::endl;
  }

  if (false) {
    print(0, 0);
    std::cout << "\n-------- Node data ----------\n";
    std::cout << util::io::printStr(d_nodes, 0) << "\n";
    std::cout << "\n-------- Element data ----------\n";
    std::cout << util::io::printStr(d_enc, 0) << "\n";
  }

  // check if we have volid element-node connectivity data for nodal volume
  // calculations
  if (d_nec.size() != d_numNodes || d_enc.empty()) {
    std::cerr << "Error: Can not compute nodal volume for given finite "
                 "element mesh as the element-node connectivity data is "
                 "invalid."
              << std::endl;
  }

  //
  // compute nodal volume
  //
  d_vol.resize(d_numNodes);
  auto f = hpx::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      this->d_numNodes, [this, quads](boost::uint64_t i) {
        double v = 0.0;

        for (auto e : this->d_nec[i]) {
          std::vector<size_t> e_ns = this->getElementConnectivity(e);

          // locate global node i in local list of element el
          int loc_i = -1;
          for (size_t l = 0; l < e_ns.size(); l++)
            if (e_ns[l] == i) loc_i = l;

          if (loc_i == -1) {
            std::cerr << "Error: Check node element connectivity.\n";
            exit(1);
          }

          // get quad data
          std::vector<util::Point3> e_nodes;
          for (auto k : e_ns) e_nodes.emplace_back(this->d_nodes[k]);

          // get volume of element
          double vol = quads->elemSize(e_nodes);
          double factor = 1.;
          if (vol < 0.) factor = -1.;

          std::vector<fe::QuadData> qds = quads->getQuadDatas(e_nodes);

          // compute V_e and add it to volume
          for (auto qd : qds) v += qd.d_shapes[loc_i] * factor * qd.d_w;
        }  // loop over elements

        // update
        this->d_vol[i] = v;
      });  // end of parallel for loop

  f.get();
}

void fe::Mesh::computeBBox() {
  std::vector<double> p1(3, 0.);
  std::vector<double> p2(3, 0.);
  for (const auto &x : d_nodes) {
    if (util::compare::definitelyLessThan(x.d_x, p1[0])) p1[0] = x.d_x;
    if (util::compare::definitelyLessThan(x.d_y, p1[1])) p1[1] = x.d_y;
    if (util::compare::definitelyLessThan(x.d_z, p1[2])) p1[2] = x.d_z;
    if (util::compare::definitelyLessThan(p2[0], x.d_x)) p2[0] = x.d_x;
    if (util::compare::definitelyLessThan(p2[1], x.d_y)) p2[1] = x.d_y;
    if (util::compare::definitelyLessThan(p2[2], x.d_z)) p2[2] = x.d_z;
  }

  d_bbox = std::make_pair(p1, p2);
}

void fe::Mesh::computeMeshSize() {
  double guess = 0.;
  if (d_nodes.size() < 2) {
    d_h = 0.;
    return;
  }

  guess = (d_nodes[0] - d_nodes[1]).length();
  for (size_t i = 0; i < d_nodes.size(); i++)
    for (size_t j = 0; j < d_nodes.size(); j++)
      if (i != j) {
        double val = d_nodes[i].dist(d_nodes[j]);

        if (util::compare::definitelyLessThan(val, 1.0E-12)) {
          std::cout << "Check nodes are too close = "
                    << util::io::printStr<util::Point3>(
                           {d_nodes[i], d_nodes[j]})
                    << "\n";
          std::cout << "Distance = " << val << ", guess = " << guess << "\n";
        }
        if (util::compare::definitelyLessThan(val, guess)) guess = val;
      }

  d_h = guess;
}

//
// Setter functions
//

void fe::Mesh::setMeshData(const size_t &dim, std::vector<util::Point3> &nodes,
                           std::vector<double> &volumes) {
  d_spatialDiscretization = "finite_difference";

  // set dimension
  d_dim = dim;

  // set node coordinates
  d_nodes = nodes;
  d_numNodes = d_nodes.size();
  d_numDofs = d_numNodes * d_dim;

  // set nodal volume
  d_vol = volumes;

  // initialize fixity
  d_fix = std::vector<uint8_t>(d_nodes.size(), uint8_t(0));

  // compute bounding box and mesh size
  computeBBox();
  computeMeshSize();
}

void fe::Mesh::setNodes(std::vector<util::Point3> &nodes) {
  // set node data
  d_nodes = nodes;
  d_numNodes = nodes.size();
}

void fe::Mesh::setNodalVolumes(std::vector<double> &volumes) {
  // set data
  d_vol = volumes;
}

void fe::Mesh::setFixity(std::vector<uint8_t> &fixity) {
  // set data
  d_fix = fixity;
}

void fe::Mesh::setFixity(const size_t &i, const unsigned int &dof,
                         const bool &flag) {
  // to set i^th bit as true of integer a,
  // a |= 1UL << (i % 8)

  // to set i^th bit as false of integer a,
  // a &= ~(1UL << (i % 8))

  flag ? (d_fix[i] |= 1UL << dof) : (d_fix[i] &= ~(1UL << dof));
}

void fe::Mesh::setMeshSize(const double &h) { d_h = h; }

void fe::Mesh::clearElementData() {
  if (!d_enc.empty()) d_enc.shrink_to_fit();
  d_numElems = 0;
  if (!d_nec.empty()) d_nec.shrink_to_fit();
}

void fe::Mesh::readFromFile(inp::MeshDeck *deck, const std::string &filename) {
  d_h = deck->d_h;
  d_dim = deck->d_dim;
  d_spatialDiscretization = deck->d_spatialDiscretization;
  d_filename = filename;

  // read file
  createData(filename, false, deck->d_isCentroidBasedDiscretization,
             deck->d_loadPUMData);

  // check if we need to compute mesh size
  if (deck->d_computeMeshSize) computeMeshSize();
}

std::string fe::Mesh::printStr(int nt, int lvl) const {
  auto tabS = util::io::getTabS(nt);
  std::ostringstream oss;
  oss << tabS << "------- Mesh --------" << std::endl << std::endl;
  oss << tabS << "Dimension = " << d_dim << std::endl;
  oss << tabS << "Spatial discretization type = " << d_spatialDiscretization
      << std::endl;
  oss << tabS << "Mesh size = " << d_h << std::endl;
  oss << tabS << "Num nodes = " << d_numNodes << std::endl;
  oss << tabS << "Num elements = " << d_numElems << std::endl;
  oss << tabS << "Element type = " << d_eType << std::endl;
  oss << tabS << "Num nodes per element = " << d_eNumVertex << std::endl;
  oss << tabS << "Num nodal vol = " << d_vol.size() << std::endl;
  oss << tabS << "Bounding box: " << std::endl;
  oss << util::io::printBoxStr(d_bbox, nt + 1);
  oss << tabS << std::endl;

  return oss.str();
}

void fe::Mesh::nodesAtCentroid() {
  // store node data temporarily
  auto fe_nodes = std::vector<util::Point3>(d_numElems, util::Point3());

  // resize vol data
  d_vol = std::vector<double>(d_numElems, 0.);

  // loop over all elements
  for (size_t i = 0; i < d_numElems; i++) {
    // get nodes of this element
    auto el_nodes = getElementConnectivityNodes(i);

    // compute centroid of the element and volume of the element
    auto center_and_vol = util::geometry::getCenterAndVol(el_nodes, d_eType);

    // add to node data
    fe_nodes[i] = center_and_vol.first;
    d_vol[i] = center_and_vol.second;
  }

  // delete old data and replace with new node data
  d_nodes = fe_nodes;
}
