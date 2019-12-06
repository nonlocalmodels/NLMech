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
#include "util/feElementDefs.h"

fe::Mesh::Mesh()
    : d_numNodes(0),
      d_numElems(0),
      d_eType(1),
      d_eNumVertex(0),
      d_numDofs(0),
      d_h(0.),
      d_dim(0) {}

fe::Mesh::Mesh(inp::MeshDeck *deck)
    : d_numNodes(0),
      d_numElems(0),
      d_eType(1),
      d_eNumVertex(0),
      d_numDofs(0),
      d_h(deck->d_h),
      d_dim(deck->d_dim),
      d_spatialDiscretization(deck->d_spatialDiscretization),
      d_filename(deck->d_filename) {
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

  if (d_dim > 2) {
    std::cerr << "Error: Check Dimension in input data. Currently we only "
                 "support dimension 2.\n";
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
  createData(d_filename);

  // check if we need to compute mesh size
  if (deck->d_computeMeshSize) computeMeshSize();
}

//
// Utility functions
//
void fe::Mesh::createData(const std::string &filename, bool ref_config) {
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
    if (!is_fd || !found_volume_data)
      rw::reader::readVtuFileCells(filename, d_dim, d_eType, d_numElems, &d_enc,
                                   &d_nec);

    // check if file has fixity data
    rw::reader::readVtuFilePointData(filename, "Fixity", &d_fix);
  }

  // compute data from mesh data
  d_numNodes = d_nodes.size();
  d_eNumVertex = util::vtk_map_element_to_num_nodes[d_eType];
  d_numDofs = d_numNodes * d_dim;

  //
  // assign default values to fixity
  //
  if (d_fix.size() != d_numNodes)
    d_fix = std::vector<uint8_t>(d_nodes.size(), uint8_t(0));

  //
  // compute nodal volume if required
  //
  bool compute_vol = false;
  if (is_fd and d_vol.empty()) compute_vol = true;

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
  d_vol.reserve(d_numNodes);
  auto f = hpx::parallel::for_loop(
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

          std::vector<fe::QuadData> qds = quads->getQuadDatas(e_nodes);

          // compute V_e and add it to volume
          for (auto qd : qds) v += qd.d_shapes[loc_i] * qd.d_w;
        }  // loop over elements

        // update
        this->d_vol[i] = v;
      });  // end of parallel for loop

  f.get();
}

void fe::Mesh::computeBBox() {
  std::vector<double> p1(3, 0.);
  std::vector<double> p2(3, 0.);
  for (auto x : d_nodes) {
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
  double guess = std::abs(d_bbox.second[0] - d_bbox.first[0]);
  if (d_dim > 1 && guess > std::abs(d_bbox.second[1] - d_bbox.first[1]))
    guess = std::abs(d_bbox.second[1] - d_bbox.first[1]);

  if (d_dim > 2 && guess > std::abs(d_bbox.second[2] - d_bbox.first[2]))
    guess = std::abs(d_bbox.second[2] - d_bbox.first[2]);

  for (size_t i = 0; i < d_nodes.size(); i++)
    for (size_t j = 0; j < d_nodes.size(); j++)
      if (i != j) {
        double val = d_nodes[i].dist(d_nodes[j]);
        if (val < guess) guess = val;
      }

  d_h = guess;
}

//
// Accessor functions
//
size_t fe::Mesh::getDimension() { return d_dim; }
size_t fe::Mesh::getDimension() const { return d_dim; }

size_t fe::Mesh::getNumNodes() { return d_numNodes; }
size_t fe::Mesh::getNumNodes() const { return d_numNodes; }

size_t fe::Mesh::getNumElements() { return d_enc.size() / d_eNumVertex; }
size_t fe::Mesh::getNumElements() const { return d_enc.size() / d_eNumVertex; }

size_t fe::Mesh::getNumDofs() { return d_numDofs; }
size_t fe::Mesh::getNumDofs() const { return d_numDofs; }

size_t fe::Mesh::getElementType() { return d_eType; }
size_t fe::Mesh::getElementType() const { return d_eType; }

double fe::Mesh::getMeshSize() { return d_h; }
double fe::Mesh::getMeshSize() const { return d_h; }

util::Point3 fe::Mesh::getNode(const size_t &i) { return d_nodes[i]; }
util::Point3 fe::Mesh::getNode(const size_t &i) const { return d_nodes[i]; }

double fe::Mesh::getNodalVolume(const size_t &i) { return d_vol[i]; }
double fe::Mesh::getNodalVolume(const size_t &i) const { return d_vol[i]; }

const std::vector<util::Point3> *fe::Mesh::getNodesP() { return &d_nodes; }
const std::vector<util::Point3> *fe::Mesh::getNodesP() const {
  return &d_nodes;
}

const std::vector<uint8_t> *fe::Mesh::getFixityP() { return &d_fix; }
const std::vector<uint8_t> *fe::Mesh::getFixityP() const { return &d_fix; }

const std::vector<double> *fe::Mesh::getNodalVolumeP() { return &d_vol; }
const std::vector<double> *fe::Mesh::getNodalVolumeP() const { return &d_vol; }

const std::vector<size_t> fe::Mesh::getElementConnectivity(const size_t &i) {
  return std::vector<size_t>(d_enc.begin() + d_eNumVertex * i,
                             d_enc.begin() + d_eNumVertex * i + d_eNumVertex);
}
const std::vector<size_t> fe::Mesh::getElementConnectivity(
    const size_t &i) const {
  return std::vector<size_t>(d_enc.begin() + d_eNumVertex * i,
                             d_enc.begin() + d_eNumVertex * i + d_eNumVertex);
}

const std::vector<util::Point3> fe::Mesh::getElementConnectivityNodes(
    const size_t &i) {
  std::vector<util::Point3> nds;
  for (size_t k = 0; k < d_eNumVertex; k++)
    nds.emplace_back(d_nodes[d_enc[d_eNumVertex * i + k]]);
  return nds;
}

const std::vector<util::Point3> fe::Mesh::getElementConnectivityNodes(
    const size_t &i) const {
  std::vector<util::Point3> nds;
  for (size_t k = 0; k < d_eNumVertex; k++)
    nds.emplace_back(d_nodes[d_enc[d_eNumVertex * i + k]]);
  return nds;
}

const std::vector<size_t> *fe::Mesh::getElementConnectivitiesP() {
  return &d_enc;
}
const std::vector<size_t> *fe::Mesh::getElementConnectivitiesP() const {
  return &d_enc;
}

const std::vector<size_t> fe::Mesh::getNodeElementConnectivity(
    const size_t &i) {
  return d_nec[i];
}
const std::vector<size_t> fe::Mesh::getNodeElementConnectivity(
    const size_t &i) const {
  return d_nec[i];
}

std::pair<std::vector<double>, std::vector<double>> fe::Mesh::getBoundingBox() {
  return d_bbox;
}
std::pair<std::vector<double>, std::vector<double>> fe::Mesh::getBoundingBox()
    const {
  return d_bbox;
}

bool fe::Mesh::isNodeFree(const size_t &i, const unsigned int &dof) {
  // below checks if d_fix has 1st bit (if dof=0), 2nd bit (if dof=1), 3rd
  // bit (if dof=2) is set to 1 or 0. If set to 1, then it means it is fixed,
  // and therefore it returns false
  return !(d_fix[i] >> dof & 1UL);
}
bool fe::Mesh::isNodeFree(const size_t &i, const unsigned int &dof) const {
  // below checks if d_fix has 1st bit (if dof=0), 2nd bit (if dof=1), 3rd
  // bit (if dof=2) is set to 1 or 0. If set to 1, then it means it is fixed,
  // and therefore it returns false
  return !(d_fix[i] >> dof & 1UL);
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
  createData(filename);

  // check if we need to compute mesh size
  if (deck->d_computeMeshSize) computeMeshSize();
}
