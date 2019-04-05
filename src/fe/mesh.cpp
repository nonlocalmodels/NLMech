// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "mesh.h"
#include "../inp/decks/meshDeck.h"
#include "../inp/policy.h"
#include "../rw/reader.h"
#include "../util/feElementDefs.h"
#include <fstream>
#include <iostream>
#include <stdint.h>

fe::Mesh::Mesh(inp::MeshDeck *deck)
    : d_numNodes(0), d_numElems(0), d_eType(0), d_eNumVertex(0), d_numDofs(0) {

  //  d_meshDeck_p = deck;
  d_dim = deck->d_dim;
  d_spatialDiscretization = deck->d_spatialDiscretization;
  d_filename = deck->d_filename;

  // perform check on input data
  if (d_spatialDiscretization != "finite_difference" and
      d_spatialDiscretization != "weak_finite_element" and
      d_spatialDiscretization != "nodal_finite_element" and
      d_spatialDiscretization != "truss_finite_element") {
    std::cerr << "Error: Spatial discretization type not known. Check input "
                 "data.\n";
    exit(1);
  }

  if (d_dim != 2) {
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
};

//
// Utility functions
//
void fe::Mesh::createData(std::string filename) {

  int file_type = -1;

  // find the extension of file and call correct reader
  if (filename.substr(filename.find_last_of(".") + 1) == "csv")
    file_type = 0;
  if (filename.substr(filename.find_last_of(".") + 1) == "vtu")
    file_type = 1;
  if (filename.substr(filename.find_last_of(".") + 1) == "msh")
    file_type = 2;

  if (d_spatialDiscretization != "finite_difference" and file_type == 0) {

    std::cerr << "Error: For discretization = " << d_spatialDiscretization
              << " .vtu or .msh mesh file is required.\n";
    exit(1);
  }

  //
  bool is_fd = false;
  if (d_spatialDiscretization == "finite_difference")
    is_fd = true;

  if (file_type == 0)
    rw::reader::readCsvFile(filename, d_dim, &d_nodes, &d_vol);
  else if (file_type == 1)
    rw::reader::readVtuFile(filename, d_dim, &d_nodes, d_eType, d_numElems,
        &d_enc, &d_nec, &d_vol, is_fd);
  else if (file_type == 2)
    rw::reader::readMshFile(filename, d_dim, &d_nodes, d_eType, d_numElems,
        &d_enc, &d_nec, &d_vol, is_fd);

  // d_numElems must be set by writer above (check!)
  d_eNumVertex = util::vtk_map_element_to_num_nodes[d_eType];
  d_numElems = d_enc.size()/d_eNumVertex;

  // get number of vertex in a given element
  d_eNumVertex = util::vtk_map_element_to_num_nodes[d_eType];

  // assign default values to fixity
  d_fix = std::vector<uint8_t>(d_nodes.size(), FREE_MASK);

  // check if we need to compute nodal volumes
  bool compute_vol = false;
  if (is_fd and d_vol.empty())
    compute_vol = true;

  // if this is weak finite element simulation then check from policy if
  // volume is to be computed
  if (d_spatialDiscretization == "weak_finite_element"
      and !inp::Policy::getInstance()->populateData("Mesh_d_vol"))
    compute_vol = false;

  if (compute_vol)
    computeVol();
}

void fe::Mesh::computeVol() {};

//
// Accessor functions
//
size_t fe::Mesh::getDimension() { return d_dim; }

size_t fe::Mesh::getNumNodes() { return d_numNodes; }

size_t fe::Mesh::getNumDofs() { return d_numDofs; }

size_t fe::Mesh::getElementType() { return d_eType; }

util::Point3 fe::Mesh::getNode(size_t i) { return d_nodes[i]; }

std::vector<util::Point3> fe::Mesh::getNodes() { return d_nodes; }

const std::vector<util::Point3> *fe::Mesh::getNodesP() { return &d_nodes; }

std::vector<size_t> fe::Mesh::getElementConnectivity(size_t i) {
  return std::vector<size_t>(d_enc.begin() + d_eNumVertex * i, d_enc.begin
  () + d_eNumVertex * i + d_eNumVertex);
}

//
// Setter functions
//
