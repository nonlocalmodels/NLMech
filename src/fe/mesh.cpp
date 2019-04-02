// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "mesh.h"
#include "../inp/decks/meshDeck.h"
#include "../inp/policy.h"
#include "../rw/reader.h"
#include <fstream>
#include <iostream>

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
    inp::Policy::getInstance()->addToTags(1, "Geometry_d_vol");

  // read mesh data from file
  readFile(d_filename);
};

//
// Utility functions
//
void fe::Mesh::readFile(std::string filename) {

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

  bool is_fd = false;
  if (d_spatialDiscretization == "finite_difference")
    is_fd = true;

  if (file_type == 0)
    rw::reader::readCsvFile(filename, d_dim, &d_nodes, &d_vol);
  else if (file_type == 1)
    rw::reader::readVtuFile(filename, d_dim, &d_nodes, d_eType, &d_enc, &d_nec,
                            &d_vol, is_fd);
  else if (file_type == 2)
    rw::reader::readMshFile(filename, d_dim, &d_nodes, d_eType, &d_enc, &d_nec,
                            &d_vol, is_fd);

  //  std::cout<<"Num nodes = "<<d_nodes.size()<<" "<<d_enc.size()<<" "
  //                                                                ""<<d_nec
  //                                                                .size()<<"\n";
  //
  // std::cout<<d_nodes[0].d_x<<" "<<d_nodes[0].d_y<<" "<<d_nodes[0].d_z<<"\n";
  //  std::cout<<d_nodes[9].d_x<<" "<<d_nodes[9].d_y<<" "<<d_nodes[9].d_z<<"\n";
  //
  //    std::ofstream myfile("nodes.csv");
  //  myfile.precision(6);
  //  for (size_t i=0; i<10; i++)
  //    myfile << i <<","<< d_nodes[i].d_x << "," << d_nodes[i].d_y << "," <<
  //    d_nodes[i].d_z <<"\n";
  //  myfile.close();
}

//
// Accessor functions
//
size_t fe::Mesh::getDimension() { return d_dim; }

size_t fe::Mesh::getNumNodes() { return d_numNodes; }

size_t fe::Mesh::getNumDofs() { return d_numDofs; }

util::Point3 fe::Mesh::getNode(size_t i) { return d_nodes[i]; }

std::vector<util::Point3> fe::Mesh::getNodes() { return d_nodes; }

const std::vector<util::Point3> *fe::Mesh::getNodesP() { return &d_nodes; }

//
// Setter functions
//
