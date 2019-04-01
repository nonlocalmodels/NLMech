// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "geometry.h"
#include "../inp/decks/geometryDeck.h"
#include "../inp/policy.h"
#include <iostream>

geometry::Geometry::Geometry(inp::GeometryDeck *deck)
    : d_numNodes(0), d_numElems(0), d_eType(0), d_eNumVertex(0), d_numDofs(0) {

  d_geometryDeck_p = deck;

  // perform check on input data
  if (d_geometryDeck_p->d_spatialDiscretization != "finite_difference"
      and d_geometryDeck_p->d_spatialDiscretization != "weak_finite_element"
      and d_geometryDeck_p->d_spatialDiscretization != "nodal_finite_element"
      and d_geometryDeck_p->d_spatialDiscretization != "truss_finite_element") {
    std::cerr << "Error: Spatial discretization type not known. Check input "
                 "data.\n";
    exit(1);
  }

  if (d_geometryDeck_p->d_dim != 2) {
    std::cerr << "Error: Check Dimension in input data. Currently we only "
                 "support dimension 2.\n";
    exit(1);
  }

  if (d_geometryDeck_p->d_quadOrder == 0)
    d_geometryDeck_p->d_quadOrder = 1;
  if (d_geometryDeck_p->d_quadOrderM == 0)
    d_geometryDeck_p->d_quadOrderM = 1;

  if (d_geometryDeck_p->d_MApproxType != "exact"
      and d_geometryDeck_p->d_MApproxType != "lumped") {
    std::cerr << "Error: Check M_Matrix_Approx tag in input data. We support "
                 " exact and lumped as two options for computation of mass "
                 "matrix.\n";
    exit(1);
  }

  if (d_geometryDeck_p->d_filename.empty()) {
    std::cerr << "Error: Filename for mesh data not specified.\n";
    exit(1);
  }

  // update policy for data population
  if (d_geometryDeck_p->d_spatialDiscretization == "weak_finite_element")
    inp::Policy::getInstance()->addToTags(1, "Geometry_d_vol");

  // check memory control flag and it is 2 or higher enforce approximation of
  // mass matrix by lumping method
  if (inp::Policy::getInstance()->getMemoryControlFlag() >= 2)
    if (d_geometryDeck_p->d_MApproxType != "lumped")
      d_geometryDeck_p->d_MApproxType = "lumped";

  // read mesh data from file
  readFile(d_geometryDeck_p->d_filename);
};

//
// Utility functions
//
void geometry::Geometry::readFile(std::string filename) {

  int file_type = -1;

  // find the extension of file and call correct reader
  if (filename.substr(filename.find_last_of(".") + 1) == "csv")
    file_type = 0;
  if (filename.substr(filename.find_last_of(".") + 1) == "vtu")
    file_type = 1;
  if (filename.substr(filename.find_last_of(".") + 1) == "msh")
    file_type = 2;

  if (d_geometryDeck_p->d_spatialDiscretization != "finite_difference" and
      file_type == 0) {

    std::cerr << "Error: For discretization = "
              << d_geometryDeck_p->d_spatialDiscretization
              << " .vtu or .msh mesh file is required.\n";
    exit(1);
  }

//  if (file_type == 0)
//    readCsvFile(filename);
//  else if (file_type == 1)
//    readVtuFile(filename);
//  else if (file_type == 2)
//    readMshFile(filename);
}


//
// Accessor functions
//
size_t geometry::Geometry::getDimension() { return d_geometryDeck_p->d_dim; }

size_t geometry::Geometry::getNumNodes() { return d_numNodes; }

size_t geometry::Geometry::getNumDofs() { return d_numDofs; }

util::Point3 geometry::Geometry::getNode(size_t i) { return d_nodes[i]; }

std::vector<util::Point3> geometry::Geometry::getNodes() { return d_nodes; }

const std::vector<util::Point3> *geometry::Geometry::getNodesP() {
  return &d_nodes;
}


//
// Setter functions
//

