////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "mshReader.h"

#include <fstream>
#include <iostream>
#include <gmsh.h>

#include "util/feElementDefs.h"

rw::reader::MshReader::MshReader(const std::string &filename)
    : d_filename(filename){};

void rw::reader::MshReader::readMesh(size_t dim,
                                     std::vector<util::Point3> *nodes,
                                     size_t &element_type, size_t &num_elems,
                                     std::vector<size_t> *enc,
                                     std::vector<std::vector<size_t>> *nec,
                                     std::vector<double> *volumes, bool is_fd) {
  
  gmsh::initialize();
  gmsh::option::setNumber("General.Terminal", 1);
  gmsh::open(d_filename);

  // clear data
  nodes->clear();
  enc->clear();
  nec->clear();
  volumes->clear();

  // getting all nodes using GMSH API
  std::vector<std::size_t> nodeTags;
  std::vector<double> nodeCoords, nodeParams;
  gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, -1, -1);

  // getting all elements using GMSH API
  std::vector<int> elemTypes;
  std::vector<std::vector<std::size_t> > elemTags, elemNodeTags;
  gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, -1, -1);

  // specify type of element to read
  if (dim != 2 and dim != 3) {
    std::cerr << "Error: MshReader currently only supports reading of "
                 "triangle/quadrangle elements in dimension 2 and tetragonal "
                 "elements in 3.\n";
    exit(1);
  }


  // Convert the coordinates to the internal datastrucutre
  nodes->resize(nodeTags.size());



  for (size_t i = 0 ; i < nodeTags.size() ; i++ )
    {
    size_t index = (nodeTags[i]-1) * 3;
    (*nodes)[nodeTags[i] - 1] = util::Point3(nodeCoords[index], nodeCoords[index+1], nodeCoords[index+2]);
    }

  
 
  // Check for element types
  // 2 = 3-node triangle
  // 3 = 4-node square

  size_t type = 0;
  size_t element_id;
  
  auto f2 = std::find(elemTypes.begin(),elemTypes.end(),2);
  auto f3 = std::find(elemTypes.begin(),elemTypes.end(),3);

  if (f2 !=  elemTypes.end()){
    type = 2;
    element_id = f2 - elemTypes.begin();
    element_type = util::vtk_type_triangle; 
  }
  else if (f3 !=  elemTypes.end()){
    type = 3;
    element_id = f3 - elemTypes.begin();
    element_type = util::vtk_type_quad;
    }

  if ( f2 !=  elemTypes.end() and f3 !=  elemTypes.end() )
    std::cerr << "Error: Only mesh with one type of elements is supported!" << std::endl ;
  

  if (type = 0)
    {
      std::cerr << "Error: Only 3-node triangle or 4-node square elements are supported!" << std::endl;
      exit(1);
    }


  nec->resize(elemNodeTags[element_id].size());


    
    size_t elem_counter = 0;
    //for (size_t i = 0 ; i < elemNodeTags.size() ; i++)
    for (size_t j = 0 ; j < elemNodeTags[element_id].size() ; j++)
    {
      
      enc->push_back(elemNodeTags[element_id][j] - 1);
      (*nec)[elemNodeTags[element_id][j] - 1].push_back(elem_counter);

      if ((j + 1) % 3 == 0)
        elem_counter++;
    }
    num_elems = elemNodeTags[element_id].size(); 


    gmsh::clear();
    gmsh::finalize();
 
}

void rw::reader::MshReader::readNodes(std::vector<util::Point3> *nodes) {
  // open file
  if (!d_file) d_file = std::ifstream(d_filename);

  if (!d_file) {
    std::cerr << "Error: Can not open file = " << d_filename + ".msh.\n";
    exit(1);
  }


  gmsh::initialize();
  gmsh::option::setNumber("General.Terminal", 1);
  gmsh::open(d_filename);

  // clear data
  nodes->clear();

    // getting all nodes using GMSH API
  std::vector<std::size_t> nodeTags;
  std::vector<double> nodeCoords, nodeParams;
  gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, -1, -1);

  // Convert the coordinates to the internal datastrucutre
  nodes->resize(nodeTags.size());

  for (size_t i = 0 ; i < nodeTags.size() ; i++ )
    {
    
    std::cout << nodeTags[i] << std::endl;
    size_t index = (nodeTags[i]-1) * 3;
    (*nodes)[nodeTags[i] - 1] = util::Point3(nodeCoords[index], nodeCoords[index+1], nodeCoords[index+2]);
    }


    gmsh::clear();
    gmsh::finalize();
  
}

bool rw::reader::MshReader::readPointData(const std::string &name,
                                          std::vector<util::Point3> *data) {
  // open file
  if (!d_file) d_file = std::ifstream(d_filename);

  if (!d_file)
    if (!d_file) {
      std::cerr << "Error: Can not open file = " << d_filename + ".msh"
                << ".\n";
      exit(1);
    }

  bool found_data = false;
  std::string line;
  while (true) {
    std::getline(d_file, line);
    if (d_file) {
      // read $Nodes block
      if (line.find("$NodeData") == static_cast<std::string::size_type>(0)) {
        // get name of data
        int num_tags = 0;
        d_file >> num_tags;
        std::string tag[num_tags];
        for (size_t i = 0; i < num_tags; i++) d_file >> tag[i];

        // read dummy data
        d_file >> num_tags;
        double real_tag = 0.;
        d_file >> real_tag;

        int tag_number = 0;
        int field_type = 0;
        int num_data = 0;
        d_file >> tag_number >> field_type >> num_data;

        // check we found the data
        if (tag[0] == name) {
          // check if data is of desired field type
          if (field_type != 3) {
            std::cerr << "Error: Data " << tag[0] << " is of type "
                      << field_type << " but we expect it to be of type " << 3
                      << ".\n";
            exit(1);
          }

          found_data = true;
          data->resize(num_data);
        }

        // we read through the data irrespective of we found it or not
        for (size_t i = 0; i < num_data; i++) {
          double d[field_type];
          for (size_t j = 0; j < field_type; j++) d_file >> d[j];

          if (found_data) (*data)[i] = util::Point3(d[0], d[1], d[2]);
        }
        // read the end of data block
        std::getline(d_file, line);
      }  // end of reading nodes
    }    // if d_file

    if (found_data) break;

    // If !d_file, check to see if EOF was set.  If so, break out
    // of while loop.
    if (d_file.eof()) break;
  }  // while true

  d_file.close();
  return found_data;
}

bool rw::reader::MshReader::readPointData(const std::string &name,
                                          std::vector<double> *data) {
  // open file
  if (!d_file) d_file = std::ifstream(d_filename);

  if (!d_file)
    if (!d_file) {
      std::cerr << "Error: Can not open file = " << d_filename + ".msh"
                << ".\n";
      exit(1);
    }

  bool found_data = false;
  std::string line;
  while (true) {
    std::getline(d_file, line);
    if (d_file) {
      // read $Nodes block
      if (line.find("$NodeData") == static_cast<std::string::size_type>(0)) {
        // get name of data
        int num_tags = 0;
        d_file >> num_tags;
        std::string tag[num_tags];
        for (size_t i = 0; i < num_tags; i++) d_file >> tag[i];

        // read dummy data
        d_file >> num_tags;
        double real_tag = 0.;
        d_file >> real_tag;

        int tag_number = 0;
        int field_type = 0;
        int num_data = 0;
        d_file >> tag_number >> field_type >> num_data;

        // check we found the data
        if (tag[0] == name) {
          // check if data is of desired field type
          if (field_type != 1) {
            std::cerr << "Error: Data " << tag[0] << " is of type "
                      << field_type << " but we expect it to be of type " << 1
                      << ".\n";
            exit(1);
          }

          found_data = true;
          data->resize(num_data);
        }

        // we read through the data irrespective of we found it or not
        for (size_t i = 0; i < num_data; i++) {
          double d[field_type];
          for (size_t j = 0; j < field_type; j++) d_file >> d[j];

          if (found_data) (*data)[i] = d[0];
        }
        // read the end of data block
        std::getline(d_file, line);
      }  // end of reading nodes
    }    // if d_file

    if (found_data) break;

    // If !d_file, check to see if EOF was set.  If so, break out
    // of while loop.
    if (d_file.eof()) break;
  }  // while true

  d_file.close();
  return found_data;
}

void rw::reader::MshReader::close() {}
