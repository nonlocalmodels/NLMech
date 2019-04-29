// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "mshReader.h"
#include "util/feElementDefs.h"
#include <fstream>
#include <iostream>

rw::reader::MshReader::MshReader(const std::string &filename)
    : d_filename(filename){};

void rw::reader::MshReader::readMesh(size_t dim,
                                     std::vector<util::Point3> *nodes,
                                     size_t &element_type, size_t &num_elems,
                                     std::vector<size_t> *enc,
                                     std::vector<std::vector<size_t>> *nec,
                                     std::vector<double> *volumes, bool is_fd) {
  std::string line;

  // clear data
  nodes->clear();
  enc->clear();
  nec->clear();
  volumes->clear();

  // specify type of element to read
  unsigned int num_nodes_con = 0;
  if (dim != 2) {

    std::cerr << "Error: MshReader currently only supports reading of triangle "
                 "elements in dimension 2.\n";
    exit(1);
  }

  bool read_nodes = false;
  bool read_elements = false;

  // perform a mock read and find number of nodes, number of elements,
  // element type
  size_t num_nodes = 0;

  // open file
  std::ifstream mesh;
  mesh.open(d_filename.c_str());

  if (!mesh) {

    std::cerr << "Error: Can not read mesh file = " << d_filename << "\n";
    exit(1);
  }

  while (true) {

    std::getline(mesh, line);

    if (mesh) {
      // read $Nodes block
      if (line.find("$NOD") == static_cast<std::string::size_type>(0) ||
        line.find("$NOE") == static_cast<std::string::size_type>(0) ||
        line.find("$Nodes") == static_cast<std::string::size_type>(0)) {

        read_nodes = true;

        mesh >> num_nodes;

        std::cout << "number of nodes = " << num_nodes << "\n";

        // read in the nodal coordinates and form points.
        double x, y, z;
        unsigned int id;
        for (unsigned int i = 0; i < num_nodes; ++i)
          mesh >> id >> x >> y >> z;
        // read the $ENDNOD delimiter
        std::getline(mesh, line);
      } // end of reading nodes
      else if (line.find("$ELM") == static_cast<std::string::size_type>(0) ||
        line.find("$Elements") ==
          static_cast<std::string::size_type>(0)) {

        read_elements = true;

        unsigned int num_elem = 0;
        unsigned int node_id = 0;

        mesh >> num_elem;

        std::cout << "number of elems (all kinds) = " << num_elem << "\n";

        size_t elem_counter = 0;
        bool found_tri = false;
        bool found_quad = false;
        int init_tri = -1;
        int init_quad = -1;
        for (unsigned int iel = 0; iel < num_elem; ++iel) {
          unsigned int id;
          unsigned int type;
          unsigned int ntags;
          int tag;

          // read element id, type, and tags
          mesh >> id >> type >> ntags;

          // dummy read ntags
          for (unsigned int j = 0; j < ntags; j++)
            mesh >> tag;

          // read element type we desire and for other element type
          // perform dummy read
          if (type == util::msh_type_triangle) {

            if (init_tri == -1) {
              found_tri = true;
              element_type = util::vtk_type_triangle;
              num_nodes_con =
                util::msh_map_element_to_num_nodes[util::msh_type_triangle];
              init_tri = 0;
            }

            // read vertex of this element
            for (unsigned int i = 0; i < num_nodes_con; i++)
              mesh >> node_id;

            // increment the element counter
            elem_counter++;
          } else if (type == util::msh_type_quadrangle) {

            if (init_quad == -1) {
              found_quad = true;
              element_type = util::vtk_type_quad;
              num_nodes_con =
                util::msh_map_element_to_num_nodes[util::msh_type_quadrangle];
              init_quad = 0;
            }

            // read vertex of this element
            for (unsigned int i = 0; i < num_nodes_con; i++)
              mesh >> node_id;

            // increment the element counter
            elem_counter++;
          } else {
            // these are the type of elements we need to ignore
            for (unsigned int i = 0; i < util::msh_map_element_to_num_nodes[type]; i++)
              mesh >> node_id;
          }

          // check if mesh with both triangle and quadrangle elements
          if (found_quad and found_tri) {

            std::cerr << "Error: Check mesh file. It appears to have both "
                         "quadrangle elements and triangle elements. "
                         "Currently we only support one kind of elements.\n";
            exit(1);
          }
        } // element loop

        // write the number of elements
        num_elems = elem_counter;

        // read the $ENDELM delimiter
        std::getline(mesh, line);
      } // if $ELM
    } // if mesh

    // If !mesh, check to see if EOF was set.  If so, break out
    // of while loop.
    if (mesh.eof())
      break;

    if (read_nodes and read_elements)
      break;
  } // while true

  // close file
  mesh.close();

  // resize data
  nodes->resize(num_nodes);
  nec->resize(num_nodes);
  enc->resize(num_nodes_con*num_elems);

  std::cout << "Num nodes = " << num_nodes << ", Num elems = " << num_elems
  << ", element-node connectivity size = " << enc->size() << "\n";

  // read mesh file again
  mesh.open(d_filename.c_str());

  if (!mesh) {

    std::cerr << "Error: Can not read mesh file = " << d_filename << "\n";
    exit(1);
  }
  while (true) {

    std::getline(mesh, line);

    if (mesh) {
      // read $Nodes block
      if (line.find("$NOD") == static_cast<std::string::size_type>(0) ||
          line.find("$NOE") == static_cast<std::string::size_type>(0) ||
          line.find("$Nodes") == static_cast<std::string::size_type>(0)) {

        read_nodes = true;

        mesh >> num_nodes;

        // read in the nodal coordinates and form points.
        double x, y, z;
        unsigned int id;

        // add the nodal coordinates to the mesh
        for (unsigned int i = 0; i < num_nodes; ++i) {

          mesh >> id >> x >> y >> z;
          (*nodes)[id - 1] = util::Point3(x, y, z);
        }

        // read the $ENDNOD delimiter
        std::getline(mesh, line);
      } // end of reading nodes
      // Read the element block
      else if (line.find("$ELM") == static_cast<std::string::size_type>(0) ||
               line.find("$Elements") ==
                   static_cast<std::string::size_type>(0)) {
        read_elements = true;

        // For reading the number of elements and the node ids from the stream
        unsigned int num_elem = 0;
        unsigned int node_id = 0;

        // read how many elements are there
        // this includes point element, line element also
        mesh >> num_elem;

        // read the elements
        size_t elem_counter = 0;
        bool found_tri = false;
        bool found_quad = false;
        int init_tri = -1;
        int init_quad = -1;
        for (unsigned int iel = 0; iel < num_elem; ++iel) {

          unsigned int id;
          unsigned int type;
          unsigned int ntags;
          int tag;

          // read element id, type, and tags
          mesh >> id >> type >> ntags;
          for (unsigned int j = 0; j < ntags; j++)
            mesh >> tag;

          // read element type we desire and for other element type perform
          // dummy read
          if (type == util::msh_type_triangle) {

            if (init_tri == -1) {
              found_tri = true;
              element_type = util::vtk_type_triangle;
              num_nodes_con =
                  util::msh_map_element_to_num_nodes[util::msh_type_triangle];
              init_tri = 0;
            }

            std::vector<size_t> ids;

            // read vertex of this element
            for (unsigned int i = 0; i < num_nodes_con; i++) {
              mesh >> node_id;
              (*enc)[num_nodes_con * elem_counter + i] = node_id - 1;
              (*nec)[node_id - 1].push_back(elem_counter);
              ids.push_back(node_id - 1);
            }

            // check
            for (size_t k = 0; k < ids.size(); k++)
              if ((*enc)[num_nodes_con * elem_counter + k] != ids[k]) {
                std::cerr << "Error: Element-node connectivity data does not "
                             "pass check.\n";
                exit(1);
              }

            // increment the element counter
            elem_counter++;
          } else if (type == util::msh_type_quadrangle) {

            if (init_quad == -1) {
              found_quad = true;
              element_type = util::vtk_type_quad;
              num_nodes_con =
                  util::msh_map_element_to_num_nodes[util::msh_type_quadrangle];
              init_quad = 0;
            }

            std::vector<size_t> ids;

            // read vertex of this element
            for (unsigned int i = 0; i < num_nodes_con; i++) {
              mesh >> node_id;
              (*enc)[num_nodes_con * elem_counter + i] = node_id - 1;
              (*nec)[node_id - 1].push_back(elem_counter);
              ids.push_back(node_id - 1);
            }

            // check
            for (size_t k = 0; k < ids.size(); k++)
              if ((*enc)[num_nodes_con * elem_counter + k] != ids[k]) {
                std::cerr << "Error: Element-node connectivity data does not "
                             "pass check.\n";
                exit(1);
              }

            // increment the element counter
            elem_counter++;
          } else {
            // these are the type of elements we need to ignore.
            for (unsigned int i = 0; i < util::msh_map_element_to_num_nodes[type]; i++)
              mesh >> node_id;
          }

          // check if mesh with both triangle and quadrangle elements
          if (found_quad and found_tri) {

            std::cerr << "Error: Check mesh file. It appears to have both "
                         "quadrangle elements and triangle elements. "
                         "Currently we only support one kind of elements.\n";
            exit(1);
          }
        } // element loop

        // read the $ENDELM delimiter
        std::getline(mesh, line);
      } // if $ELM
    } // if mesh

    // If !mesh, check to see if EOF was set.  If so, break out
    // of while loop.
    if (mesh.eof())
      break;

    if (read_nodes and read_elements)
      break;
  } // while true

  // close file
  mesh.close();
}