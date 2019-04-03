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
                                     size_t &element_type,
                                     size_t &num_elem,
                                     std::vector<size_t> *enc,
                                     std::vector<std::vector<size_t>> *nec,
                                     std::vector<double> *volumes, bool is_fd) {

  // open file
  std::ifstream mesh;
  mesh.open(d_filename.c_str());

  if (!mesh) {

    std::cerr << "Error: Can not read mesh file = " << d_filename << "\n";
    exit(1);
  }

  std::string line;
  int format = 0;
  int size = 0;
  double version = 1.0;

  // clear data
  nodes->clear();
  enc->clear();
  nec->clear();
  volumes->clear();

  // specify type of element to read
  unsigned int el_type_read;
  unsigned int num_nodes_con;
  if (dim == 2) {

    // element type (we use vtk type definition globally). In 2-d we
    // currently only read triangle element
    element_type = util::vtk_type_triangle;

    // element type for local .msh file read. In 2-d we currently only read
    // triangle element
    el_type_read = util::msh_type_triangle;

    // number of vertex per element. Requires three node ids for connectivity.
    // num_nodes_con should be 3 if not then check feElementDefs.h
    num_nodes_con = util::msh_map_element_to_num_nodes[util::msh_type_triangle];
  } else {

    std::cerr << "Error: MshReader currently only supports reading of triangle "
                 "elements in dimension 2.\n";
    exit(1);
  }

  bool read_nodes = false;
  bool read_elements = false;

  while (true) {

    std::getline(mesh, line);

    if (mesh) {
      // // read $MeshFormat block
      //       if (line.find("$MeshFormat") ==
      //       static_cast<std::string::size_type>(0)) {

      //          mesh >> version >> format >> size;
      //          if ((version != 2.0) && (version != 2.1) && (version != 2.2))
      //          {
      //             std::cerr<<"Error: Unknown msh file version " <<
      //             version<<"\n"; exit(1);
      //          }

      //          if (format) {
      //            std::cerr<<"Error: Unknown data format for mesh in Gmsh
      //            reader.\n"; exit(1);
      //          }
      //       }
      // read $Nodes block
      if (line.find("$NOD") == static_cast<std::string::size_type>(0) ||
          line.find("$NOE") == static_cast<std::string::size_type>(0) ||
          line.find("$Nodes") == static_cast<std::string::size_type>(0)) {

        read_nodes = true;

        unsigned int num_nodes = 0;
        mesh >> num_nodes;

        // nodes = new std::vector<util::point>(num_nodes, util::point());
        nodes->resize(num_nodes);
        nec->resize(num_nodes);

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

        // As of version 2.2, the format for each element line is:
        // elm-number elm-type number-of-tags < tag > ... node-number-list
        // From the Gmsh docs:
        // * the first tag is the number of the
        //   physical entity to which the element belongs
        // * the second is the number of the elementary geometrical
        //   entity to which the element belongs
        // * the third is the number of mesh partitions to which the element
        //   belongs
        // * The rest of the tags are the partition ids (negative
        //   partition ids indicate ghost cells). A zero tag is
        //   equivalent to no tag. Gmsh and most codes using the
        //   MSH 2 format require at least the first two tags
        //   (physical and elementary tags).

        // read the elements
        size_t elem_counter = 0;
        for (unsigned int iel = 0; iel < num_elem; ++iel) {

          unsigned int id;
          unsigned int type;
          unsigned int physical = 1;
          unsigned int elementary = 1;
          unsigned int nnodes = 0;
          unsigned int ntags;

          // Note: tag has to be an int because it could be negative,
          // see above.
          int tag;

          // read element id, type, and tags
          mesh >> id >> type >> ntags;

          // dummy read ntags
          for (unsigned int j = 0; j < ntags; j++) {

            mesh >> tag;

            if (j == 0)
              physical = tag;
            else if (j == 1)
              elementary = tag;
          }

          // read element type we desire and for other element type
          // perform dummy read
          if (type == el_type_read) {

            std::vector<size_t> ids;

            // read vertex of this element
            for (unsigned int i = 0; i < num_nodes_con; i++) {

              mesh >> node_id;

              // add to the element-node connectivity
              // substract 1 to correct the numbering convention
              enc->push_back(node_id - 1);

              // also store it to perform check
              ids.push_back(node_id - 1);

              // fill the node-element connectivity table
              (*nec)[node_id - 1].push_back(elem_counter);
            }

            // check
            for (size_t k=0; k<ids.size(); k++)
              if ((*enc)[num_nodes_con * elem_counter + k] != ids[k]) {
                std::cerr << "Error: Element-node connectivity data does not "
                             "pass check.\n";
                exit(1);
              }

            // increment the element counter
            elem_counter++;
          } else {

            // these are the type of elements we need to ignore.

            size_t n = util::msh_map_element_to_num_nodes[type];

            // dummy read
            for (unsigned int i = 0; i < n; i++)
              mesh >> node_id;
          }
        } // element loop

        // write the number of elements
        num_elem = elem_counter;

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

    // If !mesh and !mesh.eof(), stream is in a bad state!
    // std::cerr<<"Error: Stream is bad! Perhaps the file does not exist?\n";
    // exit(1);

  } // while true

  // // for debug
  // int num_nodes = nodes->size();
  // int num_en_con = en_con->size();
  // int num_ne_con = ne_con->size();

  // std::cout<<"--> "<<num_nodes<<", "<<num_en_con<<", "<<num_ne_con<<"\n";
}