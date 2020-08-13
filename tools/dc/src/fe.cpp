////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <fe/triElem.h>
#include <util/feElementDefs.h>

#include "dcInclude.h"
#include "fe/mesh.h"
#include "inp/decks/meshDeck.h"
#include "rw/reader.h"
#include "rw/writer.h"
#include "util/compare.h"
#include "util/matrix.h"
#include "util/point.h"
#include "util/transfomation.h"
#include "util/utilGeom.h"

static int init = -1;

/*! @brief Local namespace */
namespace fe {

/*! @brief Data structure to hold simulation data in one place */
struct SimData {
  /*!
   * @brief Constructor
   */
  SimData()
      : d_isFd(true),
        d_mesh_p(nullptr),
        d_h(0.),
        d_Nt(1),
        d_dt(0.),
        d_dtOut(1){};

  /*! @brief Is this finite difference simulation */
  bool d_isFd;

  /*! @brief Mesh filename */
  std::string d_filename;

  /*! @brief Pointer to mesh data */
  fe::Mesh *d_mesh_p;

  /*! @brief Current position of all nodes */
  std::vector<util::Point3> d_y;

  /*! @brief Current displacement of all nodes */
  std::vector<util::Point3> d_u;

  /*! @brief Mesh size */
  double d_h;

  /*! @brief Number of time steps */
  size_t d_Nt;

  /*! @brief Size of time steps */
  double d_dt;

  /*! @brief Output frequency */
  size_t d_dtOut;
};

struct DataCompare {
  /*!
   * @brief Constructor
   */
  DataCompare()
      : d_isFd(true), d_numQuads(1), d_diffAtCurrent(true), d_read12(true){};

  /*! @brief Is this finite difference simulation */
  bool d_isFd;

  /*! @brief Number of quadrature points */
  size_t d_numQuads;

  /*! @brief Compute L2 difference at current configuration */
  bool d_diffAtCurrent;

  /*! @brief Path for simulation 1 data */
  std::string d_path1;

  /*! @brief Path for simulation 2 data */
  std::string d_path2;

  /*! @brief Path for output of file */
  std::string d_pathOut;

  /*! @brief Tag for output of file */
  std::string d_filenameOut;

  /*! @brief Is this run for reading 1 and 2 tags in input file */
  bool d_read12;
};

/*! @brief Read input files
 *
 * @param sim_1 Simulation 1 data
 * @param sim_2 Simulation 2 data
 * @param dc Data comparison options
 * @param config YAML input file
 */
void readInputFile(SimData *sim_1, SimData *sim_2, DataCompare *dc,
                   YAML::Node config) {
  //
  // if read_12 =  true
  //
  // read Data_1 and Data_2 and provide output filename for
  // error between Data_1 and Data_2
  //
  // if read_12 =  false
  //
  // read Data_2 and Data_3 and provide output filename for
  // error between Data_2 and Data_3
  //
  dc->d_pathOut = config["Output"]["Path"].as<std::string>();
  dc->d_filenameOut = dc->d_pathOut + "/";
  if (config["Output"]["File"])
    dc->d_filenameOut += config["Output"]["File"].as<std::string>();
  else
    dc->d_filenameOut += "dc";

  // below tag is counter intuitive but should not be changed
  // as the python scripts understand error files based on below tags
  bool triple_data = false;
  if (config["Triple_Data"]) triple_data = config["Triple_Data"].as<bool>();

  // append tags only if triple_data is true
  if (triple_data) {
    if (dc->d_read12)
      dc->d_filenameOut += "_23.out";
    else
      dc->d_filenameOut += "_12.out";
  } else
    dc->d_filenameOut += ".out";

  // string name for data 1 and 2 (in case of read_12 is false data 1 is Data_2
  // and data 2 is Data_3)
  std::string data_1;
  std::string data_2;
  if (dc->d_read12) {
    data_1 = "Data_1";
    data_2 = "Data_2";
  } else {
    data_1 = "Data_2";
    data_2 = "Data_3";
  }

  // check if data is provided
  if (!config[data_1]) {
    std::cerr << "Input file does not contain information about data = "
              << data_1 << "\n";
    exit(1);
  }

  if (!config[data_2]) {
    std::cerr << "Input file does not contain information about data = "
              << data_2 << "\n";
    exit(1);
  }

  // is this fd simulation comparison
  dc->d_isFd = false;

  //
  if (config["Diff_Find_Current"])
    dc->d_diffAtCurrent = config["Diff_Find_Current"].as<bool>();
  else
    dc->d_diffAtCurrent = true;

  if (config["Num_Quads"])
    dc->d_numQuads = config["Num_Quads"].as<size_t>();
  else
    dc->d_numQuads = 4;

  // check if time is provided globally
  bool time_global = false;
  if (config["Final_Time"]) {
    time_global = true;

    auto time = config["Final_Time"].as<double>();
    sim_1->d_Nt = config["Time_Steps"].as<size_t>();
    sim_2->d_Nt = sim_1->d_Nt;

    sim_1->d_dt = time / double(sim_1->d_Nt);
    sim_2->d_dt = sim_1->d_dt;

    sim_1->d_dtOut = config["Output_Interval"].as<size_t>();
    sim_2->d_dtOut = sim_1->d_dtOut;
  }

  // check if horizon is provided
  bool horizon_global = false;
  double horizon = 0.0;
  if (config["Horizon"]) {
    horizon_global = true;
    horizon = config["Horizon"].as<double>();

    // compute mesh size of fine mesh by reading the ratio
    auto r = config[data_1]["Horizon_Factor_m_value"].as<size_t>();
    sim_1->d_h = horizon / double(r);

    // coarse mesh
    r = config[data_2]["Horizon_Factor_m_value"].as<size_t>();
    sim_2->d_h = horizon / double(r);
  }

  // read data 1
  if (!horizon_global) {
    sim_1->d_h = config[data_1]["Mesh_Size"].as<double>();
  }

  if (!time_global) {
    sim_1->d_Nt = config[data_1]["Time_Steps"].as<size_t>();
    auto d = config[data_1]["Final_Time"].as<double>();
    sim_1->d_dt = d / double(sim_1->d_Nt);
    sim_1->d_dtOut = config[data_1]["Output_Interval"].as<size_t>();
  }

  dc->d_path1 = config[data_1]["Path"].as<std::string>();

  // check if mesh file name is provided
  if (config[data_1]["Mesh_Filename"])
    sim_1->d_filename = config[data_1]["Mesh_Filename"].as<std::string>();

  // read data 2
  if (!horizon_global) {
    sim_2->d_h = config[data_2]["Mesh_Size"].as<double>();
  }

  if (!time_global) {
    sim_2->d_Nt = config[data_2]["Time_Steps"].as<size_t>();
    auto d = config[data_2]["Final_Time"].as<double>();
    sim_2->d_dt = d / double(sim_2->d_Nt);
    sim_2->d_dtOut = config[data_2]["Output_Interval"].as<size_t>();
  }

  dc->d_path2 = config[data_2]["Path"].as<std::string>();

  // check if mesh file name is provided
  if (config[data_1]["Mesh_Filename"])
    sim_1->d_filename = config[data_1]["Mesh_Filename"].as<std::string>();
}

/*! @brief Computes displacement of coarse mesh at quadrature point of fine mesh
 *
 * @param x1 Quadrature point in element of mesh 1
 * @param e1 Element id in mesh 1
 * @param q Quadrature id
 * @param sim1 Simulation 1 data
 * @param sim2 Simulation 2 data
 * @param dc Collection of user-specified input data
 * @param read_counter Read counter used in bypassing the calculation of
 * element id in mesh 2 (if specified in the input file)
 * @param e2 Element id containing x1 in mesh 2
 * @return u Displacement at the quad point
 */
util::Point3 getDisplacementAtQuadPoint(
    const util::Point3 &x1, const size_t &e1, const size_t &q,
    const SimData *sim1, const SimData *sim2, const DataCompare *dc,
    const size_t &read_counter, size_t &e2) {
  // get the pointer to mesh 2
  const auto mesh2 = sim2->d_mesh_p;

  //
  // find the node which is closest to x1. Use reference
  // configuration of node.
  //
  // Note that x1 here is quadrature point which is again in
  // reference configuration.
  //
  // Therefore, we can save the node id and element id for
  // reuse for each quadrature point of each element.
  //

  // see if we need to compute the element id which contains
  // point x1 or reuse previous computed id
  bool compute_el_id = !(!dc->d_diffAtCurrent and read_counter > 2);

  size_t el_id = 0;
  if (compute_el_id) {
    static int debug_counter = 0;
    static bool debug_el_ids = false;

    std::string filename = dc->d_pathOut + "/debug_find_elem.txt";
    static FILE *fout = NULL;

    if (debug_el_ids) fout = fopen(filename.c_str(), "w");

    //
    double h = sim2->d_h;

    // loop over node and find the node which is closest to the point
    int i_found = -1;
    double check = std::sqrt(2.) * 0.5 * h + 0.01 * h;
    double dist = 10.0 * h;
    for (size_t i = 0; i < mesh2->getNumNodes(); i++) {
      auto diff = x1 - mesh2->getNode(i);

      if (util::compare::definitelyLessThan(diff.length(), dist)) {
        dist = diff.length();
        i_found = i;
      }
    }
    if (i_found == -1) {
      std::cerr << "Node closer to x = " << x1.printStr() << " not found."
                << " h = " << h << ", dist = " << dist << ", check = " << check
                << std::endl;
      exit(1);
    }

    // if (dist >  check) {
    //    std::cerr<<
    //       "Node closer to x = ["<<x.x<<","<<x.y<<
    //       "] not found. h = "<<
    //       h<<", dist = "<<
    //       dist<<" check = "<<
    //       check<<".\n";
    //    exit(1);
    // }

    // now that we have found the node closest to xi
    // we loop over the elements sharing the node
    //
    // We look at the intersection of line xi \to xnode and line between two
    // other nodes in triangle
    //
    // If distance of intersection point and xi is smaller than the
    // distance of intersection point and xnode than xi will be
    // inside.
    //
    // At the same time we also test if line between xnode and
    // point xi is parallel to any of the edge of
    // element

    auto xi = mesh2->getNode(i_found);
    auto yi = sim2->d_y[i_found];
    auto i_elems = mesh2->getElementConnectivity(i_found);
    int e_found = -1;
    for (unsigned long e_id : i_elems) {
      auto e_nodes = mesh2->getElementConnectivity(e_id);

      // find the local id of node i_found in element el
      int i_loc = -1;
      for (size_t j = 0; j < e_nodes.size(); j++)
        if (e_nodes[j] == i_found) i_loc = j;

      if (i_loc == -1) {
        std::cerr << "Error: Check element node connectivity.\n";
        exit(1);
      }

      // find the nodes order acw starting from i_found
      auto node_list =
          util::transformation::cyclicOrderACW(i_loc, e_nodes.size());

      // get the angle between two lines
      //
      //
      //
      //                x'
      //                +
      //             \  |  /
      //              \ | /
      //               \|/
      //                o i = i_found
      //               /|\
      //              / | \
      //          L1 /  |  \  L2
      //            /   |L0 \
      //           /    x    \
      //          /           \
      //         o             o
      //        i+1            i-1
      //
      // Line 1: (x_i, x_i+1)
      // Line 2: (x_i, x_i-1)
      //
      // Line 0: (x_i, x)
      //
      // Angle theta_12: angle between Line 1 and Line 2
      // Angle theta_01: angle between Line 1 and Line 0
      // Angle theta_02: angle between Line 2 and Line 0
      //
      // If theta_12 >= theta_01 and theta_12 >= theta_02
      // ---> Then Line 0 lies between Line 1 and Line 2
      // ---> There are two possibilities
      //      ---> Case 1: x is indeed in the interior of
      //           element formed by i,i+1,i+2,...,i-1
      //      ---> Case 2: x is in the location x'
      //           Note that x' also satisifies the similar angular
      //           condition
      //
      // ---> We eliminate Case 2, by looking at the area
      //      of triangle T0 : (x, x_i+1, x_i-1)
      //         triangle T1 : (xi, x_i+1, x_i-1)
      //
      //      If x inside the triangle then Area(T0) < Area(T1)
      //

      // get the ids of nodes
      auto i = node_list[0];
      auto ip1 = node_list[1];
      auto im1 = node_list[node_list.size() - 1];

      xi = mesh2->getNode(e_nodes[i]);
      auto xip1 = mesh2->getNode(e_nodes[ip1]);
      auto xim1 = mesh2->getNode(e_nodes[im1]);

      double theta_12 = (xip1 - xi).angle(xim1 - xi);
      double theta_01 = (xip1 - xi).angle(x1 - xi);
      double theta_02 = (x1 - xi).angle(xim1 - xi);

      std::vector<util::Point3> nodes = {x1, xip1, xim1};

      if (theta_01 <= 1.0E-4 or theta_02 <= 1.0E-4 or
          theta_01 >= M_PI - 1.0E-4 or theta_02 >= M_PI - 1.0E-4) {
        // Line 0 is either parallel to Line 1 or Line 2

       
        // we now eliminate the other possibility
        double area_T0 = std::abs(util::geometry::getTriangleArea(nodes));
        double area_T1 = std::abs(util::geometry::getTriangleArea(nodes));

        if (area_T0 <= area_T1) {
          // we have found the element

          // write the global id
          e_found = int(e_id);
        }
      }

      // proceed ahead only if condition is met
      if (util::compare::definitelyLessThan(theta_01 - theta_12, 0.0) and
          util::compare::definitelyLessThan(theta_02 - theta_12, 0.0)) {
        // we now eliminate the other possibility
        double area_T0 = std::abs(util::geometry::getTriangleArea(nodes));
        double area_T1 = std::abs(util::geometry::getTriangleArea(nodes));

        if (util::compare::definitelyLessThan(area_T0 - area_T1, 0.0) or
            util::compare::essentiallyEqual(area_T0 - area_T1, 0.0)) {
          // we have found the element

          // write the global id
          e_found = int(e_id);
        }
      }

      if (e_found != -1) break;
    }

    if (e_found == -1) {
      std::cerr << "Error: Can not locate element in the mesh for point.\n";
      std::cerr << " x = " << x1.printStr() << ", x_node = " << xi.printStr()
                << ", dist = " << dist << ", h = " << h << "\n";
      std::cerr << "Element id in mesh 1 = " << e1 << "\n";
      for (size_t j = 0; j < i_elems.size(); j++) {
        std::cerr << "Element " << j + 1 << ": global id = " << i_elems[j]
                  << "; nodes = ";

        std::vector<size_t> jel_nodes =
            mesh2->getElementConnectivity(i_elems[j]);

        for (unsigned long jel_node : jel_nodes) {
          std::cerr << "(id = " << jel_node
                    << ", x = " << mesh2->getNode(jel_node).printStr() << ")";
        }
        std::cerr << std::endl;
      }
      exit(1);
    } else {
      // write debug output
      if (debug_counter % 10 == 0 and debug_counter < 500 and debug_el_ids) {
        std::vector<size_t> jel_nodes =
            mesh2->getElementConnectivity(i_elems[e_found]);

        fprintf(fout, "----------------------\n");
        fprintf(fout,
                "x = [%s]; xnode = [%s]; dist = %6.8e; h = "
                "%6.8e;\n",
                x1.printStr().c_str(), xi.printStr().c_str(), dist, h);
        fprintf(fout,
                "element = %lu; x1 = [%s]; x2 = [%s]; x3 = "
                "[%s];\n",
                i_elems[e_found],
                mesh2->getNode(jel_nodes[0]).printStr().c_str(),
                mesh2->getNode(jel_nodes[1]).printStr().c_str(),
                mesh2->getNode(jel_nodes[2]).printStr().c_str());
      }
    }

    debug_counter++;
    if (debug_counter == 500 and debug_el_ids) fclose(fout);

    //
    el_id = size_t(e_found);

    // we have found the element id. See if we need to
    // add it to list for future reuse
    if (read_counter == 1) e2 = el_id;

    // check if this is not the first call to this function
    if (read_counter > 1 and e2 != el_id) {
      std::cout << "Error: stored element in mesh 2 = " << e2
                << " and computed element in mesh 2 = " << el_id
                << " not matching, read_counter = " << read_counter << "\n";

      if (read_counter > 2) exit(1);
    }

  } else {
    el_id = e2;
  }

  //
  // now that we have found the element, we compute the shape function
  // at the point x and then finally compute displacement at point x
  //
  auto jel_nodes = mesh2->getElementConnectivity(el_id);
  std::vector<double> shapes;
  std::vector<util::Point3> nodes_d;
  std::vector<util::Point3> ys_d;
  for (unsigned long jel_node : jel_nodes) {
    nodes_d.push_back(mesh2->getNode(jel_node));
    ys_d.push_back(sim2->d_y[jel_node]);
  }

  // currently we only support triangle element to compute shape function at
  // any arbitrary point on element
  if (mesh2->getElementType() == util::vtk_type_triangle) {
    auto tri = fe::TriElem(1);
    shapes = tri.getShapes(x1, nodes_d);
  } else {
    std::cerr << "Error: Can not compute shape function at arbitrary point "
                 "for quadrilateral elements\n";
    exit(1);
  }

  auto u = util::Point3();
  for (size_t j = 0; j < jel_nodes.size(); j++)
    for (size_t dof = 0; dof < 3; dof++)
      u[dof] += shapes[j] * (ys_d[j][dof] - nodes_d[j][dof]);

  return u;
}

/*! @brief Computes displacement of coarse mesh at quadrature point of fine mesh
 *
 * @param x1 Node in mesh 1
 * @param e1 Element id in mesh 1
 * @param q Quadrature id
 * @param sim1 Simulation 1 data
 * @param sim2 Simulation 2 data
 * @param dc Collection of user-specified input data
 */
util::Point3 getDisplacementAtQuadPointCurrentSimple(
    const util::Point3 &x1, const size_t &e1, const size_t &q,
    const SimData *sim1, const SimData *sim2, const DataCompare *dc) {
  // get the pointer to mesh 2
  const auto mesh2 = sim2->d_mesh_p;

  //
  // find the node which is closest to y. Use current
  // configuration of node.
  //
  // Note that y here is quadrature point's current position
  //

  size_t el_id = 0;

  static int debug_counter = 0;
  static bool debug_el_ids = false;

  std::string filename = dc->d_pathOut + "/debug_find_elem.txt";
  static FILE *fout = nullptr;

  if (debug_el_ids) fout = fopen(filename.c_str(), "w");

  //
  double h = sim2->d_h;

  // loop over node and find the node which is closest to the point
  int i_found = -1;
  double check = std::sqrt(2.) * 0.5 * h + 0.01 * h;
  double dist = 10.0 * h;
  for (size_t i = 0; i < mesh2->getNumNodes(); i++) {
    auto diff = x1 - (sim2->d_y[i] + mesh2->getNode(i));

    if (util::compare::definitelyLessThan(diff.length(), dist)) {
      dist = diff.length();
      i_found = i;
    }
  }
  if (i_found == -1) {
    std::cerr << "Node closer to x1 = " << x1.printStr()
              << " not found. h = " << h << ", dist = " << dist
              << " check = " << check << ".\n";
    exit(1);
  }

  // return displacement of the node that we found
  return sim2->d_y[i_found] - mesh2->getNode(i_found);
}

//
// compute error
//
void compute(bool read_12, const YAML::Node &config) {
  // create various data
  auto sim1 = SimData();
  auto sim2 = SimData();
  auto dc = DataCompare();

  dc.d_diffAtCurrent = true;

  // read file
  fe::readInputFile(&sim1, &sim2, &dc, config);

  // create output file stream
  FILE *file_out = fopen(dc.d_filenameOut.c_str(), "w");

  // write header
  // fprintf(file_out, "Time    L2_Error    Sup_Error\n");
  if (read_12)
    std::cout << "***   set 1-2   ***\n";
  else
    std::cout << "\n***   set 2-3   ***\n";

  //
  // data for mesh 1 (fine mesh)
  //
  std::vector<std::vector<fe::QuadData>> qd_data1;

  //
  // Data to store elements searched in first call. These elements are in
  // mesh 2 and contain the quadrature points of mesh 1 elements.
  //
  std::vector<std::vector<size_t>> els_cm_of_qpts;

  // loop over time (assuming dt is same for both files)
  size_t read_counter = 1;
  for (size_t k = 0; k <= sim1.d_Nt; k++) {
    double tk1 = double(k) * sim1.d_dt;

    double l2 = 0.;   // to hold l2 norm of error
    double sup = 0.;  // to hold sup norm of error

    size_t dt_check = k % sim1.d_dtOut;

    // in first call, create mesh data for both simulation 1 and 2
    if (dt_check == 0 && read_counter == 1) {
      {
        // mesh for simulation 1
        // do we read mesh from mesh file provided in the input yaml or we read
        // mesh from the first simulation output file?
        if (sim1.d_filename.empty())
          sim1.d_filename = dc.d_path1 + "/output_" +
                            std::to_string(read_counter - 1) + ".vtu";

        // create mesh deck
        auto mdeck = inp::MeshDeck();
        mdeck.d_dim = 2;
        mdeck.d_h = sim1.d_h;
        mdeck.d_filename = sim1.d_filename;
        if (dc.d_isFd)
          mdeck.d_spatialDiscretization = "finite_difference";
        else
          mdeck.d_spatialDiscretization = "finite_element";

        // create mesh
        sim1.d_mesh_p = new fe::Mesh(&mdeck);
      }

      {
        // mesh for simulation 2
        // do we read mesh from mesh file provided in the input yaml or we read
        // mesh from the first simulation output file?
        if (sim2.d_filename.empty())
          sim2.d_filename = dc.d_path2 + "/output_" +
                            std::to_string(read_counter - 1) + ".vtu";

        // create mesh deck
        auto mdeck = inp::MeshDeck();
        mdeck.d_dim = 2;
        mdeck.d_h = sim2.d_h;
        mdeck.d_filename = sim2.d_filename;
        if (dc.d_isFd)
          mdeck.d_spatialDiscretization = "finite_difference";
        else
          mdeck.d_spatialDiscretization = "finite_element";

        // create mesh
        sim2.d_mesh_p = new fe::Mesh(&mdeck);
      }
    }

    // we compare data in this step
    if (dt_check == 0) {
      {
        // read data 1
        std::string filename1 =
            dc.d_path1 + "/output_" + std::to_string(read_counter - 1) + ".vtu";

        // read current position
        rw::reader::readVtuFileNodes(filename1, sim1.d_mesh_p->getDimension(),
                                     &sim1.d_y, false);

        // read displacement (if Tag displacement is not in simulation output
        // file, then displacement can be found by subtracting current
        // position of nodes with their reference position)
        if (!rw::reader::readVtuFilePointData(filename1, "Displacement",
                                              &sim1.d_u)) {
          sim1.d_u.resize(sim1.d_mesh_p->getNumNodes());
          for (size_t i = 0; i < sim1.d_mesh_p->getNumNodes(); i++)
            sim1.d_u[i] = sim1.d_y[i] - sim1.d_mesh_p->getNode(i);
        }
      }

      {
        // read data 1
        std::string filename2 =
            dc.d_path2 + "/output_" + std::to_string(read_counter - 1) + ".vtu";

        // read current position
        rw::reader::readVtuFileNodes(filename2, sim2.d_mesh_p->getDimension(),
                                     &sim2.d_y, false);

        // read displacement (if Tag displacement is not in simulation output
        // file, then displacement can be found by subtracting current
        // position of nodes with their reference position)
        if (!rw::reader::readVtuFilePointData(filename2, "Displacement",
                                              &sim2.d_u)) {
          sim2.d_u.resize(sim2.d_mesh_p->getNumNodes());
          for (size_t i = 0; i < sim2.d_mesh_p->getNumNodes(); i++)
            sim2.d_u[i] = sim2.d_y[i] - sim2.d_mesh_p->getNode(i);
        }
      }

      // resize dummy element list
      if (read_counter == 1 && !dc.d_isFd) {
        // els_fm is of the size of fine mesh
        qd_data1.resize(sim1.d_mesh_p->getNumElements());

        // els_cm_of_qpts is also of the size of fine mesh
        els_cm_of_qpts.resize(sim1.d_mesh_p->getNumElements());

        // allocate enough space to inner vector
        // Max number of quad points is 12 (increase this if needed)
        for (auto &els_cm_of_qpt : els_cm_of_qpts)
          els_cm_of_qpt = std::vector<size_t>(12, 0);
      }

      // get alias for mesh 1 and mesh 2
      const auto mesh1 = sim1.d_mesh_p;
      const auto mesh2 = sim2.d_mesh_p;

      // use hpx loop
      std::vector<double> l2_vec(mesh1->getNumElements(), 0.0);
      std::vector<double> sup_vec(mesh1->getNumElements(), 0.0);

      // loop over nodes in elements1
      for (size_t e = 0; e < mesh1->getNumElements(); e++) {
        auto e_nodes = mesh1->getElementConnectivity(e);

        // compute quad points on copy of current element
        std::vector<fe::QuadData> e_quads;

        if (read_counter == 1) {
          // create node data for quad function
          std::vector<util::Point3> nodes_d;
          for (unsigned long n : e_nodes) {
            nodes_d.push_back(mesh1->getNode(n));
          }

          // call quad function
          if (mesh1->getElementType() == util::vtk_type_triangle) {
            auto tri = fe::TriElem(dc.d_numQuads);
            e_quads = tri.getQuadPoints(nodes_d);
          } else {
            std::cerr << "Error: Only triangle element is supported in data "
                         "comparison"
                      << std::endl;
            exit(1);
          }

          // created temp_el with quad data so now push it to
          // dummy list for reuse next time
          qd_data1[e] = e_quads;

          // resize inner size of els_cm_of_qpts
          if (els_cm_of_qpts[e].size() != e_quads.size())
            els_cm_of_qpts[e] = std::vector<size_t>(e_quads.size(), 0);
        } else
          e_quads = qd_data1[e];

        // loop over quad points
        for (size_t q = 0; q < e_quads.size(); q++) {
          auto xq = e_quads[q].d_p;

          // get displacement and current position of quad point
          auto uq1 = util::Point3();
          auto yq = xq;
          for (size_t i = 0; i < e_nodes.size(); i++) {
            size_t n = e_nodes[i];
            auto xi = mesh1->getNode(n);
            auto yi = sim1.d_y[n];

            for (size_t dof = 0; dof < 3; dof++)
              uq1[dof] += (yi[dof] - xi[dof]) * e_quads[q].d_shapes[i];
          }

          yq += uq1;

          // get element id of coarse mesh
          size_t el_in_cm = els_cm_of_qpts[e][q];

          // get displacement of mesh 2 at xq
          auto uq2 = util::Point3();
          uq2 = getDisplacementAtQuadPoint(xq, e, q, &sim1, &sim2, &dc,
                                           read_counter, el_in_cm);

          // if this is first reading than store element
          // id for reuse in next reading
          if (read_counter == 1) els_cm_of_qpts[e][q] = el_in_cm;

          // compute difference of displacement
          auto du = uq1 - uq2;

          // compute l2 error
          l2 += du.length() * du.length() * e_quads[q].d_w;

          // compute sup norm
          if (util::compare::definitelyGreaterThan(du.length(), sup))
            sup = du.length();
        }
      }

      //
      l2 = std::sqrt(l2);

      // write error to the file
      fprintf(file_out, "%8.6e %8.6e %8.6e\n", tk1, l2, sup);

      // write to the screen
      printf("Time = %8.6e, L2 error =  %8.6e, Sup error = %8.6e\n", tk1, l2,
             sup);

      //
      read_counter++;
    }
  }

  fclose(file_out);
}

}  // namespace fe

//
// main function
//
void dc::fe(YAML::Node config) {
  // check if three sets of data are provided
  bool triple_data = false;
  if (config["Triple_Data"]) triple_data = config["Triple_Data"].as<bool>();

  bool read_12 = true;
  fe::compute(read_12, config);

  read_12 = false;
  if (triple_data) fe::compute(read_12, config);
}