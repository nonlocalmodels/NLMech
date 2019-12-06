////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

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
#include <fe/triElem.h>
#include <util/feElementDefs.h>

static int init = -1;

/*! @brief Local namespace */
namespace fd {

/*! @brief Data structure to hold simulation data in one place */
struct SimData {
  /*!
   * @brief Constructor
   */
  SimData()
      : d_isFd(true), d_mesh_p(nullptr), d_h(0.), d_Nt(1), d_dt(0.),
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
  if (config["Triple_Data"])
    triple_data = config["Triple_Data"].as<bool>();

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
  if (config["Is_FD"])
    dc->d_isFd = config["Is_FD"].as<bool>();
  else
    dc->d_isFd = true;

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
 * @param x1 Node in mesh 1
 * @param sim1 Simulation 1 data
 * @param sim2 Simulation 2 data
 * @param dc Collection of user-specified input data
 * @param read_counter read counter
 * @param n2 Id of node in mesh 2 near to x1
 * @return u Displacement in mesh 2
 */
util::Point3 getDisplacementOfNode(
    const util::Point3 &x1, const SimData *sim1, const SimData *sim2, const DataCompare *dc,
    const size_t &read_counter, size_t &n2) {

  // get the pointer to mesh 2
  const auto mesh2 = sim2->d_mesh_p;

  //
  // find the node which is closest to y. Use current
  // configuration of node.
  //
  // Note that y here is quadrature point's current position
  //

  long int i_found = -1;

  if (read_counter < 1) {
    static int debug_counter = 0;
    static bool debug_el_ids = false;

    std::string filename = dc->d_pathOut + "/debug_find_elem.txt";
    static FILE *fout = nullptr;

    if (debug_el_ids)
      fout = fopen(filename.c_str(), "w");

    //
    double h = sim2->d_h;

    // loop over node and find the node which is closest to the point

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
      std::cerr << "Node closer to x1 = " << x1.printStr()
                << " not found. h = " << h << ", dist = " << dist
                << " check = " << check << ".\n";
      exit(1);
    }
  } else
    i_found = n2;

  // return displacement of the node that we found
  return sim2->d_y[i_found] - mesh2->getNode(i_found);
}

/*! @brief Computes displacement of coarse mesh at quadrature point of fine mesh
 *
 * @param x1 Node in mesh 1
 * @param sim1 Simulation 1 data
 * @param sim2 Simulation 2 data
 * @param dc Collection of user-specified input data
 */
util::Point3 getDisplacementOfNodeCurrent(
    const util::Point3 &x1, const SimData *sim1, const SimData *sim2, const DataCompare *dc) {

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

  if (debug_el_ids)
    fout = fopen(filename.c_str(), "w");

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
void compute(bool read_12, const YAML::Node& config) {

  // create various data
  auto sim1 = SimData();
  auto sim2 = SimData();
  auto dc = DataCompare();

  dc.d_diffAtCurrent = true;

  // read file
  fd::readInputFile(&sim1, &sim2, &dc, config);

  // create output file stream
  FILE *file_out = fopen(dc.d_filenameOut.c_str(), "w");

  // write header
  // fprintf(file_out, "Time    L2_Error    Sup_Error\n");
  if (read_12)
    std::cout << "***   set 1-2   ***\n";
  else
    std::cout << "\n***   set 2-3   ***\n";

  // store nodes we found in mesh 2 during search for reuse
  std::vector<size_t> nodes_mesh_2;

  // loop over time (assuming dt is same for both files)
  size_t read_counter = 1;
  for (size_t k = 0; k <= sim1.d_Nt; k++) {

    double tk1 = double(k) * sim1.d_dt;

    double l2 = 0.;  // to hold l2 norm of error
    double sup = 0.; // to hold sup norm of error

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
          for (size_t i=0; i<sim1.d_mesh_p->getNumNodes(); i++)
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
          for (size_t i=0; i<sim2.d_mesh_p->getNumNodes(); i++)
            sim2.d_u[i] = sim2.d_y[i] - sim2.d_mesh_p->getNode(i);
        }
      }

      // initialize the saved nodes data
      if (read_counter == 1 && !dc.d_diffAtCurrent)
        nodes_mesh_2 = std::vector<size_t>(sim1.d_mesh_p->getNumNodes(), 0);

      // get alias for mesh 1 and mesh 2
      const auto mesh1 = sim1.d_mesh_p;
      const auto mesh2 = sim2.d_mesh_p;

      // use hpx loop
      std::vector<double> l2_vec(mesh1->getNumNodes(), 0.0);
      std::vector<double> sup_vec(mesh1->getNumNodes(), 0.0);

      // loop over nodes in elements1
      for (size_t i = 0; i < mesh1->getNumNodes(); i++) {
        auto xi = mesh1->getNode(i);
        auto yi = sim1.d_y[i];
        auto ui = sim1.d_u[i];

        // get displacement of mesh 2 at xq
        auto uj = util::Point3();
        if (dc.d_diffAtCurrent)
          uj = getDisplacementOfNodeCurrent(yi, &sim1, &sim2, &dc);
        else {
          auto n2 = nodes_mesh_2[i];
          uj = getDisplacementOfNode(xi, &sim1, &sim2, &dc, read_counter, n2);

          if (read_counter == 1)
            nodes_mesh_2[i] = n2;
        }

        // compute difference of displacement
        auto du = ui - uj;

        // compute l2 error
        l2 += du.length() * du.length() * mesh1->getNodalVolume(i);

        // compute sup norm
        if (util::compare::definitelyGreaterThan(du.length(), sup))
          sup = du.length();
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

} // namespace fe

//
// main function
//
void dc::fd(YAML::Node config) {

  // check if three sets of data are provided
  bool triple_data = false;
  if (config["Triple_Data"])
    triple_data = config["Triple_Data"].as<bool>();

  bool read_12 = true;
  fd::compute(read_12, config);

  read_12 = false;
  if (triple_data)
    fd::compute(read_12, config);
}