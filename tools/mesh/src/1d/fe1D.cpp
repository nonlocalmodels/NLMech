// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "fe1D.hpp"

namespace {

struct InpData {
  std::string d_meshFile;
  std::string d_uBCLeftFile;
  std::string d_uBCRightFile,
  bool d_active;
  std::pair<double, double> d_domain;
  double d_horizon;
  size_t d_r;
  double d_h;
  size_t d_uICFlag;
  std::vector<double> d_uICParams;
  size_t d_vICFlag;
  std::vector<double> d_vICParams;
};

void readInputFile(InpData *dir_data, InpData *neu_data, YAML::Node config) {

  // read directory path where input files should be created
  std::string path;

  std::string dummy;

  //
  // local variable
  //
  std::vector<double> d;

  // read
  if (config["OutputFile"]["Dirichlet_Filename"]) {
    create_dirichlet_files = true;

    util::methods::checkForData(config, "OutputFile", "Dirichlet_Filename",
                                "Path");
    path = config["OutputFile"]["Dirichlet_Filename"]["Path"].as<std::string>();

    util::methods::checkForData(config, "OutputFile", "Dirichlet_Filename",
                                "Fe_Mesh_File");
    dummy = config["OutputFile"]["Dirichlet_Filename"]["Fe_Mesh_File"]
                .as<std::string>();
    fe_mesh_dir_file.append(path);
    fe_mesh_dir_file.append("/");
    fe_mesh_dir_file.append(dummy);

    util::methods::checkForData(config, "OutputFile", "Dirichlet_Filename",
                                "BC_Displ_ID_Left_File");
    dummy = config["OutputFile"]["Dirichlet_Filename"]["BC_Displ_ID_Left_File"]
                .as<std::string>();
    displ_bc_left_dir_file.append(path);
    displ_bc_left_dir_file.append("/");
    displ_bc_left_dir_file.append(dummy);

    util::methods::checkForData(config, "OutputFile", "Dirichlet_Filename",
                                "BC_Displ_ID_Right_File");
    dummy = config["OutputFile"]["Dirichlet_Filename"]["BC_Displ_ID_Right_File"]
                .as<std::string>();
    displ_bc_right_dir_file.append(path);
    displ_bc_right_dir_file.append("/");
    displ_bc_right_dir_file.append(dummy);
  }

  if (config["OutputFile"]["Neuman_Filename"]) {
    create_neuman_files = true;

    util::methods::checkForData(config, "OutputFile", "Neuman_Filename",
                                "Path");
    path = config["OutputFile"]["Neuman_Filename"]["Path"].as<std::string>();

    util::methods::checkForData(config, "OutputFile", "Neuman_Filename",
                                "Fe_Mesh_File");
    dummy = config["OutputFile"]["Neuman_Filename"]["Fe_Mesh_File"]
                .as<std::string>();
    fe_mesh_neu_file.append(path);
    fe_mesh_neu_file.append("/");
    fe_mesh_neu_file.append(dummy);

    util::methods::checkForData(config, "OutputFile", "Neuman_Filename",
                                "BC_Displ_ID_Left_File");
    dummy = config["OutputFile"]["Neuman_Filename"]["BC_Displ_ID_Left_File"]
                .as<std::string>();
    displ_bc_left_neu_file.append(path);
    displ_bc_left_neu_file.append("/");
    displ_bc_left_neu_file.append(dummy);

    util::methods::checkForData(config, "OutputFile", "Neuman_Filename",
                                "BC_Displ_ID_Right_File");
    dummy = config["OutputFile"]["Neuman_Filename"]["BC_Displ_ID_Right_File"]
                .as<std::string>();
    displ_bc_right_neu_file.append(path);
    displ_bc_right_neu_file.append("/");
    displ_bc_right_neu_file.append(dummy);
  }

  util::methods::checkForData(config, "MeshData");
  if (config["MeshData"]) {
    for (auto e : config["MeshData"]["Domain"])
      d.push_back(e.as<double>());

    Domain.first = d[0];
    Domain.second = d[1];

    Horizon = config["MeshData"]["Horizon"].as<double>();

    Ratio = config["MeshData"]["Horizon_Mesh_Ratio"].as<int>();

    // mesh size
    meshSize = Horizon / (double(Ratio));
  }

  //
  // for fe, test_flag is set to false
  //
  testFlag = false;

  util::methods::checkForData(config, "ICData", "Displacement");
  disp_ic_flag = config["ICData"]["Displacement"]["Flag"].as<size_t>();

  disp_ic_params.clear();
  for (auto e : config["ICData"]["Displacement"]["Parameters"])
    disp_ic_params.push_back(e.as<double>());

  util::methods::checkForData(config, "ICData", "Velocity");
  velo_ic_flag = config["ICData"]["Velocity"]["Flag"].as<size_t>();

  velo_ic_params.clear();
  for (auto e : config["ICData"]["Velocity"]["Parameters"])
    velo_ic_params.push_back(e.as<double>());

  return;
}

//
// apply initial condition
//
void computeIC(double xi, double &ui, double &vi,
               std::pair<double, double> domain, bool test_flag,
               size_t disp_ic_flag, std::vector<double> disp_ic_params,
               size_t velo_ic_flag, std::vector<double> velo_ic_params) {

  if (test_flag == false) {
    // displacement
    if (disp_ic_flag == 0) {
      ui = 0.0;
    } else if (disp_ic_flag == 1) {
      // Gaussian pulse
      // u(x) = a exp[-(x-xc)^2/beta]
      // parameters are as follows
      // center of pulse = xc : first element of displ_ic_params
      // beta: second element of displ_ic_params
      // amplitude = a : third element of displ_ic_params

      double xc = disp_ic_params[0];
      double a = disp_ic_params[2];
      double beta = disp_ic_params[1];
      double x_diff = xc - xi;

      ui = a * std::exp(-x_diff * x_diff / beta);
      // if (counter < 100 and std::abs(xi - 0.5)<=0.001){
      //   std::cout<<"a="<<a<<" beta="<<beta<<" ui="<<ui<<"\n";
      //   counter++;
      // }
    } else if (disp_ic_flag == 2) {
      // Gaussian pulse centered at given position
      // u(x) = a exp[-(x-xc_l)^2/beta] + a exp[-(x-xc_r)^2/beta]
      // parameters are as follows
      // center of left pulse = xc_l : first element of displ_ic_params
      // center of right pulse = xc_r : second element of displ_ic_params
      // beta: third element of displ_ic_params
      // amplitude = a : fourth element of displ_ic_params

      double xc_l = disp_ic_params[0];
      double xc_r = disp_ic_params[1];
      double a = disp_ic_params[3];
      double beta = disp_ic_params[2];

      double x_diff_l = xc_l - xi;
      double x_diff_r = xc_r - xi;

      ui = a * std::exp(-x_diff_l * x_diff_l / beta) +
           a * std::exp(-x_diff_r * x_diff_r / beta);
    } else {
      std::cout << "Mesh: Error in displacement ic flag = " << disp_ic_flag
                << "\n";
      exit(EXIT_FAILURE);
    }

    // velocity
    if (velo_ic_flag == 0) {
      vi = 0.0;
    } else if (velo_ic_flag == 1) {
      // Gaussian pulse
      // v(x) = a exp[-(x-xc)^2/beta]
      // parameters are as follows
      // center of pulse = xc : first element of vel_ic_params
      // beta : second element of vel_ic_params
      // amplitude = a : third element of vel_ic_params

      double xc = velo_ic_params[0];
      double a = velo_ic_params[2];
      double beta = velo_ic_params[1];
      double x_diff = xc - xi;

      vi = a * std::exp(-x_diff * x_diff / beta);
    } else if (velo_ic_flag == 2) {
      // Gaussian pulse centered at given position
      // v(x) = a exp[-(x-xc_l)^2/beta] + a exp[-(x-xc_r)^2/beta]
      // parameters are as follows
      // center of left pulse = xc_l : first element of vel_ic_params
      // center of right pulse = xc_r : second element of vel_ic_params
      // beta : third element of vel_ic_params
      // amplitude = a : fourth element of vel_ic_params

      double xc_l = velo_ic_params[0];
      double xc_r = velo_ic_params[1];
      double a = velo_ic_params[3];
      double beta = velo_ic_params[2];

      double x_diff_l = xc_l - xi;
      double x_diff_r = xc_r - xi;

      vi = a * std::exp(-x_diff_l * x_diff_l / beta) +
           a * std::exp(-x_diff_r * x_diff_r / beta);
    } else {
      std::cout << "Mesh: Error in velocity ic flag = " << velo_ic_flag << "\n";
      exit(EXIT_FAILURE);
    }
  } else {

    return;
  }

  return;
}

} // namespace

void mesh::fe_oned(YAML::Node config) {

  // //  local variables
  // YAML::Node config;
  // std::string config_filename = "mesh_fe_1d.yaml";

  std::pair<double, double> domain;
  double horizon;
  int ratio;
  double mesh_size;
  size_t displ_ic_flag;
  std::vector<double> displ_ic_params;
  size_t vel_ic_flag;
  std::vector<double> vel_ic_params;
  int center[3];

  bool test_flag; // flag if input file is for testing the code

  bool create_dirichlet_files = true;
  std::string fe_mesh_diri_filename;
  std::string displ_bc_left_diri_filename;
  std::string displ_bc_right_diri_filename;

  bool create_neuman_files = true;
  std::string fe_mesh_neum_filename;
  std::string displ_bc_left_neum_filename;
  std::string displ_bc_right_neum_filename;

  // read data from config file
  fe_1d::readDataFile(fe_mesh_diri_filename, displ_bc_left_diri_filename,
                      displ_bc_right_diri_filename, fe_mesh_neum_filename,
                      displ_bc_left_neum_filename, displ_bc_right_neum_filename,
                      create_dirichlet_files, create_neuman_files, domain,
                      horizon, ratio, mesh_size, displ_ic_flag, displ_ic_params,
                      vel_ic_flag, vel_ic_params, test_flag, config);

  //
  //	Given [a,b], h, and e (horizon)
  //
  //-------------------------------------------------------------------------
  //
  //	Left Boundary: [a-e,a]
  //
  //	Nonlocal boundary in left is [a-e, a] and its discretization is
  // { a-e, a-e + h, a-e + 2h, ..., a-e + N_l h}
  // where N_l = int(e/h)
  //
  // Note: If e is not multiple of h then N_l < e/h
  // and therefore a-e + N_l h < a
  //
  //	In case e is not multiple of h, we choose e_l = (int(e/h) + 1) h
  // and define nonlocal boundary as
  // [a - e_l, a]
  //
  // For modified nonlocal boundary following holds
  //
  // 1. [a-e,a] is contained in [a-e_l,a]
  // 2. a - e_l + N_new_l h = a where N_new_l = int(e/h) + 1
  //
  //-------------------------------------------------------------------------
  //
  //	Main domain: [a,b]
  //
  //	Discretization is {a, a+h, a+2h, ..., a + Nh}
  // where N = int( (b-a)/h )
  //
  //	Note: If b-a is not multiple of h, then a + Nh < b
  //
  //	In this case, we modify b slightly
  // b_new = a + (int((b-a)/h) + 1) h
  //
  // New domain is [a,b_new]
  //
  //-------------------------------------------------------------------------
  //
  //	Right boundary: [b_new, b_new + e]
  //
  //	Nonlocal boundary in right is [b_new, b_new + e]
  //	and its discretization is
  // {b_new, b_new + h, b_new + 2h, ..., b_new + N_r h}
  // where N_r = int(e/h)
  //
  // Note: If e is not multiple of h then N_l < e/h
  // and therefore b_new + N_r h < b_new + e
  //
  //	In case e is not multiple of h, we modify e as follows
  // e_r = (N_r + 1) h
  // and define nonlocal boundary as
  // [b_new, b_new + e_r]
  //
  // For modified nonlocal boundary following holds
  //
  // 1. [b_new, b_new + e] is contained in [b_new, b_new + e_r]
  // 2. b_new + N_new_r h = b_new + e_r where N_new_r = int(e/h) + 1
  //
  //-------------------------------------------------------------------------
  double a = domain.first;
  double b = domain.second;
  double bd_l = a - horizon; // boundary point on left
  double bd_r = b + horizon; // boundary point on right

  // check left
  int N_l = int(horizon / mesh_size);
  if (double(N_l) * mesh_size < horizon)
    N_l++;
  bd_l = a - N_l * mesh_size;

  // check middle
  int N_m = int((b - a) / mesh_size);
  if (double(N_m) * mesh_size < b) {
    // modify end point of domain
    b = double(N_m) * mesh_size;

    // modify end point nonlocal boundary on right side
    bd_r = b + horizon;
  }

  // check right
  int N_r = int(horizon / mesh_size);
  if (double(N_r) * mesh_size < horizon)
    N_r++;
  bd_r = b + N_r * mesh_size;

  //
  // summary:
  // Nonlocal left boundary : [bd_l, a]
  // Domain : [a,b]
  // Nonlocal right boundary : [b, bd_r]
  //
  int total_num_elems = int((bd_r - bd_l) / mesh_size);
  int total_num_nodes = total_num_elems + 1;

  // open files for corresponding data
  std::ofstream F_dir_u_l;
  std::ofstream F_dir_u_r;

  if (create_dirichlet_files == true) {
    F_dir_u_l.open(displ_bc_left_diri_filename);
    // header
    F_dir_u_l << "id\n";

    F_dir_u_r.open(displ_bc_right_diri_filename);
    // header
    F_dir_u_r << "id\n";
  }

  // open files for corresponding data Neuman
  std::ofstream F_neu_u_l;
  std::ofstream F_neu_u_r;

  if (create_neuman_files == true) {
    F_neu_u_l.open(displ_bc_left_neum_filename);
    // header
    F_neu_u_l << "id\n";

    F_neu_u_r.open(displ_bc_right_neum_filename);
    // header
    F_neu_u_r << "id\n";
  }

  // node counter
  int n_nodes = 0;

  // node and element data
  std::vector<util::node> nodes;
  std::vector<util::element> elements;
  std::vector<util::point> velocity;
  std::vector<util::point> displacement;

  // create nodes first (also add the id of nodes to
  // File_disp_left and File_disp_right)

  //
  // first we create data corresponding to Dirichlet boundary condition
  //
  double tol = 1.0E-12;
  if (create_dirichlet_files == true) {
    for (int i_node = 0; i_node < total_num_nodes; i_node++) {

      double xi = bd_l + i_node * mesh_size;

      if (std::abs(xi - a) <= tol)
        F_dir_u_l << n_nodes << "\n";
      if (std::abs(xi - b) <= tol)
        F_dir_u_r << n_nodes << "\n";

      // create node data
      util::node n_node = util::node(xi);
      n_node.n = n_nodes;

      // set fixity mask of n_node
      if (xi <= a) {
        n_node.fixity |= NONLOCAL_MASK;
        n_node.fixity |= NONLOCAL_X_BEGIN_MASK;
      }

      if (xi >= b) {
        n_node.fixity |= NONLOCAL_MASK;
        n_node.fixity |= NONLOCAL_X_END_MASK;
      }

      if (std::abs(xi - a) <= tol) {
        n_node.fixity |= SURFACE_MASK;
        n_node.fixity |= SURFACE_X_BEGIN_MASK;
      }

      if (std::abs(xi - b) <= tol) {
        n_node.fixity |= SURFACE_MASK;
        n_node.fixity |= SURFACE_X_END_MASK;
      }

      // compute the initial displacement and velocity
      double ui = 0.;
      double vi = 0.;

      fe_1d::computeIC(xi, ui, vi, domain, test_flag, displ_ic_flag,
                       displ_ic_params, vel_ic_flag, vel_ic_params);

      // add to the node and velocity
      n_node.x[0] = n_node.X[0] + ui;

      util::point p = util::point();
      p[0] = vi;
      n_node.v = p;

      // push the node data to nodes
      nodes.push_back(n_node);

      // increment the node counter
      n_nodes++;
    }

    // now create elements
    int element_counter = 0;
    // another loop for elements
    for (int j = 0; j < total_num_elems; j++) {

      // element number
      // int n = j;
      int n = element_counter;

      // element node connectivity (put it in clockwise order)
      // correct formula for node number
      int n1 = j;
      int n2 = j + 1;

      // location of two nodes
      double x1 = bd_l + double(n1) * mesh_size;
      double x2 = bd_l + double(n2) * mesh_size;

      // create element data
      util::element n_elem = util::element(size_t(n));

      // set fixity data of element

      // element is nonlocal if all of it is inside nonlocal boundary
      if (x2 <= a) {
        n_elem.fixity |= NONLOCAL_MASK;
        n_elem.fixity |= NONLOCAL_X_BEGIN_MASK;
      } else if (x1 >= b) {
        n_elem.fixity |= NONLOCAL_MASK;
        n_elem.fixity |= NONLOCAL_X_END_MASK;
      }

      // update connectivity in element in clockwise order
      // (for one d left-to-right)
      n_elem.nodes.push_back(size_t(n1));
      n_elem.nodes.push_back(size_t(n2));

      // assign element type (e.g. triangle, square, line)
      n_elem.type = VTK_TYPE_LINE;

      // update in node data too
      nodes[n1].elems.push_back(n);
      nodes[n2].elems.push_back(n);

      // write length to element data
      n_elem.len = mesh_size;

      // // call function to compute quadrature points
      // std::vector<std::pair<double, double>> quads;
      // util::methods::getOneDQuadPoints(quad_order,
      //    (*nodes)[n1].X[0],
      //    (*nodes)[n2].X[0],
      //    quads);

      // n_elem.quads.resize(quads.size());
      // for (size_t i=0;i<quads.size();i++) {
      //    n_elem.quads[i].first = quads[i].first;
      //    n_elem.quads[i].second[0] = quads[i].second;
      // }

      //
      // n_elem data is complete so add it to elements vector
      //
      elements.push_back(n_elem);

      // increment the element counter
      element_counter++;
    }

    std::cout << "Total number of nodes (Dirichlet Data) = " << n_nodes << "\n";

    // call VtkWriter to output the data
    IO::VtkWriter writer = IO::VtkWriter();
    writer.open(fe_mesh_diri_filename);
    writer.appendMesh(&nodes, &elements);
    // displacement and velocities is in nodes data
    // writer.appendData("Displacement", displacement);
    // writer.appendData("Velocity", velocity);
    writer.addTimeStep(0.0);
    writer.close();

    // close csv files
    F_dir_u_l.close();
    F_dir_u_r.close();
  }

  //
  // now we create data corresponding to Neuman boundary condition
  //
  if (create_neuman_files == true) {

    total_num_elems = int((b - a) / mesh_size);
    total_num_nodes = total_num_elems + 1;

    n_nodes = 0;

    for (int i_node = 0; i_node < total_num_nodes; i_node++) {

      double xi = a + i_node * mesh_size;

      if (i_node)
        F_neu_u_l << i_node << "\n";
      if (i_node == total_num_nodes)
        F_neu_u_r << i_node << "\n";

      // create node data
      util::node n_node = util::node(xi);
      n_node.n = n_nodes;

      // set fixity mask of n_node
      if (i_node == 0) {
        n_node.fixity |= SURFACE_MASK;
        n_node.fixity |= SURFACE_X_BEGIN_MASK;
      }

      if (i_node == total_num_nodes) {
        n_node.fixity |= SURFACE_MASK;
        n_node.fixity |= SURFACE_X_END_MASK;
      }

      // compute the initial displacement and velocity
      double ui = 0.;
      double vi = 0.;

      fe_1d::computeIC(xi, ui, vi, domain, test_flag, displ_ic_flag,
                       displ_ic_params, vel_ic_flag, vel_ic_params);

      // add to the node and velocity
      n_node.x[0] = n_node.X[0] + ui;

      util::point p = util::point();
      p[0] = vi;
      n_node.v = p;

      // push the node data to nodes
      nodes.push_back(n_node);

      // increment the node counter
      n_nodes++;
    }

    // now create elements
    int element_counter = 0;
    // another loop for elements
    for (int j = 0; j < total_num_elems; j++) {

      // element number
      // int n = j;
      int n = element_counter;

      // element node connectivity (put it in clockwise order)
      // correct formula for node number
      int n1 = j;
      int n2 = j + 1;

      // location of two nodes
      double x1 = a + double(n1) * mesh_size;
      double x2 = a + double(n2) * mesh_size;

      // create element data
      util::element n_elem = util::element(size_t(n));

      // set fixity data of element

      // // element is nonlocal if all of it is inside nonlocal boundary
      // if (x2<= a) {
      //    n_elem.fixity |= NONLOCAL_MASK;
      //    n_elem.fixity |= NONLOCAL_X_BEGIN_MASK;
      // }
      // else if (x1 >= b) {
      //    n_elem.fixity |= NONLOCAL_MASK;
      //    n_elem.fixity |= NONLOCAL_X_END_MASK;
      // }

      // update connectivity in element in clockwise order
      // (for one d left-to-right)
      n_elem.nodes.push_back(size_t(n1));
      n_elem.nodes.push_back(size_t(n2));

      // assign element type (e.g. triangle, square, line)
      n_elem.type = VTK_TYPE_LINE;

      // update in node data too
      nodes[n1].elems.push_back(n);
      nodes[n2].elems.push_back(n);

      // write length to element data
      n_elem.len = mesh_size;

      // // call function to compute quadrature points
      // std::vector<std::pair<double, double>> quads;
      // util::methods::getOneDQuadPoints(quad_order,
      //    (*nodes)[n1].X[0],
      //    (*nodes)[n2].X[0],
      //    quads);

      // n_elem.quads.resize(quads.size());
      // for (size_t i=0;i<quads.size();i++) {
      //    n_elem.quads[i].first = quads[i].first;
      //    n_elem.quads[i].second[0] = quads[i].second;
      // }

      //
      // n_elem data is complete so add it to elements vector
      //
      elements.push_back(n_elem);

      // increment the element counter
      element_counter++;
    }

    std::cout << "Total number of nodes (Neumann Data) = " << n_nodes << "\n";

    // call VtkWriter to output the data
    IO::VtkWriter writer = IO::VtkWriter();
    writer.open(fe_mesh_neum_filename);
    writer.appendMesh(&nodes, &elements);
    // displacement and velocities is in nodes data
    // writer.appendData("Displacement", displacement);
    // writer.appendData("Velocity", velocity);
    writer.addTimeStep(0.0);
    writer.close();

    // close csv files
    F_neu_u_l.close();
    F_neu_u_r.close();
  }
}