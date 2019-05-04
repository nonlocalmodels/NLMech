// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "fe2D.h"
#include "rw/writer.h"  // definition of VtkWriterInterface
#include "util/point.h" // definition of Point3
#include "util/feElementDefs.h"   // definition of fe element type
#include "util/compare.h"         // compare real numbers
#include <cmath>
#include <iostream>
#include <yaml-cpp/yaml.h> // YAML reader

namespace {

struct InpData {
  std::string d_pathFile;
  std::string d_meshFile;
  std::pair<std::vector<double>, std::vector<double>> d_domain;
  double d_horizon;
  size_t d_r;
  double d_h;
  std::string d_meshType;
  bool d_isFd;

  InpData()
      : d_domain(std::make_pair(std::vector<double>(2, 0.),
                                std::vector<double>(2, 0.))),
        d_horizon(0.), d_r(1), d_h(0.), d_isFd(false){};
};

void readInputFile(InpData *data, YAML::Node config) {

  // read output path, filenames
  if (config["Output"]) {
    if (config["Output"]["Path"])
      data->d_pathFile = config["Output"]["Path"].as<std::string>();
    else
      data->d_pathFile = "./";

    data->d_meshFile = data->d_pathFile;
    if (config["Output"]["Mesh"])
      data->d_meshFile += config["Output"]["Mesh"].as<std::string>();
    else
      data->d_meshFile += "mesh.vtu";
  }

  if (config["Domain"]) {
    std::vector<double> d;
    for (auto e : config["Domain"])
      d.push_back(e.as<double>());

    data->d_domain.first.resize(2);
    data->d_domain.first[0] = d[0];
    data->d_domain.first[1] = d[1];

    data->d_domain.second.resize(2);
    data->d_domain.second[0] = d[2];
    data->d_domain.second[1] = d[3];
  } else {
    std::cerr << "Error: Domain data is not provided in mesh input file.\n";
    exit(1);
  }

  if (config["Horizon"])
    data->d_horizon = config["Horizon"].as<double>();
  else {
    std::cerr << "Error: Horizon is not provided in mesh input file.\n";
    exit(1);
  }

  if (config["Horizon_h_Ratio"]) {
    data->d_r = config["Horizon_h_Ratio"].as<size_t>();
    data->d_h = data->d_horizon / double(data->d_r);
  }

  if (config["Mesh_Size"])
    data->d_h = config["Mesh_Size"].as<double>();
  else {
    if (!config["Horizon_h_Ratio"]) {
      std::cerr << "Error: Can not calculate mesh size. Either provide "
                   "Horizon_h_Ratio or Mesh_Size.\n";
      exit(1);
    }
  }

  if (config["Mesh_Type"])
    data->d_meshType = config["Mesh_Type"].as<std::string>();
  else
    data->d_meshType = "Uniform_Triangle";

  if (config["Is_FD"])
    data->d_isFd = config["Is_FD"].as<bool>();
  else
    data->d_isFd = false;
}

} // namespace

void tools::mesh::uniformSquare(const std::string &filename) {

  // read input file
  YAML::Node config = YAML::LoadFile(filename);
  InpData data;
  readInputFile(&data, config);

  // set tolerance
  const double tol = 1.0E-12;

  // modify boundary so that discretization is good
  // modify domain's x end
  size_t ba_h = (data.d_domain.second[0] - data.d_domain.first[0]) / data.d_h;
  auto B = data.d_domain.first[0] + ba_h * data.d_h;
  if (util::compare::definitelyLessThan(B, data.d_domain.second[0]))
    data.d_domain.second[0] = data.d_domain.first[0] + ba_h * data.d_h;

  // modify domain's y end
  ba_h = (data.d_domain.second[1] - data.d_domain.first[1]) / data.d_h;
  B = data.d_domain.first[1] + ba_h * data.d_h;
  if (util::compare::definitelyLessThan(B, data.d_domain.second[1]))
    data.d_domain.second[1] = data.d_domain.first[1] + ba_h * data.d_h;

  // number of divisions in x and y directions and total number of nodes
  size_t nx = (data.d_domain.second[0] - data.d_domain.first[0]) / data.d_h;
  size_t ny = (data.d_domain.second[1] - data.d_domain.first[1]) / data.d_h;
  size_t num_nodes = (nx + 1) * (ny + 1);
  size_t num_elems = nx * ny;

  // local nodal datas
  std::vector<util::Point3> nodes(num_nodes, util::Point3());
  std::vector<double> nodal_vols(num_nodes, data.d_h * data.d_h);

  // create nodal data
  for (size_t j = 0; j <= ny; j++)
    for (size_t i = 0; i <= nx; i++) {
      // node number
      size_t n = j * nx + i;
      nodes[n] =
          util::Point3(data.d_domain.first[0] + double(i) * data.d_h,
                       data.d_domain.first[1] + double(j) * data.d_h, 0.);

      if (i == 0 || i == nx)
        nodal_vols[n] *= 0.5;
      if (j == 0 || j == ny)
        nodal_vols[n] *= 0.5;
    } // loop over j

  // if this mesh is being generated for finite difference simulation
  // we do not require element-node connectivity
  if (data.d_isFd) {

    // write data to vtu file
    auto writer = rw::writer::VtkWriterInterface(data.d_meshFile);
    writer.appendNodes(&nodes);
    writer.appendPointData("Node_Volume", &nodal_vols);
    writer.addTimeStep(0.);
    writer.close();

    return;
  }

  // we reached here means we also compute element-node connectivity
  size_t element_type = util::vtk_type_quad;
  std::vector<size_t> en_con(4 * num_elems, 0);

  // create element-node connectivity data
  for (size_t j = 0; j < ny; j++)
    for (size_t i = 0; i < nx; i++) {
      // element number
      auto n = j * nx + i;

      // element node connectivity (put it in anti clockwise order)
      en_con[4 * n + 0] = j * (nx + 1) + i;
      en_con[4 * n + 1] = j * (nx + 1) + i + 1;
      en_con[4 * n + 2] = (j + 1) * (nx + 1) + i + 1;
      en_con[4 * n + 3] = (j + 1) * (nx + 1) + i;
    }

  // dummy displacement field
  std::vector<util::Point3> u(num_nodes, util::Point3());

  // write to vtu file
  auto writer = rw::writer::VtkWriterInterface(data.d_meshFile);
  writer.appendMesh(&nodes, element_type, &en_con, &u);
  writer.appendPointData("Node_Volume", &nodal_vols);
  writer.addTimeStep(0.);
  writer.close();
}

//
// We create triangular mesh by periodic arrangement of following
// base unit cell consisting of four triangles
//
//
//        n4=n1+nx+1    n5=n4+1       n6=n4+2
//         o-------------o-------------o
//         |    h     /  |  \          |
//         |        /    |    \        |
//    h    |      /      |      \      |
//         |    /        |        \    |
//         |  /          |          \  |
//         o-------------o-------------o
//        n1           n2=n1+1        n3=n1+2
//
// Above is a unit cell of twice of mesh size in x-direction and mesh size in
// y-direction.
//
// We loop over discrete intervals in x and y direction. Suppose we have nx-1
// intervals in x-direction and ny-1 interval in y-direction, and suppose we
// are at (i,j) interval, where 0\leq i \leq nx-1 and 0\leq j \leq ny-1
//
// Then (i,j) has two possibilities:
// 1. It belongs to left half of above unit cell. In this case two triangles
// T1(n1,n5,n4) and T2(n1, n2, n5) are created.
// 2. It belongs to right half of above unit cell. In this case two triangles
// T1(n2, n3, n5) and T2(n3, n6, n5) are created.
//
void tools::mesh::uniformTri(const std::string &filename) {

  // read input file
  YAML::Node config = YAML::LoadFile(filename);
  InpData data;
  readInputFile(&data, config);

  // set tolerance
  const double tol = 1.0E-12;

  // modify boundary so that discretization is good
  // modify domain's x end
  size_t ba_h = (data.d_domain.second[0] - data.d_domain.first[0]) / data.d_h;
  auto B = data.d_domain.first[0] + ba_h * data.d_h;
  if (util::compare::definitelyLessThan(B, data.d_domain.second[0]))
    data.d_domain.second[0] = data.d_domain.first[0] + ba_h * data.d_h;

  // modify domain's y end
  ba_h = (data.d_domain.second[1] - data.d_domain.first[1]) / data.d_h;
  B = data.d_domain.first[1] + ba_h * data.d_h;
  if (util::compare::definitelyLessThan(B, data.d_domain.second[1]))
    data.d_domain.second[1] = data.d_domain.first[1] + ba_h * data.d_h;

  // number of divisions in x and y directions and total number of nodes
  size_t nx = (data.d_domain.second[0] - data.d_domain.first[0]) / data.d_h;
  size_t ny = (data.d_domain.second[1] - data.d_domain.first[1]) / data.d_h;
  size_t num_nodes = (nx + 1) * (ny + 1);

  // local nodal datas
  std::vector<util::Point3> nodes(num_nodes, util::Point3());
  std::vector<double> nodal_vols(num_nodes, data.d_h * data.d_h);

  // create nodal data
  for (size_t j = 0; j <= ny; j++)
    for (size_t i = 0; i <= nx; i++) {
      // node number
      size_t n = j * nx + i;
      nodes[n] =
          util::Point3(data.d_domain.first[0] + double(i) * data.d_h,
                       data.d_domain.first[1] + double(j) * data.d_h, 0.);

      if (i == 0 || i == nx)
        nodal_vols[n] *= 0.5;
      if (j == 0 || j == ny)
        nodal_vols[n] *= 0.5;
    } // loop over j

  // if this mesh is being generated for finite difference simulation
  // we do not require element-node connectivity
  if (data.d_isFd) {

    // write data to vtu file
    auto writer = rw::writer::VtkWriterInterface(data.d_meshFile);
    writer.appendNodes(&nodes);
    writer.appendPointData("Node_Volume", &nodal_vols);
    writer.addTimeStep(0.);
    writer.close();

    return;
  }

  // we reached here means we also compute element-node connectivity
  size_t element_type = util::vtk_type_triangle;
  size_t num_elems = 2 * nx * ny;
  std::vector<size_t> en_con(3 * num_elems, 0);

  // create element-node connectivity
  for (size_t j = 0; j < ny; j++)
    for (size_t i = 0; i < nx; i++) {

      // get node numbers
      auto n1 = j * (nx + 1) + i;
      auto n2 = j * (nx + 1) + i + 1;
      auto n3 = (j + 1) * (nx + 1) + i;
      auto n4 = (j + 1) * (nx + 1) + i + 1;

      // get element numbers
      auto T1 = j * 2 * nx + 2 * i;
      auto T2 = T1 + 1;

      // element-node connectivity
      if (i % 2 == 0) {
        // T1 (anticlockwise order)
        en_con[3 * T1 + 0] = n1;
        en_con[3 * T1 + 1] = n2;
        en_con[3 * T1 + 2] = n4;

        // T2 (anticlockwise order)
        en_con[3 * T2 + 0] = n1;
        en_con[3 * T2 + 1] = n4;
        en_con[3 * T2 + 2] = n3;
      } else {
        // T1 (anticlockwise order)
        en_con[3 * T1 + 0] = n1;
        en_con[3 * T1 + 1] = n2;
        en_con[3 * T1 + 2] = n3;

        // T2 (anticlockwise order)
        en_con[3 * T2 + 0] = n2;
        en_con[3 * T2 + 1] = n4;
        en_con[3 * T2 + 2] = n3;
      }
    } // loop over i

  // dummy displacement field
  std::vector<util::Point3> u(num_nodes, util::Point3());

  // write to vtu file
  auto writer = rw::writer::VtkWriterInterface(data.d_meshFile);
  writer.appendMesh(&nodes, element_type, &en_con, &u);
  writer.appendPointData("Node_Volume", &nodal_vols);
  writer.addTimeStep(0.);
  writer.close();
}

//
// We create symmetric and uniform triangular mesh
//
// Mesh:
//
//         o---------o---------o---------o---------o
//         | \     / | \     / | \     / | \     / |
//         |  \  /   |  \  /   |  \  /   |  \  /   |
//         |    o    |    o    |    o    |    o    |
//         |  /  \   |  /  \   |  /  \   |  /  \   |
//         | /     \ | /     \ | /     \ | /     \ |
//         o---------o---------o---------o---------o
//         | \     / | \     / | \     / | \     / |
//         |  \  /   |  \  /   |  \  /   |  \  /   |
//         |    o    |    o    |    o    |    o    |
//         |  /  \   |  /  \   |  /  \   |  /  \   |
//         | /     \ | /     \ | /     \ | /     \ |
//         o---------o---------o---------o---------o
//         | \     / | \     / | \     / | \     / |
//         |  \  /   |  \  /   |  \  /   |  \  /   |
//         |    o    |    o    |    o    |    o    |
//         |  /12 \  |  /  \   |  /  \   |  /  \   |
//       11| /     \ | /     \ | /     \ | /     \ |
//     --- o---------o---------o---------o---------o
//      |  | \     / | \     / | \     / | \     / |
//      |  |  \  /   |  \  /   |  \  /   |  \  /   |
//     2h  |    o    |    o    |    o    |    o    |
//      |  |  / 2\   |  / 5\   |  / 7\   |  / 9\   |
//      |  | /     \ | /     \ | /     \ | /     \ |
//     --- o---------o---------o---------o---------o
//         1         3         6         8         10
//
//         |--- 2h --|
//
//
// Unit cell of this mesh consists of four triangles
//    n4         n5
//     o---------o
//     | \     / |
//     |  \  /   |
//     |    o n3 |
//     |  /  \   |
//     | /     \ |
//     o---------o
//    n1        n2
//
// n1 - n2 : Distance is 2h
// n1 - n4 : Distance is 2h
//
// 1. We need to create 4 triangles T1(n1,n2,n3), T2(n2,n5,n3), T3(n5,n4,n3), T5
// (n4,n1,n3) for each (i,j) where 0\leq i \leq Lx/2h, 0\leq j \leq Ly/2h
//
// 2. We also need to create node at the center of each unit cell
//
// It can be shown easily that the nodal volume of nodes n1,n2,n3,n4,n5 all are
// 2*h*h when they all are inside sufficiently away from boundary.
//
// For interior node n3, nodal volume is alway 2*h*h
//
// For nodes on the boundary, the correction factor of 0.5 is applied for
// either on x-boundary or y-boundary. If it is in corner then factor 0.25 is
// applied.
//
void tools::mesh::uniformTriSym(const std::string &filename) {

  // read input file
  YAML::Node config = YAML::LoadFile(filename);
  InpData data;
  readInputFile(&data, config);

  // set tolerance
  const double tol = 1.0E-12;

  // modify boundary so that discretization is good
  double h_cell = 2.0 * data.d_h;

  // modify domain's x end
  size_t ba_h = (data.d_domain.second[0] - data.d_domain.first[0]) / h_cell;
  auto B = data.d_domain.first[0] + ba_h * h_cell;
  if (util::compare::definitelyLessThan(B, data.d_domain.second[0]))
    data.d_domain.second[0] = data.d_domain.first[0] + ba_h * h_cell;

  // modify domain's y end
  ba_h = (data.d_domain.second[1] - data.d_domain.first[1]) / h_cell;
  B = data.d_domain.first[1] + ba_h * h_cell;
  if (util::compare::definitelyLessThan(B, data.d_domain.second[1]))
    data.d_domain.second[1] = data.d_domain.first[1] + ba_h * h_cell;

  // number of divisions in x and y directions and total number of nodes
  size_t nx_cell = (data.d_domain.second[0] - data.d_domain.first[0]) / h_cell;
  size_t ny_cell = (data.d_domain.second[1] - data.d_domain.first[1]) / h_cell;
  size_t num_nodes = (nx_cell + 1) * (ny_cell + 1) + nx_cell * ny_cell;

  // local nodal datas
  std::vector<util::Point3> nodes(num_nodes, util::Point3());
  std::vector<double> nodal_vols(num_nodes, 2. * data.d_h * data.d_h);

  int node_counter = 0;
  for (size_t j = 0; j <= ny_cell; j++)
    for (size_t i = 0; i <= nx_cell; i++) {

      // first node (lattice cite / corner of unit cell)

      // numbering above is as follows
      //
      //               n2         n4          n6
      //                o          o          o        <-- these are center
      //                                                   of lattice
      //
      //         o-----------o---------o------         <--- these are lattice
      //        n1          n3         n5                   cites

      // for each j, we have nx+1 lattice cites and nx cites at the center
      // of lattice. nx because for the lattice cite at the right vertical
      // boundary, we do not have center of lattice cite.
      // Special case: When j == ny_cell. In this case we do not have any
      // nodes at the lattice center as we have reached the upper boundary.
      size_t n = j * (nx_cell + 1 + nx_cell) + 2 * i;
      if (j == ny_cell)
        n = j * (nx_cell + 1 + nx_cell) + i;

      // create node
      nodes[n] = util::Point3(data.d_domain.first[0] + double(i) * h_cell,
                              data.d_domain.first[1] + double(j) * h_cell, 0.);

      // modify the volume
      if (i == 0 || i == nx_cell)
        nodal_vols[n] *= 0.5;
      if (j == 0 || j == ny_cell)
        nodal_vols[n] *= 0.5;

      // second node (center of lattice)
      // create only if i is not at the boundary and j is not at the boundary
      if (i < nx_cell && j < ny_cell) {
        n += 1;
        nodes[n] = util::Point3(
            data.d_domain.first[0] + double(i) * h_cell + data.d_h,
            data.d_domain.first[1] + double(j) * h_cell + data.d_h, 0.);
      }
    } // loop over i

  // if this mesh is being generated for finite difference simulation
  // we do not require element-node connectivity
  if (data.d_isFd) {

    // write data to vtu file
    auto writer = rw::writer::VtkWriterInterface(data.d_meshFile);
    writer.appendNodes(&nodes);
    writer.appendPointData("Node_Volume", &nodal_vols);
    writer.addTimeStep(0.);
    writer.close();

    return;
  }

  // we reached here means we also compute element-node connectivity
  size_t element_type = util::vtk_type_triangle;
  size_t num_elems = 4 * nx_cell * ny_cell;
  std::vector<size_t> en_con(3 * num_elems, 0);

  // create element-node connectivity
  for (int j = 0; j < ny_cell; j++)
    for (int i = 0; i < nx_cell; i++) {

      //
      // Each cell consists of 4 triangles and 5 nodes
      //       n4          n3
      //         o---------o
      //         | \ T3  / |
      //         |  \  /   |
      //         |T0  o T2 |
      //         |  /  \   |
      //         | / T1  \ |
      //         o---------o
      //        n0        n2
      //
      // center node: n1
      //
      // size of cell is 2h
      //
      //
      // get node numbers
      std::vector<size_t> ns(5, 0);
      ns[0] = j * (nx_cell + 1 + nx_cell) + 2 * i;
      ns[1] = ns[0] + 1;
      ns[2] = ns[0] + 2;
      // handle special case when j+1 == ny_cell
      if (j + 1 < ny_cell) {
        ns[4] = (j + 1) * (nx_cell + 1 + nx_cell) + 2 * i;
        ns[3] = ns[4] + 2;
      } else {
        ns[4] = (j + 1) * (nx_cell + 1 + nx_cell) + i;
        ns[3] = ns[4] + 1;
      }

      // T1
      auto T = j * 4 * nx_cell + 4 * i + 0;
      en_con[3 * T + 0] = ns[0];
      en_con[3 * T + 1] = ns[1];
      en_con[3 * T + 2] = ns[4];

      // T2
      T = j * 4 * nx_cell + 4 * i + 1;
      en_con[3 * T + 0] = ns[0];
      en_con[3 * T + 1] = ns[2];
      en_con[3 * T + 2] = ns[1];

      // T3
      T = j * 4 * nx_cell + 4 * i + 2;
      en_con[3 * T + 0] = ns[1];
      en_con[3 * T + 1] = ns[2];
      en_con[3 * T + 2] = ns[3];

      // T4
      T = j * 4 * nx_cell + 4 * i + 3;
      en_con[3 * T + 0] = ns[1];
      en_con[3 * T + 1] = ns[3];
      en_con[3 * T + 2] = ns[4];
    } // loop over i

  // dummy displacement field
  std::vector<util::Point3> u(num_nodes, util::Point3());

  // write to vtu file
  auto writer = rw::writer::VtkWriterInterface(data.d_meshFile);
  writer.appendMesh(&nodes, element_type, &en_con, &u);
  writer.appendPointData("Node_Volume", &nodal_vols);
  writer.addTimeStep(0.);
  writer.close();
}

void tools::mesh::fe2D(const std::string &filename) {

  // read input file
  YAML::Node config = YAML::LoadFile(filename);
  InpData data;
  readInputFile(&data, config);

  if (data.d_meshType == "uniform_square")
    tools::mesh::uniformSquare(filename);
  else if (data.d_meshType == "uniform_tri")
    tools::mesh::uniformTri(filename);
  else if (data.d_meshType == "uniform_tri_sym")
    tools::mesh::uniformTriSym(filename);
  else {
    std::cerr << "Error: Check Mesh_Type data. Currently only uniform_square"
                 ", uniform_tri and uniform_tri_sym is implemented.\n";
    exit(1);
  }
}