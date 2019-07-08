// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "fe1D.h"
#include "rw/writer.h"    // definition of vtk and msh writer interface
#include "util/point.h" // definition of Point3
#include "util/feElementDefs.h"   // definition of fe element type
#include <cmath>
#include <iostream>
#include <yaml-cpp/yaml.h> // YAML reader

namespace {

struct InpData {
  std::string d_pathFile;
  std::string d_meshFile;
  std::string d_outFormat;
  std::pair<double, double> d_domain;
  double d_horizon;
  size_t d_r;
  double d_h;
  std::string d_compressType;

  InpData()
      : d_domain(std::make_pair(0., 0.)), d_horizon(0.), d_r(1), d_h(0.){};
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
      data->d_meshFile += "mesh";

    if (config["Output"]["File_Format"])
      data->d_outFormat = config["Output"]["File_Format"].as<std::string>();
    else
      data->d_outFormat = "vtu";

    if (data->d_outFormat == "msh") {
      std::cerr << "Error: For 1-d msh format for mesh output is not "
                   "implemented.\n";
      exit(1);
    }
  }

  if (config["Domain"]) {
    std::vector<double> d;
    for (auto e : config["Domain"])
      d.push_back(e.as<double>());

    data->d_domain.first = d[0];
    data->d_domain.second = d[1];
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

  if (config["Compress_Type"])
    data->d_compressType = config["Compress_Type"].as<std::string>();
}

} // namespace

void tools::mesh::fe1D(const std::string &filename) {

  // read input file
  YAML::Node config = YAML::LoadFile(filename);
  InpData data;
  readInputFile(&data, config);

  //
  size_t num_elems =
      size_t((data.d_domain.second - data.d_domain.first) / data.d_h);
  size_t num_nodes = num_elems + 1;

  // node and element-node connectivity data
  std::vector<util::Point3> nodes(num_nodes, util::Point3());
  std::vector<double> nodal_vols(num_nodes, data.d_h);
  nodal_vols[0] *= 0.5;
  nodal_vols[num_nodes - 1] *= 0.5;
  size_t element_type = util::vtk_type_line;
  std::vector<size_t> en_con(2 * num_elems, 0);

  // create nodes
  for (size_t i = 0; i < num_nodes; i++)
    nodes[i].d_x = data.d_domain.first + i * data.d_h;

  // create element-node connectivity
  for (size_t i = 0; i < num_elems; i++) {
    en_con[2 * i] = i;
    en_con[2 * i + 1] = i + 1;
  }

  // write to vtu file
  auto writer = rw::writer::Writer(data.d_meshFile, data.d_outFormat,
                                   data.d_compressType);
  writer.appendMesh(&nodes, element_type, &en_con);
  writer.appendPointData("Node_Volume", &nodal_vols);
  writer.addTimeStep(0.);
  writer.close();
}