////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <fe/triElem.h>
#include <util/feElementDefs.h>

#include <limits>

#include "fe/mesh.h"
#include "inp/decks/meshDeck.h"
#include "rw/reader.h"
#include "rw/util.h"
#include "rw/writer.h"
#include "srcInclude.h"
#include "util/compare.h"
#include "util/matrix.h"
#include "util/point.h"
#include "util/transfomation.h"
#include "util/utilGeom.h"

static int init = -1;

static bool error = false;

/*! @brief Local namespace */
namespace {

/*! @brief Read input files
 *
 * @param dim Dimension
 * @param filename1 Simulation data 1
 * @param filename2 Simulation data 2
 * @param out_filename Output filename
 * @param print_screen Data should or should not be printed to screen
 * @param compare_tags List of tags data which should be compared
 * @param tolerance
 * @param config YAML input file
 */
void readInputFile(size_t &dim, std::string &in_filename,
                   std::string &out_filename, double &scale_factor,
                   const YAML::Node &config) {

  if (config["Dimension"])
    dim = config["Dimension"].as<size_t>();
  else
    dim = 2;

  in_filename = config["Input_Filename"].as<std::string>();
  out_filename = config["Out_Filename"].as<std::string>();

  if (config["Scale_Factor"])
    scale_factor = config["Scale_Factor"].as<double>();
}

//
// compute error
//
void compute(const YAML::Node &config) {

  size_t dim = 2;
  std::string in_filename;
  std::string out_filename;
  double scale_factor = 1.;

  // read file
  readInputFile(dim, in_filename, out_filename, scale_factor, config);

  //  std::cout << in_filename << ", " << out_filename << ", " << scale_factor
  //            << "\n";

  // read input vtu file
  // nodes current position
  std::vector<util::Point3> nodes_current;
  rw::reader::readVtuFileNodes(in_filename, dim, &nodes_current, false);

  std::vector<util::Point3> nodes_u;
  rw::reader::readVtuFilePointData(in_filename, "Displacement", &nodes_u);

  auto data_tags = rw::reader::readVtuFilePointTags(in_filename);

  // loop over nodes and modify current position
  size_t counter = 0;
  for (const auto &u : nodes_u) {

    auto temp = nodes_current[counter] - u;
    nodes_current[counter] = temp + u * scale_factor;

    counter++;
  }

  // open a new vtu file and write the data
  auto writer = rw::writer::Writer(out_filename, "vtu", "");
  writer.appendNodes(&nodes_current);

  // loop over tag in data_tags, read data from input vtu file, and
  // write the data to output vtu file
  for (const auto &tag : data_tags) {

    // case when data is of double type
    auto type = rw::getDataType(tag);
    if (type == "double") {
      std::vector<double> data;

      // read
      rw::reader::readVtuFilePointData(in_filename, tag, &data);

      // write
      writer.appendPointData(tag, &data);
    } else if (type == "point") {

      std::vector<util::Point3> data;

      // read
      rw::reader::readVtuFilePointData(in_filename, tag, &data);

      // write
      writer.appendPointData(tag, &data);
    }
  }

  // get time from input vtu file
  writer.addTimeStep(0.);
  writer.close();
}

} // namespace

//
// main function
//
void uscale::run(YAML::Node config) { compute(config); }
