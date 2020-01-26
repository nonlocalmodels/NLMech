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
void readInputFile(std::string &in_filename,
                   std::string &out_filename, double &scale_factor,
                   YAML::Node config) {

  in_filename = config["Input_Filename"].as<std::string>();
  out_filename = config["Out_Filename"].as<std::string>();

  if (config["Scale_Factor"]) scale_factor = config["Scale_Factor"].as<double>();  
}

//
// compute error
//
void compute(const YAML::Node &config) {
  
  std::string in_filename;
  std::string out_filename;
  double scale_factor = 1.;

  // read file
  fdSimple::readInputFile(in_filename, out_filename, scale_factor, config);

  // create output file stream
  FILE *file_out = fopen(out_filename.c_str(), "w");

  // write header
  size_t s_counter = 0;

  std::ostringstream oss;
  for (const auto &s : compare_tags) {
    oss << s.c_str() << "_L2_Error, " << s.c_str() << "_Sup_Error";

    // handle special cases
    if (s_counter < compare_tags.size() - 1) oss << ", ";

    s_counter++;
  }
  oss << "\n";

  fprintf(file_out, "%s", oss.str().c_str());
  if (print_screen) {
    std::cout << "------------------------\n";
    std::cout << oss.str();
  }

  // reset oss
  oss.str("");
  oss.clear();

  // read input filename
  // nodes current position
  std::vector<util::Point3> nodes_current;
  rw::reader::readVtuFileNodes(in_filename, dim, &nodes_current, false);
  
  std::vector<util::Point3> nodes_u;
  rw::reader::readVtuFilePointData(in_filename, "Displacement", &nodes_u);

  // loop over nodes and modify current position
  counter = 0;
  for (const auto &u : nodes_u) {

    auto temp = nodes_current[counter] - u;
    nodes_current[counter] = temp + scale_factor * u;

    counter++;
  }

  // open a new vtu file and write the data
    
}

}  // namespace fdSimple

//
// main function
//
bool dc::fdSimple(YAML::Node config) {
  fdSimple::compute(config);

  return error;
}
