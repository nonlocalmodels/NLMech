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
namespace fdSimple {

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
void readInputFile(size_t &dim, std::string &filename1, std::string &filename2,
                   std::string &out_filename, bool &print_screen,
                   std::vector<std::string> &compare_tags, double &tolerance, YAML::Node config) {
  dim = config["Dimension"].as<size_t>();

  filename1 = config["Filename_1"].as<std::string>();
  filename2 = config["Filename_2"].as<std::string>();

  tolerance = config["Tolerance"].as<double>();


  out_filename = config["Out_Filename"].as<std::string>();
  if (config["Print_Screen"])
    print_screen = config["Print_Screen"].as<bool>();
  else
    print_screen = true;

  compare_tags.clear();
  auto ct = config["Compare_Tags"];
  for (auto f : ct) compare_tags.push_back(f.as<std::string>());

  if (compare_tags.empty()) {
    std::cerr << "Error: Please specify at least one data which should be "
                 "compared.\n";
    exit(1);
  }
}

//
// compute error
//
void compute(const YAML::Node &config) {
  size_t dim;

  std::string filename1;
  std::string filename2;

  std::string out_filename;

  bool print_screen;

  double tolerance = std::numeric_limits<double>::max();

  std::vector<std::string> compare_tags;

  // read file
  fdSimple::readInputFile(dim, filename1, filename2, out_filename, print_screen,
                          compare_tags, tolerance, config);

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

  // read filename 1
  // nodes current position
  std::vector<util::Point3> nodes_current_1;
  rw::reader::readVtuFileNodes(filename1, dim, &nodes_current_1, false);

  // read filename 2
  // nodes current position
  std::vector<util::Point3> nodes_current_2;
  rw::reader::readVtuFileNodes(filename2, dim, &nodes_current_2, false);

  // check if size match
  if (nodes_current_1.size() != nodes_current_2.size()) {
    std::cerr << "Error: Two files have different number of nodes.\n";
    exit(1);
  }

  // loop over tags and compute the error for each tag
  s_counter = 0;
  for (const auto &tag : compare_tags) {
    // implement specific cases

    if (tag == "Displacement") {
      // read displacement from file 1
      std::vector<util::Point3> nodes_u_1;

      if (!rw::reader::readVtuFilePointData(filename1, tag, &nodes_u_1)) {
        std::cerr << "Error: " << tag
                  << " data can not be found in the file ="
                     " "
                  << filename1 << std::endl;
        exit(1);
      }

      // read displacement from file 2
      std::vector<util::Point3> nodes_u_2;

      if (!rw::reader::readVtuFilePointData(filename2, tag, &nodes_u_2)) {
        std::cerr << "Error: " << tag
                  << " data can not be found in the file ="
                     " "
                  << filename2 << std::endl;
        exit(1);
      }

      // compute L2 and sup norm
      double l2 = 0.;
      double sup = 0.;
      for (size_t i = 0; i < nodes_u_1.size(); i++) {
        auto du = nodes_u_1[i] - nodes_u_2[i];
        l2 += du.length() * du.length();

        if (util::compare::definitelyGreaterThan(du.length(), sup))
          sup = du.length();
      }

      l2 = std::sqrt(l2);

      // append data
      oss << l2 << ", " << sup;

    } else if (tag == "Velocity") {
      // read displacement from file 1
      std::vector<util::Point3> nodes_v_1;

      if (!rw::reader::readVtuFilePointData(filename1, tag, &nodes_v_1)) {
        std::cerr << "Error: " << tag
                  << " data can not be found in the file ="
                     " "
                  << filename1 << std::endl;
        exit(1);
      }

      // read displacement from file 2
      std::vector<util::Point3> nodes_v_2;

      if (!rw::reader::readVtuFilePointData(filename2, tag, &nodes_v_2)) {
        std::cerr << "Error: " << tag
                  << " data can not be found in the file ="
                     " "
                  << filename2 << std::endl;
        exit(1);
      }

      // compute L2 and sup norm
      double l2 = 0.;
      double sup = 0.;
      for (size_t i = 0; i < nodes_v_1.size(); i++) {
        auto dv = nodes_v_1[i] - nodes_v_2[i];
        l2 += dv.length() * dv.length();

        if (util::compare::definitelyGreaterThan(dv.length(), sup))
          sup = dv.length();
      }

      l2 = std::sqrt(l2);

      // append data
      oss << l2 << ", " << sup;

    } else if (tag == "Force") {
      // read displacement from file 1
      std::vector<util::Point3> nodes_f_1;

      if (!rw::reader::readVtuFilePointData(filename1, tag, &nodes_f_1)) {
        std::cerr << "Error: " << tag
                  << " data can not be found in the file ="
                     " "
                  << filename1 << std::endl;
        exit(1);
      }

      // read displacement from file 2
      std::vector<util::Point3> nodes_f_2;

      if (!rw::reader::readVtuFilePointData(filename2, tag, &nodes_f_2)) {
        std::cerr << "Error: " << tag
                  << " data can not be found in the file ="
                     " "
                  << filename2 << std::endl;
        exit(1);
      }

      // compute L2 and sup norm
      double l2 = 0.;
      double sup = 0.;
      for (size_t i = 0; i < nodes_f_1.size(); i++) {
        auto df = nodes_f_1[i] - nodes_f_2[i];
        l2 += df.length() * df.length();

        if (util::compare::definitelyGreaterThan(df.length(), sup))
          sup = df.length();
      }

      l2 = std::sqrt(l2);

      // append data
      oss << l2 << ", " << sup;
    } else if (tag == "Strain_Energy") {
      // read strain energy from file 1
      std::vector<double> nodes_e_1;

      if (!rw::reader::readVtuFilePointData(filename1, tag, &nodes_e_1)) {
        std::cerr << "Error: " << tag
                  << " data can not be found in the file ="
                     " "
                  << filename1 << std::endl;
        exit(1);
      }

      // read strain energy from file 2
      std::vector<double> nodes_e_2;

      if (!rw::reader::readVtuFilePointData(filename2, tag, &nodes_e_2)) {
        std::cerr << "Error: " << tag
                  << " data can not be found in the file ="
                     " "
                  << filename2 << std::endl;
        exit(1);
      }

      // compute L2 and sup norm
      double l2 = 0.;
      double sup = 0.;
      for (size_t i = 0; i < nodes_e_1.size(); i++) {
        auto df = nodes_e_1[i] - nodes_e_2[i];
        l2 += df * df;

        if (util::compare::definitelyGreaterThan(df, sup)) sup = df;
      }

      l2 = std::sqrt(l2);

      // append data
      oss << l2 << ", " << sup;

    } else if (tag == "Strain_Tensor") {
      // read strain tensor from file 1
      std::vector<util::Matrix33> nodes_strain_1;

      if (!rw::reader::readVtuFilePointData(filename1, tag, &nodes_strain_1)) {
        std::cerr << "Error: " << tag
                  << " data can not be found in the file ="
                     " "
                  << filename1 << std::endl;
        exit(1);
      }

      // read displacement from file 2
      std::vector<util::Matrix33> nodes_strain_2;

      if (!rw::reader::readVtuFilePointData(filename2, tag, &nodes_strain_2)) {
        std::cerr << "Error: " << tag
                  << " data can not be found in the file ="
                     " "
                  << filename2 << std::endl;
        exit(1);
      }

      // compute L2 and sup norm
      double l2 = 0.;
      double sup = 0.;
      for (size_t i = 0; i < nodes_strain_1.size(); i++) {
        auto df = blaze::sum(nodes_strain_1[i] - nodes_strain_2[i]);
        l2 += df * df;

        if (util::compare::definitelyGreaterThan(df, sup)) sup = df;
      }

      l2 = std::sqrt(l2);

      // append data
      oss << l2 << ", " << sup;

    } else if (tag == "Stress_Tensor") {
      // read strain tensor from file 1
      std::vector<util::Matrix33> nodes_strain_1;

      if (!rw::reader::readVtuFilePointData(filename1, tag, &nodes_strain_1)) {
        std::cerr << "Error: " << tag
                  << " data can not be found in the file ="
                     " "
                  << filename1 << std::endl;
        exit(1);
      }

      // read displacement from file 2
      std::vector<util::Matrix33> nodes_strain_2;

      if (!rw::reader::readVtuFilePointData(filename2, tag, &nodes_strain_2)) {
        std::cerr << "Error: " << tag
                  << " data can not be found in the file ="
                     " "
                  << filename2 << std::endl;
        exit(1);
      }

      // compute L2 and sup norm
      double l2 = 0.;
      double sup = 0.;
      for (size_t i = 0; i < nodes_strain_1.size(); i++) {
        auto df = blaze::sum(nodes_strain_1[i] - nodes_strain_2[i]);
        l2 += df * df;

        if (util::compare::definitelyGreaterThan(df, sup)) sup = df;
      }

      l2 = std::sqrt(l2);

      // append data
      oss << l2 << ", " << sup;

    } else {
      std::cerr << "Error: Comparison for tag = " << tag
                << " has not yet been implemented.\n";
      exit(1);
    }

    // handle special cases
    if (s_counter < compare_tags.size() - 1) oss << ", ";

    s_counter++;
  }
  oss << "\n";

  fprintf(file_out, "%s", oss.str().c_str());
  if (print_screen) std::cout << oss.str();

  fclose(file_out);
}

}  // namespace fdSimple

//
// main function
//
void dc::fdSimple(YAML::Node config) { fdSimple::compute(config); }