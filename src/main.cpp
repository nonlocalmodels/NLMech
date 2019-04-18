// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include <hpx/hpx_main.hpp>                     // Need main source file
#include <hpx/util/high_resolution_clock.hpp>
#include "inp/input.h"                          // Input class
#include "model/fd/fDModel.h"                   // Model class
#include <iostream>
#include <boost/program_options.hpp> // program options

int main(int argc, char *argv[]) {

  boost::program_options::options_description desc("Allowed options");
  desc.add_options()("help", "produce help message")(
      "input-file,i", boost::program_options::value<std::string>(),
      "Configuration file");

  boost::program_options::variables_map vm;
  boost::program_options::store(
      boost::program_options::parse_command_line(argc, argv, desc), vm);
  boost::program_options::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 1;
  }

  // read input file
  std::string filename;
  if (vm.count("input-file"))
    filename = vm["input-file"].as<std::string>();

  if (filename.empty()) {
    std::cerr << argv[0] << " -i input.yaml --hpx:threads=n \n";
    exit(1);
  }

  // record current time
  std::uint64_t begin = hpx::util::high_resolution_clock::now();

  // read input data
  auto *deck = new inp::Input(filename);

  // check which model to run
  if (deck->getSpatialDiscretization() == "finite_difference")
    model::FDModel fdModel(deck);

  // get time elapsed
  std::uint64_t end = hpx::util::high_resolution_clock::now();
  double elapsed_secs = double(end - begin) / 1.0e9;

  std::cout << " Time elapsed = " << elapsed_secs << " sec \n";

  return EXIT_SUCCESS;
}