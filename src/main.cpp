////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <Config.h>
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
    std::cerr << argv[0] << " (Version " << MAJOR_VERSION <<
    		"." << MINOR_VERSION << "." << UPDATE_VERSION <<
			") -i input.yaml --hpx:threads=n" << std::endl;
    exit(1);
  }
  // Print program version
  std::cout << argv[0] << " (Version " << MAJOR_VERSION <<
  		"." << MINOR_VERSION << "." << UPDATE_VERSION <<
			")" << std::endl;
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
