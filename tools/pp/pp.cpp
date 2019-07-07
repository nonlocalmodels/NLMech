// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "src/compute.h"
#include <algorithm>
#include <hpx/hpx_main.hpp>
#include <boost/program_options.hpp>
#include <iostream>

int main(int argc, char *argv[]) {

  boost::program_options::options_description desc("Allowed options");
  desc.add_options()("help", "produce help message")(
      "input-file,i", boost::program_options::value<std::string>(),
      "Configuration file")(
      "kind,k", boost::program_options::value<std::string>(),
      "Dimension");

  boost::program_options::variables_map vm;
  boost::program_options::store(
      boost::program_options::parse_command_line(argc, argv, desc), vm);
  boost::program_options::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 1;
  }

  std::string filename;
  std::string kind;
  size_t found = 0;

  if (vm.count("input-file")) {

    filename = vm["input-file"].as<std::string>();
    found++;
  }

  if (vm.count("kind")) {

    kind = vm["kind"].as<std::string>();
    found++;
  }

  if (found != 1) {
    std::cerr << argv[0] << " -i input.yaml --hpx:threads=n" << std::endl;
    exit(1);
  }

  // call compute
  auto compute = tools::pp::Compute(filename);

  return EXIT_SUCCESS;
}
