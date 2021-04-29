// Copyright (c)     2017
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

#include <algorithm>
#include <hpx/config.hpp>
#include <hpx/modules/program_options.hpp>
#include <iostream>

#include "dcInclude.h"

int main(int argc, char *argv[]) {
  hpx::program_options::options_description desc("Allowed options");
  desc.add_options()("help", "produce help message")(
      "input-file,i", hpx::program_options::value<std::string>(),
      "Configuration file")("kind,k",
                            hpx::program_options::value<std::string>(),
                            "Kind of simulation");

  hpx::program_options::variables_map vm;
  hpx::program_options::store(
      hpx::program_options::parse_command_line(argc, argv, desc), vm);
  hpx::program_options::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 1;
  }

  std::string filename = "dummy";
  std::string type = "dummy";
  size_t found = 0;

  if (vm.count("input-file")) {
    filename = vm["input-file"].as<std::string>();
    found++;
  }

  if (vm.count("kind")) {
    type = vm["kind"].as<std::string>();
    found++;
  }

  if (filename == "dummy") {
    std::cerr << argv[0] << " -i input.yaml -k type" << std::endl;
    exit(1);
  }

  //
  // read YAML file
  //
  YAML::Node config = YAML::LoadFile(filename);

  bool error = false;

  if (type == "fe") dc::fe(config);

  if (type == "fd") dc::fd(config);

  if (type == "fd_simple") error = dc::fdSimple(config);

  if (error) return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
