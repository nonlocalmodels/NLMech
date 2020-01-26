// Copyright (c)     2017
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

#include <algorithm>
#include <hpx/config.hpp>
#include <boost/program_options.hpp>
#include <iostream>

#include "srcInclude.h"

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

  std::string filename = "dummy";
  size_t found = 0;

  if (vm.count("input-file")) {
    filename = vm["input-file"].as<std::string>();
    found++;
  }

  if (filename = "dummy") {
    std::cerr << argv[0] << " -i input.yaml -k type" << std::endl;
    exit(1);
  }


  //
  // read YAML file
  //
  YAML::Node config = YAML::LoadFile(filename);

  uscale::run(config);

  return EXIT_SUCCESS;
}
