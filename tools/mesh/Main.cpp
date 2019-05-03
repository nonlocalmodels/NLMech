// Copyright (c)     2017
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

#include "src/1d/mesh_oned.hpp"
#include "src/2d/mesh_twod.hpp"

#include <algorithm>

#include <boost/program_options.hpp>

int main(int argc, char *argv[]) {

	boost::program_options::options_description desc("Allowed options");
	desc.add_options()("help", "produce help message")("input-file,i",
			boost::program_options::value<std::string>(), "Configuration file");

	boost::program_options::variables_map vm;
	boost::program_options::store(
			boost::program_options::parse_command_line(argc, argv, desc), vm);
	boost::program_options::notify(vm);

	if (vm.count("help")) {
		cout << desc << "\n";
		return 1;
	}

	std::string filename;
	size_t found = 0;

	if (vm.count("input-file")) {

		filename = vm["input-file"].as<std::string>();
		found++;
	}

	if (found != 1) {
		std::cerr << argv[0] << " -i input.yaml" << std::endl;
		exit(1);
	}

	//
	// read YAML file
	//
	YAML::Node config = YAML::LoadFile(filename);

	std::string type = config["Type"].as<std::string>();

	if (type == "fd_twod")
		mesh::fd_twod(config);

	if (type == "fd_fracture_twod")
		mesh::fd_fracture_twod(config);

	if (type == "fe_twod")
		mesh::fe_twod(config);

	if (type == "fd_oned")
		mesh::fd_oned(config);

	if (type == "fe_oned")
		mesh::fe_oned(config);

	return EXIT_SUCCESS;
}
