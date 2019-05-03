// Copyright (c)		2017
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

#ifndef IO_MESH_TWOD_INCLUDE_HPP
#define IO_MESH_TWOD_INCLUDE_HPP

#include <hpx/config.hpp>

#include "../../../../src/includes/IO.hpp"
#include "../../../../src/IO/VtkWriter.hpp"
#include "../../../../src/util/Miscel.hpp" // to check data in config file
#include "../../../../src/external/csv.h" // csv functions
#include "../../../../src/util/Point.hpp"
#include "../../../../src/util/types.hpp"
#include "../../../../src/util/Crack.hpp"

#include <yaml-cpp/yaml.h>

#include <fstream> // to write to csv file


namespace mesh {

void fd_fracture_twod(YAML::Node config);

void fd_twod(YAML::Node config);

void fe_twod(YAML::Node config);

}

#endif