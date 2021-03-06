////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef DC_INCLUDE_HPP
#define DC_INCLUDE_HPP

#include <hpx/config.hpp>
#include <yaml-cpp/yaml.h>

namespace dc {

bool fdSimple(YAML::Node config);

void fd(YAML::Node config);

void fe(YAML::Node config);
}

#endif