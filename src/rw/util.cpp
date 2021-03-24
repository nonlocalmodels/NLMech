////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "util.h"

std::vector<std::string> get_double_tags() {
  return {"Damage",       "Damage_Phi", "Damage_Z",  "Node_Volume",
          "Nodal_Volume", "Fixity",     "Work_Done", "Strain_Energy"};
}

std::vector<std::string> get_point_tags() {
  return {"Displacement", "Velocity", "Force", "Force_Density"};
}

std::string rw::getDataType(const std::string &data_name) {
  for (const auto &tag : get_double_tags())
    if (data_name == tag) return "double";

  for (const auto &tag : get_point_tags())
    if (data_name == tag) return "point";

  return "";
}

}