////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef RW_UTIL_H
#define RW_UTIL_H

#include "util/point.h"           // definition of Point3
#include "util/matrix.h"           // definition of matrices
#include <vector>

namespace rw {

std::string getDataType(const std::string &data_name);

} // namespace rw

#endif // RW_UTIL_H
