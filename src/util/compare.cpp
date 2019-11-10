////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "compare.h"

#include <cmath>

bool util::compare::approximatelyEqual(const double &a, const double &b) {
  return std::abs(a - b) <=
         ((std::abs(a) < std::abs(b) ? std::abs(b) : std::abs(a)) *
          COMPARE_EPS);
}

bool util::compare::essentiallyEqual(const double &a, const double &b) {
  return std::abs(a - b) <=
         ((std::abs(a) > std::abs(b) ? std::abs(b) : std::abs(a)) *
          COMPARE_EPS);
}

bool util::compare::definitelyGreaterThan(const double &a, const double &b) {
  return (a - b) > ((std::abs(a) < std::abs(b) ? std::abs(b) : std::abs(a)) *
                    COMPARE_EPS);
}

bool util::compare::definitelyLessThan(const double &a, const double &b) {
  return (b - a) > ((std::abs(a) < std::abs(b) ? std::abs(b) : std::abs(a)) *
                    COMPARE_EPS);
}
