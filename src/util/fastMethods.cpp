// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "fastMethods.h"
#include "compare.h"
#include <hpx/include/parallel_minmax.hpp>
#include <hpx/include/parallel_reduce.hpp>

static bool compare_point(const util::Point3 &a, const util::Point3 &b) {

  return util::compare::definitelyLessThan(a.length(), b.length());
}

double util::methods::add(const std::vector<double> &data) {

  return hpx::parallel::reduce(hpx::parallel::execution::par, data.begin(),
                               data.end());
}

double util::methods::max(const std::vector<double> &data, size_t *i) {

  auto max_i = hpx::parallel::max_element(hpx::parallel::execution::par,
                                          data.begin(), data.end());

  if (i != nullptr)
    *i = std::distance(data.begin(), max_i);
  return data[std::distance(data.begin(), max_i)];
}

double util::methods::min(const std::vector<double> &data, size_t *i) {

  auto min_i = hpx::parallel::min_element(hpx::parallel::execution::par,
                                          data.begin(), data.end());

  if (i != nullptr)
    *i = std::distance(data.begin(), min_i);
  return data[std::distance(data.begin(), min_i)];
}

float util::methods::add(const std::vector<float> &data) {

  return hpx::parallel::reduce(hpx::parallel::execution::par, data.begin(),
                               data.end());
}

float util::methods::max(const std::vector<float> &data, size_t *i) {

  auto max_i = hpx::parallel::max_element(hpx::parallel::execution::par,
                                          data.begin(), data.end());

  if (i != nullptr)
    *i = std::distance(data.begin(), max_i);
  return data[std::distance(data.begin(), max_i)];
}

float util::methods::min(const std::vector<float> &data, size_t *i) {

  auto min_i = hpx::parallel::min_element(hpx::parallel::execution::par,
                                          data.begin(), data.end());

  if (i != nullptr)
    *i = std::distance(data.begin(), min_i);
  return data[std::distance(data.begin(), min_i)];
}

util::Point3 util::methods::maxLength(const std::vector<util::Point3> &data) {

  auto max_i = hpx::parallel::max_element(
      hpx::parallel::execution::par, data.begin(), data.end(), &compare_point);

  return data[std::distance(data.begin(), max_i)];
}