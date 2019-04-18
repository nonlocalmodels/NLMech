// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "influenceFn.h"
#include <cmath>

material::pd::ConstInfluenceFn::ConstInfluenceFn(
    const std::vector<double> &params, const size_t &dim)
    : BaseInfluenceFn(), d_a0(0.) {

  d_a0 = params.empty() ? double(dim + 1) : params[0];
}

double material::pd::ConstInfluenceFn::getInfFn(const double &r) {
  return d_a0;
}

double material::pd::ConstInfluenceFn::getMoment(const size_t &i) {
  return d_a0 / double(i + 1);
}

material::pd::LinearInfluenceFn::LinearInfluenceFn(
    const std::vector<double> &params, const size_t &dim)
    : BaseInfluenceFn(), d_a0(0.), d_a1(0.) {

  if (params.empty()) {
    // choose a0, a1 = -a0 such that \int_0^1 J(r) r^d dr = 1
    // and J(r) = a0 (1 - r)
    if (dim == 1) {
      d_a0 = 6.;
      d_a1 = -d_a0;
    } else if (dim == 2) {
      d_a0 = 12.;
      d_a1 = -d_a0;
    } else if (dim == 3) {
      d_a0 = 20.;
      d_a1 = -d_a0;
    }
  } else {
    d_a0 = params[0];
    if (params.size() < 2)
      d_a1 = -d_a0;
    else
      d_a1 = params[1];
  }
}

double material::pd::LinearInfluenceFn::getInfFn(const double &r) {
  return d_a0 + d_a1 * r;
}

double material::pd::LinearInfluenceFn::getMoment(const size_t &i) {
  return (d_a0 / double(i + 1)) - (d_a1 / double(i + 2));
}

material::pd::GaussianInfluenceFn::GaussianInfluenceFn(
    const std::vector<double> &params, const size_t &dim)
    : BaseInfluenceFn(), d_alpha(0.), d_beta(0.) {

  if (params.empty()) {
    // beta = 0.2 (default value)
    // choose alpha such that \int_0^1 J(r) r^d dr = 1
    d_beta = 0.2;
    if (dim == 1)
      d_alpha = 2. / (d_beta * (1. - std::exp(-1. / d_beta)));
    else if (dim == 2)
      d_alpha = (4.0 / d_beta) * 1.0 /
                (std::sqrt(M_PI * d_beta) * std::erf(1.0 / std::sqrt(d_beta)) -
                 2.0 * std::exp(-1.0 / d_beta));
    else if (dim == 3)
      d_alpha = (2.0 / d_beta) * 1.0 /
                (d_beta - (d_beta + 1.) * std::exp(-1.0 / d_beta));
  } else {
    d_alpha = params[0];
    d_beta = params[1];
  }
}

double material::pd::GaussianInfluenceFn::getInfFn(const double &r) {
  return d_alpha * std::exp(-r * r / d_beta);
}

double material::pd::GaussianInfluenceFn::getMoment(const size_t &i) {

  double sq1 = std::sqrt(d_beta);
  double sq2 = std::sqrt(M_PI);
  // M_i = \int_0^1 alpha exp(-r^2/beta) r^i dr

  if (i == 0) {
    // M0 = 0.5 * \alpha (\beta)^(1/2) * (pi)^(1/2) * erf((1/beta)^(1/2))

    return 0.5 * d_alpha * sq1 * sq2 * std::erf(1. / sq1);
  } else if (i == 1) {
    // M1 = 0.5 * \alpha \beta (1 - exp(-1/beta))

    return 0.5 * d_alpha * d_beta * (1. - std::exp(-1. / d_beta));
  } else if (i == 2) {
    // M2 = 0.5 * \alpha (\beta)^(3/2) * [0.5 * (pi)^(1/2) erf((1/beta)^(1/2)
    // ) - (1/beta)^(1/2) * exp(-1/beta) ]

    return 0.5 * d_alpha * d_beta * sq1 *
           (0.5 * sq2 * std::erf(1. / sq1) -
            (1. / sq1) * std::exp(-1. / d_beta));
  } else if (i == 3) {
    // M3 = 0.5 * \alpha (\beta)^(2) * [1 - ((1/beta) + 1) * exp(-1/beta)]

    return 0.5 * d_alpha * d_beta * d_beta *
           (1. - (1. + 1. / d_beta) * std::exp(-1. / d_beta));
  }
}