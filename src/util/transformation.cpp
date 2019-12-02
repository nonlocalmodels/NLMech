////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <cmath> // definition of sin, cosine etc

#include "transfomation.h"

std::vector<double>
util::transformation::rotateCW2D(const std::vector<double> &x,
                                 const double &theta) {
  return std::vector<double>{x[0] * std::cos(theta) + x[1] * std::sin(theta),
                             -x[0] * std::sin(theta) + x[1] * std::cos(theta),
                             0.0};
}

util::Point3 util::transformation::rotateCW2D(const util::Point3 &x,
                                              const double &theta) {
  return {x.d_x * std::cos(theta) + x.d_y * std::sin(theta),
          -x.d_x * std::sin(theta) + x.d_y * std::cos(theta), 0.0};
}

std::vector<double>
util::transformation::rotateACW2D(const std::vector<double> &x,
                                  const double &theta) {
  return rotateCW2D(x, -theta);
}

util::Point3 util::transformation::rotateACW2D(const util::Point3 &x,
                                               const double &theta) {
  return rotateCW2D(x, -theta);
}

std::vector<size_t> util::transformation::cyclicOrderACW(const size_t &i,
                                                         const size_t &n) {

  std::vector<size_t> d;

  for (size_t j = 0; j < n; j++)
    d.push_back((i + j) % n);

  return d;
}

std::vector<size_t> util::transformation::cyclicOrderACW(const size_t &i,
                                                         const size_t &j,
                                                         const size_t &n) {

  // data to return
  std::vector<size_t> d;

  // exhaust both possibilities
  if (j == (i + 1) % n) {

    // d.push_back(i%n);
    // d.push_back((i+1)%n);
    // d.push_back((i+2)%n);
    for (size_t k = 0; k < n; k++)
      d.push_back((i + k) % n);

    return d;
  } else if (i == (j + 1) % n) {

    // d.push_back(j%n);
    // d.push_back((j+1)%n);
    // d.push_back((j+2)%n);

    for (size_t k = 0; k < n; k++)
      d.push_back((j + k) % n);
    return d;
  }

  std::cerr << "Error: i and j are not connected in cyclic list.\n";
  std::cerr << "Error: i = " << i << ", j = " << j << ", size n = " << n
            << "\n";
  exit(1);
}

std::vector<size_t> util::transformation::cyclicOrderACW(const size_t &i,
                                                         const size_t &j,
                                                         const size_t &k,
                                                         const size_t &n) {

  // TODO Check this implementation

  // data to return
  std::vector<size_t> d;

  // only case of n = 4 is implemented
  if (n != 4) {

    std::cerr << "Error: currently cyclicOrderACW with i,j,k are implemented "
                 "for n = 4 only.\n";
    exit(1);
  }

  // sort indices in increasing order
  std::vector<size_t> s;
  s.push_back(i);
  s.push_back(j);
  s.push_back(k);

  std::sort(s.begin(), s.end());

  size_t s1 = s[0];
  size_t s2 = s[1];
  size_t s3 = s[2];

  // exhaust all possibilities
  if (s1 == 0 and s2 == 1 and s3 == 2) {

    for (size_t m = 0; m < n; m++)
      d.push_back((i + m) % n);

    return d;
  }

  if (s1 == 1 and s2 == 2 and s3 == 3) {

    for (size_t m = 0; m < n; m++)
      d.push_back((i + m) % n);

    return d;
  }

  if (s1 == 0 and s2 == 2 and s3 == 3) {

    for (size_t m = 0; m < n; m++)
      d.push_back((j + m) % n);

    return d;
  }

  if (s1 == 0 and s2 == 1 and s3 == 3) {

    for (size_t m = 0; m < n; m++)
      d.push_back((k + m) % n);

    return d;
  }

  // if reached here then something is wrong
  std::cerr << "Error: indices i, j, and k are not connected in cyclic list.\n";
  std::cerr << "Error: i = " << i << ", j = " << j << ", k = " << k
            << ", size n = " << n << "\n";
  exit(1);
}
