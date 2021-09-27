////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTIL_PARALLEL_H
#define UTIL_PARALELL_H

#include<hpx/include/parallel_copy.hpp>

namespace util {

namespace parallel {

/*
 *  @brief Copies the values in the container in to the container out using HPX's parallel algorithms
 *  @param in The values to be copied
 *  @param out The container the values are copied to
 *  @note The container has to be a container form the c++ standard (std::vector, std::array, or std::list).
 */
template<typename T>
inline void copy(T in, T &out) {

	hpx::for_loop(hpx::execution::par, 0, in.size(),
			[&out, in](boost::uint64_t i) {

				out[i] = in[i];

			});

}

/*
 *  @brief Add the container b in place to the container a
 *  @param a The container data is added to
 *  @param b The container added to the container a
 *  @note The container has to be a container form the c++ standard (std::vector, std::array, or std::list).
 */
template<typename T>
inline void addInplace(T &a, T b) {


	hpx::for_loop(hpx::execution::par, 0, a.size(),
			[&a, b](boost::uint64_t i) {

				a[i] = a[i] + b[i];

			});
}

/*
 *  @brief Subtracts the container b in place to the container a
 *  @param a The container data is added to
 *  @param b The container subtracted to the container a
 *  @note The container has to be a container form the c++ standard (std::vector, std::array, or std::list).
 */
template<typename T>
inline void subInplace(T &a, T b) {

	hpx::for_loop(hpx::execution::par, 0, a.size(),
			[&a, b](boost::uint64_t i) {

				a[i] = a[i] - b[i];

			});
}

}

}

#endif

