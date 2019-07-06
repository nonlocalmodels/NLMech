////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef LOADING_LOADING_H
#define LOADING_LOADING_H

#include <string>
#include <vector>

// forward declaration of loading deck
namespace inp {
struct LoadingDeck;
struct BCData;
} // namespace inp

/*!
 * @brief Collection of methods and database related to loading
 *
 * This namespace provides methods and data members specific to application
 * of displacement and force boundary condition and also initial condition.
 *
 * @sa Loading, InitialCondition
 */
namespace loading {

/*!
 * @brief A base class to apply displacement and force boundary condition
 *
 * Base class which provides method and database for application of boundary
 * conditions in the form of displacement or force. Later temperature
 * boundary condition or other type of boundary condition can also be
 * implemented using this base class.
 */
class Loading {

public:
  /*! @brief Constructor */
  Loading() = default;

protected:
  /**
   * @name Internal data
   */
  /**@{*/

  /*! @brief List of displacement bcs */
  std::vector<inp::BCData> d_bcData;

  /*! @brief List of nodal ids on which bc is to be applied */
  std::vector<std::vector<size_t>> d_bcNodes;

  /** @}*/
};

} // namespace loading

#endif // LOADING_LOADING_H
