////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef INP_INTERIORFLAGSDECK_H
#define INP_INTERIORFLAGSDECK_H

#include "util/utilIO.h"

namespace inp {

/**
 * \ingroup Input
 */
/**@{*/

/*!
 * @brief Structure to read and store inout data for interior flags (**no-fail
 * region**)
 */
struct InteriorFlagsDeck {

  /*! @brief Flag which indicates if no-fail region is active */
  bool d_noFailActive;

  /*! @brief Flag which indicates if we compute the flag inline instead of
   * storing it. This is effective only if d_noFailActive is set to true. */
  bool d_computeAndNotStoreFlag;

  /*! @brief Tolerance to decide if the point is in interior/exterior */
  double d_noFailTol;

  /*! @brief Specify multiple regions in which we set no-fail flag to true */
  std::vector<std::pair<std::string, std::vector<double>>> d_noFailRegions;

  /*!
   * @brief Constructor
   */
  InteriorFlagsDeck()
      : d_noFailActive(false), d_computeAndNotStoreFlag(false),
        d_noFailTol(0.){};

  /*!
   * @brief Returns the string containing information about the instance of
   * the object
   *
   * @param nt Number of tabs to append before each line of string
   * @param lvl Level of information sought (higher level means more
   * information)
   * @return string String containing information about this object
   * */
  std::string printStr(int nt = 0, int lvl = 0) const {
    auto tabS = util::io::getTabS(nt);
    std::ostringstream oss;
    oss << tabS << "------- InteriorFlagsDeck --------" << std::endl << std::endl;
    oss << tabS << "No-fail method active = " << d_noFailActive << std::endl;
    oss << tabS << "Compute flags on-the-fly = " << d_computeAndNotStoreFlag << std::endl;
    oss << tabS << "No-fail region tol = " << d_noFailTol <<
        std::endl;
    oss << tabS << std::endl;

    return oss.str();
  };

  /*!
   * @brief Prints the information about the instance of the object
   *
   * @param nt Number of tabs to append before each line of string
   * @param lvl Level of information sought (higher level means more
   * information)
   * */
  void print(int nt = 0, int lvl = 0) const { std::cout << printStr(nt, lvl); };
};

/** @}*/

} // namespace inp

#endif // INP_INTERIORFLAGSDECK_H
