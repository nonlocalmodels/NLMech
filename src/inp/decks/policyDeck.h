////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef INP_POLICYDECK_H
#define INP_POLICYDECK_H

#include <vector>
#include "util/utilIO.h"

namespace inp {

/**
 * \ingroup Input
 */
/**@{*/

/*! @brief Structure to read and store policy related input data */
struct PolicyDeck {
  /*!
   * @brief Flag which indicates level of memory control to be enforced
   *
   * Default is 0 which means no control. Max at present is 2 which means as
   * much control as possible.
   */
  int d_memControlFlag;

  /*!
   * @brief Enable post-processing calculation
   *
   * Default is true.
   */
  bool d_enablePostProcessing;

  /*!
   * @brief Constructor
   */
  PolicyDeck() : d_memControlFlag(0), d_enablePostProcessing(true){};

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
    oss << tabS << "------- PolicyDeck --------" << std::endl << std::endl;
    oss << tabS << "Memory control flag = " << d_memControlFlag << std::endl;
    oss << tabS << "Post-processing active = " << d_enablePostProcessing << std::endl;
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

#endif // INP_POLICYDECK_H
