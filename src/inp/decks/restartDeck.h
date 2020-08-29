////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef INP_RESTARTDECK_H
#define INP_RESTARTDECK_H

#include <string>
#include "util/utilIO.h"

namespace inp {

/**
 * \ingroup Input
 */
/**@{*/

/*! @brief Structure to read and store restart related data input */
struct RestartDeck {

  /*! @brief restart filename */
  std::string d_file;

  /*! @brief Restart time step */
  size_t d_step;

  /*!
   * @brief Constructor
   */
  RestartDeck() : d_step(0){};

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
    oss << tabS << "------- RestartDeck --------" << std::endl << std::endl;
    oss << tabS << "Restart file = " << d_file << std::endl;
    oss << tabS << "Restart step = " << d_step << std::endl;
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

#endif // INP_RESTARTDECK_H
