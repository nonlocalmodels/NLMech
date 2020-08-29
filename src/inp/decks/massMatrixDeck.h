////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef INP_MASSMATRIXDECK_H
#define INP_MASSMATRIXDECK_H

#include <string>
#include "util/utilIO.h"

namespace inp {

/**
 * \ingroup Input
 */
/**@{*/

/*! @brief Structure to read and store mass matrix related input data */
struct MassMatrixDeck {
  /*!
   * @brief Mass matrix approximation type
   *
   * List of allowed values are:
   * - \a exact -- no approximation
   * - \a lumped -- lumping of mass matrix
   */
  std::string d_MApproxType;

  /*!
   * @brief Constructor
   */
  MassMatrixDeck() = default;

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
    oss << tabS << "------- MassMatrixDeck --------" << std::endl << std::endl;
    oss << tabS << "Mass matrix type = " << d_MApproxType << std::endl;
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

#endif // INP_MASSMATRIXDECK_H
