////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef INP_QUADRATUREDECK_H
#define INP_QUADRATUREDECK_H


#include "util/utilIO.h"

namespace inp {

/**
 * \ingroup Input
 */
/**@{*/

/*! @brief Structure to read and store quadrature point related input data */
struct QuadratureDeck {

  /*! @brief Order of quadrature point integration approximation */
  size_t d_quadOrder;

  /*!
   * @brief Order of quadrature point integration approximation for
   * mass matrix
   */
  size_t d_quadOrderM;

  /*!
   * @brief Constructor
   */
  QuadratureDeck() : d_quadOrder(0), d_quadOrderM(0){};

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
    oss << tabS << "------- QuadratureDeck --------" << std::endl << std::endl;
    oss << tabS << "Quadrature approximation order = " << d_quadOrder <<
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

#endif // INP_QUADRATUREDECK_H
