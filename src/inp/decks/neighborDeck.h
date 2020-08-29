////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef INP_NEIGHBORDECK_H
#define INP_NEIGHBORDECK_H

#include "util/utilIO.h"

namespace inp {

/**
 * \ingroup Input
 */
/**@{*/

/*! @brief Structure to read and store neighbor list related input data */
struct NeighborDeck {

  /*! @brief Safety factor for neighbor list calculation */
  double d_safetyFactor;

  /*!
   * @brief Flag to include partially inside nodes in neighbor list
   *
   *
   * Flag which specifies if the partially inside nodes (elements if
   * the discretization is weak finite element) should be included in the
   * neighbor list
   */
  bool d_addPartialElems;

  /*!
   * @brief Constructor
   */
  NeighborDeck() : d_safetyFactor(1.0), d_addPartialElems(false){};

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
    oss << tabS << "------- NeighborDeck --------" << std::endl << std::endl;
    oss << tabS << "Safety factor = " << d_safetyFactor << std::endl;
    oss << tabS << "Add partially inside elements = " << d_addPartialElems
        << std::endl;
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

#endif // INP_NEIGHBORDECK_H
