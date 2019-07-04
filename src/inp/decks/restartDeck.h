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

namespace inp {

/**
 * \ingroup Input
 */
/**@{*/

/*! @brief Structure to read and store restart related data input */
struct RestartDeck {

  /**
   * @name Data members
   */
  /**@{*/

  /*! @brief restart filename */
  std::string d_file;

  /*! @brief Restart time step */
  size_t d_step;

  /** @}*/

  /*!
   * @brief Constructor
   */
  RestartDeck() : d_step(0){};
};

/** @}*/

} // namespace inp

#endif // INP_RESTARTDECK_H
