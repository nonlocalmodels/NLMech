// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef INITIALCONDITIONDECK_H
#define INITIALCONDITIONDECK_H

#include <vector>
#include <string>

namespace inp {

/*! @brief Structure for initial condition data */
struct ICData {

  /**
   * \defgroup Data members
   */
  /**@{*/

  /*! @brief Filename from which initial velocity/displacement will
   * be read
   */
  std::string d_file;

  /*! @brief Tag which specifies the type of formula to be used in
   * computing initial condition.
   * E.g. d_type = "constant" (for constant function), d_type = "file" (if
   * read from file)
   */
  std::string d_type;

  /*! @brief List of parameters */
  std::vector<double> d_params;

  /** @}*/

  /*!
   * @brief Constructor
   */
  ICData() = default;
};

/*! @brief Structure to read and store policy data */
struct InitialConditionDeck {

  /**
   * \defgroup Data members
   */
  /**@{*/

  /*! @brief Initial condition data for displacement */
  inp::ICData d_uICData;

  /*! @brief Initial condition data for velocity */
  inp::ICData d_vICData;

  /** @}*/

  /*!
   * @brief Constructor
   */
  InitialConditionDeck() : d_uICData(), d_vICData(){};
};
} // namespace inp
#endif // INITIALCONDITIONDECK_H
