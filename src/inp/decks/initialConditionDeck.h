// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef INP_INITIALCONDITIONDECK_H
#define INP_INITIALCONDITIONDECK_H

#include <string>
#include <vector>

namespace inp {

/*! @brief Structure to store initial condition input data */
struct ICData {

  /*! @brief Filename (if any) to read velocity/displacement */
  std::string d_file;

  /*!
   * @brief Tag specifying the formula to compute initial displacement and
   * velocity
   *
   * Example:
   *
   * - \a constant -- for constant function
   */
  std::string d_type;

  /*! @brief List of parameters */
  std::vector<double> d_params;

  /*!
   * @brief Constructor
   */
  ICData() = default;
};

/**
 * \ingroup Input
 */
/**@{*/

/*! @brief Structure to read and store policy data */
struct InitialConditionDeck {

  /*! @brief Initial condition data for displacement */
  inp::ICData d_uICData;

  /*! @brief Initial condition data for velocity */
  inp::ICData d_vICData;

  /*!
   * @brief Constructor
   */
  InitialConditionDeck() : d_uICData(), d_vICData(){};
};

/** @}*/

} // namespace inp

#endif // INP_INITIALCONDITIONDECK_H
