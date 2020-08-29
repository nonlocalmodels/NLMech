////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////
#ifndef INP_INITIALCONDITIONDECK_H
#define INP_INITIALCONDITIONDECK_H

#include <string>
#include <vector>
#include "util/utilIO.h"

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
    oss << tabS << "------- ICData --------" << std::endl << std::endl;
    oss << tabS << "Initial condition file = " << d_file << std::endl;
    oss << tabS << "Initial condition type = " << d_type << std::endl;
    oss << tabS << "Number of parameters = " << d_params.size() << std::endl;
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
    oss << tabS << "------- InitialConditionDeck --------" << std::endl << std::endl;
    oss << tabS << "Displacement initial condition data " << std::endl;
    oss << d_uICData.printStr(nt+1, lvl);
    oss << tabS << "Velocity initial condition data " << std::endl;
    oss << d_vICData.printStr(nt+1, lvl);
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

#endif // INP_INITIALCONDITIONDECK_H
