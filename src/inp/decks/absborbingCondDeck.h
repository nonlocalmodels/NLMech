////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef INP_ABSORBING_COND_DECK_H
#define INP_ABSORBING_COND_DECK_H

#include <string>
#include "util/utilIO.h"

namespace inp {


/*! Struct for the damping geometry */
struct DampingGeomData {

  /*! @brief Relative location type */
  std::string d_relativeLoc;

  /*! @brief Boolean for checking the x-direction */
  bool d_checkX;
  /*! @brief Boolean for checking the y-direction */
  bool d_checkY;
  /*! @brief Boolean for checking the z-direction */
  bool d_checkZ;

  /*! Thickness of the layer in x-direction */
  double d_layerThicknessX;
  /*! Thickness of the layer in y-direction */
  double d_layerThicknessY;
  /*! Thickness of the layer in z-direction */
  double d_layerThicknessZ;

  /*! Point 1 to define bounding box */
  util::Point3 d_p1;
  /*! Point 2 to define bounding box */
  util::Point3 d_p2;

  DampingGeomData() : d_checkX(false), d_checkY(false), d_checkZ(false),
      d_layerThicknessX(0.), d_layerThicknessY(0.), d_layerThicknessZ(0.),
      d_p1(util::Point3()), d_p2(util::Point3()) {};

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
    oss << tabS << "------- DampingGeomData --------" << std::endl << std::endl;
    oss << tabS << "Relative loc = " << d_relativeLoc << std::endl;
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

/*! @brief Structure to read and store mesh related input data */
struct AbsorbingCondDeck {

  /*!
   * @brief Damping type
   * E.g. "viscous" and "non_viscous"
   */
  bool d_isViscousDamping;

  /*! @brief Damping active */
  bool d_dampingActive;

  /*! @brief Damping coefficient type
   *
   * E.g. "linear", "quadratic", "exponential"
   */
  std::string d_dampingCoeffType;

  /*! @brief Damping coefficient parameters */
  std::vector<double> d_dampingCoeffParams;

  /*! @brief Damping region */
  std::vector<DampingGeomData> d_dampingGeoms;

  /*!
   * @brief Constructor
   */
  AbsorbingCondDeck() : d_isViscousDamping(false), d_dampingActive(false) {};

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
    oss << tabS << "------- AbsorbingCondDeck --------" << std::endl << std::endl;
    oss << tabS << "Is damping active = " << d_dampingActive << std::endl;
    oss << tabS << "Is viscous damping = " << d_isViscousDamping << std::endl;
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

#endif // INP_ABSORBING_COND_DECK_H
