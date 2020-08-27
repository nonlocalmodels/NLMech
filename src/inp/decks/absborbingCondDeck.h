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

namespace inp {


/*! Struct for the damping geometry */
struct DampingGeomData {
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

  util::Point3 d_p1;
  util::Point3 d_p2;

  DampingGeomData() : d_checkX(false), d_checkY(false), d_checkZ(false),
      d_layerThicknessX(0.), d_layerThicknessY(0.), d_layerThicknessZ(0.),
      d_p1(util::Point3()), d_p2(util::Point3()) {};
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
};

/** @}*/

} // namespace inp

#endif // INP_ABSORBING_COND_DECK_H
