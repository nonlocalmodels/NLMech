// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef INP_MATERIALDECK_H
#define INP_MATERIALDECK_H

#include <cmath>
#include <string>
#include <vector>

namespace inp {

/*! @brief Structure for elastic properties and fracture properties */
struct MatData {

  /**
   * @name Elastic material properties
   */
  /**@{*/

  /*! @brief Young's elastic modulus */
  double d_E;

  /*! @brief Shear modulus or Lame second parameter */
  double d_G;

  /*! @brief Bulk modulus */
  double d_K;

  /*! @brief Poisson's ratio */
  double d_nu;

  /*! @brief Lame first parameter */
  double d_lambda;

  /*! @brief Lame second parameter */
  double d_mu;

  /** @}*/

  /**
   * @name Fracture properties
   */
  /**@{*/

  /*! @brief Critical stress intensity factor */
  double d_KIc;

  /*! @brief Critical energy release rate */
  double d_Gc;

  /** @}*/

  /*!
   * @brief Constructor
   */
  MatData()
      : d_E(-1.), d_G(-1.), d_K(-1.), d_nu(-1.), d_lambda(-1.), d_mu(-1.),
        d_KIc(-1.), d_Gc(-1.){};

  /**
   * @name Conversion methods
   */
  /**@{*/

  /*!
   * @brief Compute Poisson's ratio from Lame parameters
   * @param lambda Lame first parameter
   * @param mu Lame second parameter
   * @return nu Poisson's ratio
   */
  double toNu(double lambda, double mu) { return lambda * 0.5 / (lambda + mu); }

  /*!
   * @brief Compute Young's modulus E from Bulk modulus K and Poisson's ratio nu
   * @param K Bulk modulus
   * @param nu Poisson's ratio
   * @return E Young's modulus
   */
  double toE(double K, double nu) { return K * (3. * (1. - 2. * nu)); }

  /*!
   * @brief Compute Bulk modulus K from Young's modulus K and Poisson's ratio nu
   * @param E Young's modulus
   * @param nu Poisson's ratio
   * @return K Bulk modulus
   */
  double toK(double E, double nu) { return E / (3. * (1. - 2. * nu)); }

  /*!
   * @brief Compute Lame first parameter lambda from Young's modulus E
   * and Poisson's ratio nu
   * @param E Young's modulus
   * @param nu Poisson's ratio
   * @return lambda Lame first parameter
   */
  double toLambdaE(double E, double nu) {
    return E * nu / ((1. + nu) * (1. - 2. * nu));
  }

  /*!
   * @brief Compute Lame first parameter lambda from Bulk modulus K and
   * Poisson's ratio nu
   * @param K Bulk modulus
   * @param nu Poisson's ratio
   * @return lambda Lame first parameter
   */
  double toLambdaK(double K, double nu) { return 3. * K * nu / (1. + nu); }

  /*!
   * @brief Compute shear modulus from Young's modulus E and Poisson's ratio nu
   * @param E Young's modulus
   * @param nu Poisson's ratio
   * @return G Shear modulus
   */
  double toGE(double E, double nu) { return E / (2. * (1. + nu)); }

  /*!
   * @brief Compute shear modulus from Bulk modulus K and Poisson's ratio nu
   * @param K Bulk modulus
   * @param nu Poisson's ratio
   * @return G Shear modulus
   */
  double toGK(double K, double nu) {
    return 3. * K * (1. - 2. * nu) / (2. * (1. + nu));
  }

  /*!
 * @brief Compute Young's modulus E from Lame first parameter lambda and
 * Poisson's ratio nu
 * @param lambda Lame first parameter
 * @param nu Poisson's ratio
 * @return E Young's modulus
 */
  double toELambda(double lambda, double nu) {
    return lambda * (1. + nu) * (1. - 2. * nu) / nu;
  }

  /*!
   * @brief Compute critical energy release rate Gc from critical
   * stress-intensity factor KIc, Poisson's ratio nu, and Young's modulus E
   *
   * Below conversion from KIc to Gc assumes **plane-stress** condition. For
   * **plane-stress** condition, we need to modify the Young's modulus \f$
   * E\f$ to \f$ \frac{E}{1 - \nu^2} \f$ where \f$ \nu\f$ is the Poisson's
   * ratio.
   * @param KIc Critical stress-intensity factor
   * @param nu Poisson's ratio
   * @param E Young's modulus
   * @return Gc Critical energy release rate
   */
  double toGc(double KIc, double nu, double E) { return KIc * KIc / E; }

  /*!
   * @brief Compute critical stress-intensity factor KIc from critical energy
   * release rate Gc, Poisson's ratio \f$ nu\f$, and Young's modulus E
   *
   * Below conversion from Gc to KIc assumes **plane-stress** condition. For
   * **plane-stress** condition, we need to modify the Young's modulus \f$
   * E\f$ to \f$ \frac{E}{1 - \nu^2} \f$ where \f$ \nu\f$ is the Poisson's
   * ratio.
   * @param Gc Critical energy release rate
   * @param nu Poisson's ratio
   * @param E Young's modulus
   * @return KIc Critical stress-intensity factor
   */
  double toKIc(double Gc, double nu, double E) { return std::sqrt(Gc * E); }

  /** @}*/
};

/**
 * \ingroup Input
 */
/**@{*/

/*! @brief Structure to read and store material related data */
struct MaterialDeck {

  /**
   * @name Data members
   */
  /**@{*/

  /*!
   * @brief Indicates if the 2-d simulation is of plane-strain type (thick
   * material) or plane-stress type (thin material)
   */
  bool d_isPlaneStrain;

  /*! @brief Material type */
  std::string d_materialType;

  /*! @brief Type of pairwise (bond-based) potential */
  size_t d_bondPotentialType;

  /*! @brief Type of hydrostatic (state-based) potential */
  size_t d_statePotentialType;

  /*! @brief Type of influence function */
  size_t d_influenceFnType;

  /*! @brief List of parameters for pairwise potential */
  std::vector<double> d_bondPotentialParams;

  /*! @brief List of parameters for hydrostatic potential */
  std::vector<double> d_statePotentialParams;

  /*! @brief List of parameters for influence function */
  std::vector<double> d_influenceFnParams;

  /*!
   * @brief Flag for irreversible breaking of bonds.
   *
   * True means bond breaking is irreversible.
   */
  bool d_irreversibleBondBreak;

  /*! @brief Flag for contribution to hydrostatic force from the broken bond */
  bool d_stateContributionFromBrokenBond;

  /*! @brief Factor to check if bond is broken */
  double d_checkScFactor;

  /*! @brief Compute Peridynamic material properties from elastic properties */
  bool d_computeParamsFromElastic;

  /*! @brief List of elastic and fracture properties */
  inp::MatData d_matData;

  /*! @brief Density of material */
  double d_density;

  /** @}*/

  /*!
   * @brief Constructor
   */
  MaterialDeck()
      : d_isPlaneStrain(false), d_bondPotentialType(0), d_statePotentialType(0),
        d_influenceFnType(0), d_irreversibleBondBreak(true),
        d_stateContributionFromBrokenBond(true), d_checkScFactor(1.),
        d_computeParamsFromElastic(true), d_matData(inp::MatData()),
        d_density(1.){};
};

/** @}*/

} // namespace inp

#endif // INP_MATERIALDECK_H
