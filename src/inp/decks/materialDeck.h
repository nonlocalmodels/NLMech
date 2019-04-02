// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef MATERIALDECK_H
#define MATERIALDECK_H

#include <cmath>
#include <vector>

namespace inp {

/*! @brief Structure for elastic properties and fracture properties */
struct MatData {

  /**
   * \defgroup Elastic material properties
   */
  /**@{*/

  /*! @brief Young's elastic modulus */
  double d_E;

  /*! @brief Shear modulus or Lam'e second parameter */
  double d_G;

  /*! @brief Bulk modulus */
  double d_K;

  /*! @brief Poisson's ratio */
  double d_nu;

  /*! @brief Lam'e first parameter */
  double d_lambda;

  /*! @brief Lam'e second parameter */
  double d_mu;

  /** @}*/

  /**
   * \defgroup Fracture properties
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
      : d_E(0.), d_G(0.), d_K(0.), d_nu(0.), d_lambda(0.), d_mu(0.), d_KIc(0.),
        d_Gc(0.){};

  /**
   * \defgroup Methods to compute elastic property using given elastic
   * properties
   */
  /**@{*/

  /*!
   * @brief Compute Young's modulus E from Bulk modulus K and Poisson's ratio nu
   * @param K Bulk modulus
   * @param nu Poisson's ratio
   * @return E Young's modulus
   */
  double toE(double K, double nu) { return K * (3.0 * (1 - 2.0 * nu)); }

  /*!
   * @brief Compute Bulk modulus K from Young's modulus K and Poisson's ratio nu
   * @param E Young's modulus
   * @param nu Poisson's ratio
   * @return K Bulk modulus
   */
  double toK(double E, double nu) { return E / (3.0 * (1 - 2.0 * nu)); }

  /*!
   * @brief Compute Lam'e first parameter lambda from Young's modulus K and
   * Poisson's ratio nu
   * @param E Young's modulus
   * @param nu Poisson's ratio
   * @return lambda Lam'e first parameter
   */
  double toLambdaE(double E, double nu) {
    return E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  }

  /*!
   * @brief Compute Lam'e first parameter lambda from Bulk modulus K and
   * Poisson's ratio nu
   * @param K Bulk modulus
   * @param nu Poisson's ratio
   * @return lambda Lam'e first parameter
   */
  double toLambdaK(double K, double nu) { return 3.0 * K * nu / (1.0 + nu); }

  /*!
   * @brief Compute critical energy release rate Gc from critical
   * stress-intensity factor KIc, Poisson's ratio nu, and Young's modulus E
   * @param KIc Critical stress-intensity factor
   * @param nu Poisson's ratio
   * @param E Young's modulus
   * @return Gc Critical energy release rate
   */
  double toGc(double KIc, double nu, double E) {
    return KIc * KIc * (1. - nu * nu) / E;
  }

  /*!
   * @brief Compute critical stress-intensity factor KIc from critical energy
   * release rate Gc, Poisson's ratio
   * nu, and Young's modulus E
   * @param Gc Critical energy release rate
   * @param nu Poisson's ratio
   * @param E Young's modulus
   * @return KIc Critical stress-intensity factor
   */
  double toKIc(double Gc, double nu, double E) {
    return std::sqrt(Gc * E / (1.0 - nu * nu));
  }

  /** @}*/
};

/*! @brief Structure to read and store material related data */
struct MaterialDeck {

  /**
   * \defgroup Data members
   */
  /**@{*/

  /*! @brief 2D type, either plane-stress (thin material) or plane-strain
   * (very thick material)
   */
  std::string d_2DType;

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

  /*! @brief Flag for irreversible breaking of bonds. True means bond
   * breaking is irreversible
   */
  bool d_irreversibleBondBreak;

  /*! @brief Flag for contribution to hydrostatic force from the broken bond */
  bool d_stateContributionFromBrokenBond;

  /*! @brief Factor to check if bond is broken
   */
  double d_checkScFactor;

  /*! @brief Specify if parameters should be computed from provided elastic
   * material properties
   */
  bool d_computeParamsFromElastic;

  /*! @brief List of elastic and fracture properties */
  inp::MatData d_matData;

  /** @}*/

  /*!
   * @brief Constructor
   */
  MaterialDeck()
      : d_bondPotentialType(0), d_statePotentialType(0), d_influenceFnType(0),
        d_irreversibleBondBreak(true), d_stateContributionFromBrokenBond(true),
        d_checkScFactor(1.), d_computeParamsFromElastic(true), d_matData(){};
};

} // namespace inp
#endif // MATERIALDECK_H
