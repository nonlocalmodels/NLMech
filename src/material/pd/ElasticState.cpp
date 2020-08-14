// Copyright (c) 2019    Patrick Diehl
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "ElasticState.h"

#include <iostream>

#include "data/DataManager.h"
#include "inp/decks/materialDeck.h"
#include "inp/decks/modelDeck.h"
#include "inp/decks/outputDeck.h"
#include "util/compare.h"
#include "util/stateBasedHelperFunctions.h"

material::pd::ElasticState::ElasticState(inp::MaterialDeck *deck,
                                         data::DataManager *dataManager)
    : BaseMaterial(dim, horizon) {
  d_dataManager_p = dataManager;

  horizon = dataManager->getModelDeckP()->d_horizon;

  dim = dataManager->getModelDeckP()->d_dim;

  if (dataManager->getStateBasedHelperFunctionsP() != nullptr)
    delete dataManager->getStateBasedHelperFunctionsP();

  // Compute the PD properties from CCM
  if (deck->d_computeParamsFromElastic)
    computeParameters(deck, dim);
  else {
    std::cerr << "Parameters from classical continuum mechanics are necessary "
                 "to compute the elastic state-based peridynamic properties."
              << std::endl;
    exit(1);
  };
  dataManager->setStateBasedHelperFunctionsP(
      new util::StateBasedHelperFunctions(d_dataManager_p, this->d_factor2D));
  d_deck = deck;
  strainEnergy =
      d_dataManager_p->getOutputDeckP()->isTagInOutput("Strain_Energy");
}

void material::pd::ElasticState::computeParameters(inp::MaterialDeck *deck,
                                                   size_t dim) {
  //
  // Need following elastic and fracture properties
  // 1. E or K
  // 2. Poisson's ratio
  // 3. Gc or KIc
  // For bond-based, Poisson's ratio is fixed to 1/4
  //

  // Check for the 1D case (Only the Young's modulus E is mandatory)

  if (util::compare::definitelyLessThan(deck->d_matData.d_E, 0.) && dim == 1) {
    std::cerr
        << "Error: Require the Young's modulus E"
           " to compute the 1D elastic state-based peridynamic parameters.\n";
    exit(1);
  }

  // Check for the 2D case (We need the shear modulus and the Bulk modulus)

  // Check if Young's moduls and Bulk modulus are provided
  if (util::compare::definitelyGreaterThan(deck->d_matData.d_E, 0.) &&
      util::compare::definitelyGreaterThan(deck->d_matData.d_K, 0.) &&
      dim == 2) {
    std::cout << "Warning: Both Young's modulus E and Bulk modulus K are "
                 "provided.\n";
    std::cout << "Warning: To compute the elastic state-based peridynamic "
                 "parameters, we only require the Bulk modulus K.\n";

    exit(1);
  }

  if (util::compare::approximatelyEqual(deck->d_matData.d_K, -1.) && dim == 2) {
    std::cout << "Warning: To compute the elastic state-based peridynamic "
                 "parameters, we require the Bulk modulus K.\n";
    exit(1);
  }

  if (util::compare::approximatelyEqual(deck->d_matData.d_G, -1.) && dim == 2) {
    std::cout << "Warning: To compute the elastic state-based peridynamic "
                 "parameters, we require the Shear modulus G.\n";
    exit(1);
  }

  deck->d_matData.d_nu =
      deck->d_matData.toNuClassical(deck->d_matData.d_K, deck->d_matData.d_G);

  deck->d_matData.d_mu = deck->d_matData.d_G;

  // Compute the correction for 2D plain strain
  if (deck->d_isPlaneStrain == false and dim == 2)

    this->d_factor2D =
        (2. * deck->d_matData.d_nu - 1.) / (deck->d_matData.d_nu - 1.);

  else {
    this->d_factor2D = 1;
  }

  /* Todo Add damage to elastic model
   if (util::compare::definitelyLessThan(deck->d_matData.d_Gc, 0.) &&
   util::compare::definitelyLessThan(deck->d_matData.d_KIc, 0.)) {
   std::cerr << "Error: Require either critical energy release rate Gc or "
   "critical stress intensity factor KIc to compute the RNP "
   "state-based peridynamic parameters.\n";
   exit(1);
   } else if (util::compare::definitelyGreaterThan(deck->d_matData.d_Gc, 0.) &&
   util::compare::definitelyGreaterThan(deck->d_matData.d_KIc, 0.)) {
   std::cout << "Warning: Both critical energy release rate Gc and critical "
   "stress intensity factor KIc are provided.\n";
   std::cout << "Warning: To compute the RNP state-based peridynamic "
   "parameters, we only require one of those.\n";
   std::cout << "Warning: Selecting critical energy release rate Gc to "
   "compute parameters.\n";
   }

   if (util::compare::definitelyLessThan(deck->d_matData.d_nu, 0.)) {
   std::cerr << "Error: Require Poisson's ratio to compute the RNP "
   "state-based peridynamic parameters.\n";
   exit(1);
   }

   */

  // compute E if not provided or K if not provided
  /*
  if (deck->d_matData.d_E > 0.)
    deck->d_matData.d_K =
        deck->d_matData.toK(deck->d_matData.d_E, deck->d_matData.d_nu);

  if (deck->d_matData.d_K > 0. && deck->d_matData.d_E < 0.)
    deck->d_matData.d_E =
        deck->d_matData.toE(deck->d_matData.d_K, deck->d_matData.d_nu);
*/
  /* Todo Add damage to elastic model
   if (deck->d_matData.d_Gc > 0.)
   deck->d_matData.d_KIc = deck->d_matData.toKIc(
   deck->d_matData.d_Gc, deck->d_matData.d_nu, deck->d_matData.d_E);

   if (deck->d_matData.d_KIc > 0. && deck->d_matData.d_Gc < 0.)
   deck->d_matData.d_Gc = deck->d_matData.toGc(
   deck->d_matData.d_KIc, deck->d_matData.d_nu, deck->d_matData.d_E);

   // compute lame parameter
   deck->d_matData.d_lambda =
   deck->d_matData.toLambdaE(deck->d_matData.d_E, deck->d_matData.d_nu);
   deck->d_matData.d_G =
   deck->d_matData.toGE(deck->d_matData.d_E, deck->d_matData.d_nu);
   deck->d_matData.d_mu = deck->d_matData.d_G;
   */
}

std::pair<util::Point3, double> material::pd::ElasticState::getBondEF(
    size_t i, size_t j) {
  double w = 1;

  double t = 0.;
  double alpha_s = 0.;
  double alpha_d = 0.;
  double e_s = 0.;
  double e_d = 0.;
  double alpha = 0.;
  double t_s = 0.;
  double t_d = 0.;

  util::Point3 Y = ((*d_dataManager_p->getMeshP()->getNodesP())[j] +
                    (*d_dataManager_p->getDisplacementP())[j]) -
                   ((*d_dataManager_p->getMeshP()->getNodesP())[i] +
                    (*d_dataManager_p->getDisplacementP())[i]);
  util::Point3 X = (*d_dataManager_p->getMeshP()->getNodesP())[j] -
                   (*d_dataManager_p->getMeshP()->getNodesP())[i];

  util::Point3 M = Y / Y.length();

  auto it =
      std::find(d_dataManager_p->getNeighborP()->getNeighbors(i).begin(),
                d_dataManager_p->getNeighborP()->getNeighbors(i).end(), j);

  size_t k = std::distance(
      d_dataManager_p->getNeighborP()->getNeighbors(i).begin(), it);

  switch (dim) {
    case 1:

      alpha = d_deck->d_matData.d_E /
              (*d_dataManager_p->getVolumeCorrectionP()->d_weightedVolume_p)[i];
      // Scalar force state
      t = alpha * w * (*d_dataManager_p->getExtensionP())[i][k];
      break;
    case 2:
      // PD material parameter

      if (d_deck->d_isPlaneStrain == false)
        alpha_s = (9. / (*d_dataManager_p->getVolumeCorrectionP()
                              ->d_weightedVolume_p)[i]) *
                  (d_deck->d_matData.d_K +
                   std::pow((d_deck->d_matData.d_nu + 1.) /
                                (2. * d_deck->d_matData.d_nu - 1.),
                            2) *
                       d_deck->d_matData.d_mu / 9.);

      if (d_deck->d_isPlaneStrain == true)
        alpha_s = (9. / (*d_dataManager_p->getVolumeCorrectionP()
                              ->d_weightedVolume_p)[i]) *
                  (d_deck->d_matData.d_K + d_deck->d_matData.d_mu / 9.);

      alpha_d =
          (8. /
           (*d_dataManager_p->getVolumeCorrectionP()->d_weightedVolume_p)[i]) *
          d_deck->d_matData.d_mu;

      // Scalar extension states
      e_s = (*d_dataManager_p->getDilatationP())[i] * X.length() / 3.;
      e_d = (*d_dataManager_p->getExtensionP())[i][k] - e_s;

      // Scalar force states
      t_s = (2. * d_factor2D * alpha_s - (3. - 2. * d_factor2D) * alpha_d) * w *
            e_s / 3.;
      t_d = alpha_d * w * e_d;
      t = t_s + t_d;

      break;
    case 3:

      // PD material parameter
      alpha_s =
          (9. /
           (*d_dataManager_p->getVolumeCorrectionP()->d_weightedVolume_p)[i]) *
          d_deck->d_matData.d_K;
      alpha_d =
          (15. /
           (*d_dataManager_p->getVolumeCorrectionP()->d_weightedVolume_p)[i]) *
          d_deck->d_matData.d_mu;
      // Scalar extension states
      e_s = (*d_dataManager_p->getDilatationP())[i] * X.length() / 3.;
      // Scalar force states
      t_s = alpha_s * w * e_s;
      t_d = alpha_d * w * e_d;
      t = t_s + t_d;

      break;
  }

  double strainE = 0;

  if (strainEnergy) {
    if (dim == 1)

      strainE = 0.5 * alpha * w * (*d_dataManager_p->getExtensionP())[i][k] *
                (*d_dataManager_p->getExtensionP())[i][k] *
                (*d_dataManager_p->getVolumeCorrectionP()
                      ->d_volumeCorrection_p)[i][k] *
                (*d_dataManager_p->getMeshP()->getNodalVolumesP())[j];

    else

      strainE = 0.5 * w * (alpha_s * e_s * e_s + alpha_d * e_d * e_d) *
                (*d_dataManager_p->getVolumeCorrectionP()
                      ->d_volumeCorrection_p)[i][k] *
                (*d_dataManager_p->getMeshP()->getNodalVolumesP())[j];
  }

  return std::make_pair<util::Point3, double>(
      std::move(M * t) * (*d_dataManager_p->getVolumeCorrectionP()
                               ->d_volumeCorrection_p)[i][k],
      std::move(strainE));
}

double material::pd::ElasticState::getSc(size_t i, size_t j) { return 0; }

util::Point3 material::pd::ElasticState::Y_vector_state(size_t i, size_t j) {
  return ((*d_dataManager_p->getMeshP()->getNodesP())[j] +
          (*d_dataManager_p->getDisplacementP())[j]) -
         ((*d_dataManager_p->getMeshP()->getNodesP())[i] +
          (*d_dataManager_p->getDisplacementP())[i]);
}

util::Point3 material::pd::ElasticState::X_vector_state(size_t i, size_t j) {
  return (*d_dataManager_p->getMeshP()->getNodesP())[j] -
         (*d_dataManager_p->getMeshP()->getNodesP())[i];
}

util::Matrix33 material::pd::ElasticState::K_shape_tensor(size_t i) {
  double w = 1;

  util::Matrix33 K = util::Matrix33(0.);

  size_t n = 0;
  for (auto j : d_dataManager_p->getNeighborP()->getNeighbors(i)) {
    util::Point3 X = this->X_vector_state(i, j);

    K +=
        X.toMatrix() * w *
        (*d_dataManager_p->getVolumeCorrectionP()->d_volumeCorrection_p)[i][n] *
        (*d_dataManager_p->getMeshP()->getNodalVolumesP())[j];
  }
  return K;
}

util::Matrix33 material::pd::ElasticState::deformation_gradient(size_t i) {
  double w = 1;

  util::Matrix33 tmp = util::Matrix33(0.);

  size_t n = 0;
  for (auto j : d_dataManager_p->getNeighborP()->getNeighbors(i)) {
    util::Point3 X = this->X_vector_state(i, j);
    util::Point3 Y = this->Y_vector_state(i, j);

    tmp +=
        Y.toMatrix(X) * w *
        (*d_dataManager_p->getVolumeCorrectionP()->d_volumeCorrection_p)[i][n] *
        (*d_dataManager_p->getMeshP()->getNodalVolumesP())[j];

    n++;
  }

  util::Matrixij inv = util::Matrixij(dim, dim, 0.);

  inv =
      blaze::inv(blaze::submatrix(this->K_shape_tensor(i), 0UL, 0UL, dim, dim));

  util::Matrix33 shape_tensor_inv;

  blaze::submatrix(shape_tensor_inv, 0UL, 0UL, dim, dim) = blaze::trans(inv);

  util::Matrix33 deformation = tmp * shape_tensor_inv;

  return deformation;
}

util::Matrix33 material::pd::ElasticState::getStrain(size_t i) {
  util::Matrix33 F = deformation_gradient(i);

  util::Matrix33 FTrans = blaze::trans(F);

  util::Matrix33 identity = util::Matrix33(0.);
  blaze::submatrix(identity, 0UL, 0UL, dim, dim) = util::IdentityMatrix(dim);

  return (F + FTrans) * 0.5 - identity;
}

double material::pd::ElasticState::dirac_delta(util::Point3 x, size_t i,
                                               size_t j, size_t m) {
  double delta = 0.;

  if (util::compare::essentiallyEqual(x.length(), 0))
    delta =
        1. / (*d_dataManager_p->getMeshP()->getNodalVolumesP())[j] *
        (*d_dataManager_p->getVolumeCorrectionP()->d_volumeCorrection_p)[i][m];

  return delta;
}

util::Matrix33 material::pd::ElasticState::K_modulus_tensor(size_t i, size_t j,
                                                            size_t k,
                                                            size_t m) {
  util::Point3 Xj = this->X_vector_state(i, j);
  util::Point3 M = Xj / Xj.length();
  util::Point3 Xk = this->X_vector_state(i, k);

  util::Matrix33 K;

  double alpha = 0.;
  double alpha_s = 0.;
  double alpha_d = 0.;
  double alpha_sb = 0.;

  double w = 1;

  if (dim == 1) {
    alpha = d_deck->d_matData.d_E /
            (*d_dataManager_p->getVolumeCorrectionP()->d_weightedVolume_p)[i];
    K = alpha * w * M.toMatrix() * this->dirac_delta(Xk - Xj, i, k, m);
  }

  if (dim == 2) {
    double Nu = (3. * d_deck->d_matData.d_K - 2. * d_deck->d_matData.d_mu) /
                (2. * (3. * d_deck->d_matData.d_K + d_deck->d_matData.d_mu));
    /*
     double factor2d = 0.;
     if (deck->get2DType() == "Plane_Stress")
     factor2d = (2. * Nu - 1.) / (Nu - 1.);
     if (deck->get2DType() == "Plane_Strain")
     factor2d = 1.;
     */

    // Plane stress
    alpha_s =
        (9. /
         (*d_dataManager_p->getVolumeCorrectionP()->d_weightedVolume_p)[i]) *
        (d_deck->d_matData.d_K +
         std::pow((Nu + 1.) / (2. * Nu - 1.), 2) * d_deck->d_matData.d_mu / 9.);

    // Plane strain
    alpha_d =
        (8. /
         (*d_dataManager_p->getVolumeCorrectionP()->d_weightedVolume_p)[i]) *
        d_deck->d_matData.d_mu;
    alpha_sb =
        (2. * d_factor2D * alpha_s - (3. - 2. * d_factor2D) * alpha_d) / 3.;

    K = ((alpha_sb - alpha_d) /
         (*d_dataManager_p->getVolumeCorrectionP()->d_weightedVolume_p)[i]) *
            w * w * Xj.toMatrix(Xk) +
        alpha_d * w * M.toMatrix() * this->dirac_delta(Xk - Xj, i, k, m);
  }

  if (dim == 3) {
    // PD material parameter
    alpha_s =
        (9. /
         (*d_dataManager_p->getVolumeCorrectionP()->d_weightedVolume_p)[i]) *
        d_deck->d_matData.d_K;
    alpha_d =
        (15. /
         (*d_dataManager_p->getVolumeCorrectionP()->d_weightedVolume_p)[i]) *
        d_deck->d_matData.d_mu;
    K = ((alpha_s - alpha_d) /
         (*d_dataManager_p->getVolumeCorrectionP()->d_weightedVolume_p)[i]) *
            w * w * Xj.toMatrix(Xk) +
        alpha_d * w * M.toMatrix() * this->dirac_delta(Xk - Xj, i, k, m);
  }

  return K;
}

util::Matrix33 material::pd::ElasticState::getStress(size_t i) {
  util::Matrix33 stress = util::Matrix33(0.);

  size_t n = 0;

  for (auto j : d_dataManager_p->getNeighborP()->getNeighbors(i)) {
    util::Point3 Xj = this->X_vector_state(i, j);

    size_t m = 0;

    for (auto k : d_dataManager_p->getNeighborP()->getNeighbors(i)) {
      util::Point3 Xk = this->X_vector_state(i, k);

      double volume = (*d_dataManager_p->getVolumeCorrectionP()
                            ->d_volumeCorrection_p)[i][m] *
                      (*d_dataManager_p->getMeshP()->getNodalVolumesP())[k] *
                      (*d_dataManager_p->getVolumeCorrectionP()
                            ->d_volumeCorrection_p)[i][n] *
                      (*d_dataManager_p->getMeshP()->getNodalVolumesP())[j];

      util::Vector3 res = (this->K_modulus_tensor(i, j, k, m) *
                           (this->getStrain(i) * Xk.toVector()));

      util::Matrix33 tmp = util::Matrix33(0.);

      for (size_t o = 0; o < 3; o++) {
        tmp(o, 0) = res[o] * Xj[0];
        tmp(o, 1) = res[o] * Xj[1];
        tmp(o, 2) = res[o] * Xj[2];
      }
      stress += tmp * volume;

      m++;
    }

    n++;
  }

  return stress;
}

double material::pd::ElasticState::getFactor2D() { return d_factor2D; }

void material::pd::ElasticState::update() {
  if (d_dataManager_p->getStateBasedHelperFunctionsP() != nullptr)
    delete d_dataManager_p->getStateBasedHelperFunctionsP();

  d_dataManager_p->setStateBasedHelperFunctionsP(
      new util::StateBasedHelperFunctions(d_dataManager_p, this->d_factor2D));
}
