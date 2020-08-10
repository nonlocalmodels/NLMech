////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
///////////////////////////////////////////////////////////////////////////////

#include "DataManager.h"

#include <cassert>
#include <iostream>

#include "fe/mesh.h"
#include "geometry/neighbor.h"
#include "geometry/volumeCorrection.h"
#include "inp/decks/materialDeck.h"
#include "inp/decks/modelDeck.h"
#include "inp/decks/outputDeck.h"
#include "inp/input.h"
#include "loading/fLoading.h"
#include "loading/initialCondition.h"
#include "loading/uLoading.h"
#include "material/pd/ElasticState.h"
#include "material/pdMaterial.h"
#include "util/stateBasedHelperFunctions.h"

data::DataManager::DataManager() {
  // if (instances != 0)
  // std::cerr << "Warning you should have exactly one data manager instance"
  //        << std::endl;
  // instances++;
}

void data::DataManager::setModelDeckP(inp::ModelDeck* pointer) {
  d_modelDeck_p = pointer;
}

inp::ModelDeck* data::DataManager::getModelDeckP() { return d_modelDeck_p; }

void data::DataManager::setOutputDeckP(inp::OutputDeck* pointer) {
  d_outputDeck_p = pointer;
}

inp::OutputDeck* data::DataManager::getOutputDeckP() { return d_outputDeck_p; }

void data::DataManager::setBodyForceP(std::vector<util::Point3>* pointer) {
  d_b_p = pointer;
}

std::vector<util::Point3>* data::DataManager::getBodyForceP() { return d_b_p; }

void data::DataManager::setForceP(std::vector<util::Point3>* pointer) {
  d_f_p = pointer;
}

std::vector<util::Point3>* data::DataManager::getForceP() { return d_f_p; }

void data::DataManager::setVelocityP(std::vector<util::Point3>* pointer) {
  d_v_p = pointer;
}

std::vector<util::Point3>* data::DataManager::getVelocityP() { return d_v_p; }

void data::DataManager::setDisplacementP(std::vector<util::Point3>* pointer) {
  d_u_p = pointer;
}

std::vector<util::Point3>* data::DataManager::getDisplacementP() {
  return d_u_p;
}

void data::DataManager::setMeshP(fe::Mesh* pointer) { d_mesh_p = pointer; }

fe::Mesh* data::DataManager::getMeshP() { return d_mesh_p; }

void data::DataManager::setNeighborP(geometry::Neighbor* pointer) {
  d_neighbor_p = pointer;
}

geometry::Neighbor* data::DataManager::getNeighborP() { return d_neighbor_p; }

void data::DataManager::setVolumeCorrectionP(
    geometry::VolumeCorrection* pointer) {
  d_volumeCorrection_p = pointer;
}

geometry::VolumeCorrection* data::DataManager::getVolumeCorrectionP() {
  return d_volumeCorrection_p;
}

void data::DataManager::setStateBasedHelperFunctionsP(
    util::StateBasedHelperFunctions* pointer) {
  d_sbhelper_p = pointer;
}

util::StateBasedHelperFunctions*
data::DataManager::getStateBasedHelperFunctionsP() {
  return d_sbhelper_p;
}

void data::DataManager::setDisplacementLoadingP(loading::ULoading* pointer) {
  d_uLoading_p = pointer;
}

loading::ULoading* data::DataManager::getDisplacementLoadingP() {
  return d_uLoading_p;
}

void data::DataManager::setForceLoadingP(loading::FLoading* pointer) {
  d_fLoading_p = pointer;
}

loading::FLoading* data::DataManager::getForceLoadingP() {
  return d_fLoading_p;
}

void data::DataManager::setExtensionP(
    std::vector<std::vector<double>>* pointer) {
  d_extension_p = pointer;
}

std::vector<std::vector<double>>* data::DataManager::getExtensionP() {
  return d_extension_p;
}

void data::DataManager::setStrainEnergyP(std::vector<float>* pointer) {
  d_e_p = pointer;
}

std::vector<float>* data::DataManager::getStrainEnergyP() { return d_e_p; }

void data::DataManager::setStressTensorP(std::vector<util::Matrix33>* pointer) {
  d_stress_p = pointer;
}

std::vector<util::Matrix33>* data::DataManager::getStressTensorP() {
  return d_stress_p;
}

void data::DataManager::setStrainTensorP(std::vector<util::Matrix33>* pointer) {
  d_strain_p = pointer;
}

std::vector<util::Matrix33>* data::DataManager::getStrainTensorP() {
  return d_strain_p;
}

void data::DataManager::setDilatationP(std::vector<double>* pointer) {
  d_dilatation_p = pointer;
}

std::vector<double>* data::DataManager::getDilatationP() {
  return d_dilatation_p;
}

void data::DataManager::setReactionForceP(std::vector<util::Point3> *pointer){
    d_reaction_force_p = pointer;

}

	std::vector<util::Point3>* data::DataManager::getReactionForceP(){

    return d_reaction_force_p;

  }

  

