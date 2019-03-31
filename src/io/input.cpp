// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)
#include <iostream>
#include <yaml-cpp/yaml.h>

#include "decks/fractureDeck.h"
#include "decks/geometryDeck.h"
#include "decks/initialConditionDeck.h"
#include "decks/interiorFlagsDeck.h"
#include "decks/loadingDeck.h"
#include "decks/materialDeck.h"
#include "decks/neighborDeck.h"
#include "decks/outputDeck.h"
#include "decks/policyDeck.h"
#include "decks/problemDeck.h"
#include "decks/solverDeck.h"
#include "input.h"

io::Input::Input(const std::string& filename)
    : d_fractureDeck_p(nullptr), d_geometryDeck_p(nullptr),
      d_initialConditionDeck_p(nullptr), d_interiorFlagsDeck_p(nullptr),
      d_loadingDeck_p(nullptr), d_materialDeck_p(nullptr),
      d_neighborDeck_p(nullptr), d_outputDeck_p(nullptr),
      d_policyDeck_p(nullptr), d_problemDeck_p(nullptr),
      d_solverDeck_p(nullptr) {

  d_inputFilename = filename;

  // follow the order of reading
  setGeometryDeck();
  setProblemDeck();
  setNeighborDeck();
  setFractureDeck();
  setInteriorFlagsDeck();
  setInitialConditionDeck();
  setLoadingDeck();
  setMaterialDeck();
  setOutputDeck();
  setPolicyDeck();
  setSolverDeck();
}

io::FractureDeck *io::Input::getFractureDeck() { return d_fractureDeck_p; }

io::GeometryDeck *io::Input::getGeometryDeck() { return d_geometryDeck_p; }

io::InitialConditionDeck *io::Input::getInitialConditionDeck() {
  return d_initialConditionDeck_p;
}

io::InteriorFlagsDeck *io::Input::getInteriorFlagsDeck() {
  return d_interiorFlagsDeck_p;
}

io::LoadingDeck *io::Input::getLoadingDeck() { return d_loadingDeck_p; }

io::MaterialDeck *io::Input::getMaterialDeck() { return d_materialDeck_p; }

io::NeighborDeck *io::Input::getNeighborDeck() { return d_neighborDeck_p; }

io::OutputDeck *io::Input::getOutputDeck() { return d_outputDeck_p; }

io::PolicyDeck *io::Input::getPolicyDeck() { return d_policyDeck_p; }

io::ProblemDeck *io::Input::getProblemDeck() { return d_problemDeck_p; }

io::SolverDeck *io::Input::getSolverDeck() { return d_solverDeck_p; }

void io::Input::setGeometryDeck() { d_geometryDeck_p = new io::GeometryDeck(); }

void io::Input::setFractureDeck() { d_fractureDeck_p = new io::FractureDeck(); }

void io::Input::setInitialConditionDeck() {
  d_initialConditionDeck_p = new io::InitialConditionDeck();
}

void io::Input::setInteriorFlagsDeck() {
  d_interiorFlagsDeck_p = new io::InteriorFlagsDeck();
}

void io::Input::setLoadingDeck() { d_loadingDeck_p = new io::LoadingDeck(); }

void io::Input::setMaterialDeck() { d_materialDeck_p = new io::MaterialDeck(); }

void io::Input::setNeighborDeck() { d_neighborDeck_p = new io::NeighborDeck(); }

void io::Input::setOutputDeck() { d_outputDeck_p = new io::OutputDeck(); }

void io::Input::setPolicyDeck() { d_policyDeck_p = new io::PolicyDeck(); }

void io::Input::setProblemDeck() { d_problemDeck_p = new io::ProblemDeck(); }

void io::Input::setSolverDeck() { d_solverDeck_p = new io::SolverDeck(); }

const std::string io::Input::getSpatialDiscretization() {
  return d_problemDeck_p->d_spatialDiscretization;
}
