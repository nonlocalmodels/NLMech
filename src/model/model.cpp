// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "model.h"

size_t model::Model::currentStep() { return d_n; }

float model::Model::getEnergy() { return d_te - d_tw + d_tk; }