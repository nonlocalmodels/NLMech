// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "material.h"
#include "../inp/decks/materialDeck.h"

material::Material::Material(inp::MaterialDeck *deck){}

bool material::Material::isStateActive() { return false; }
