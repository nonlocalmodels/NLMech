// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "quadrature.h"
#include "../inp/decks/quadratureDeck.h"

fe::Quadrature::Quadrature(inp::QuadratureDeck *deck) {
  d_quadratureDeck_p = deck;

  if (d_quadratureDeck_p->d_quadOrder == 0)
    d_quadratureDeck_p->d_quadOrder = 1;
  if (d_quadratureDeck_p->d_quadOrderM == 0)
    d_quadratureDeck_p->d_quadOrderM = d_quadratureDeck_p->d_quadOrder;
};