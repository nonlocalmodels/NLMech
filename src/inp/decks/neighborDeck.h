// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef NEIGHBORDECK_H
#define NEIGHBORDECK_H

namespace inp {

/*! @brief Structure to read and store neighbor list related input data */
struct NeighborDeck {

  /**
   * \defgroup Data members
   */
  /**@{*/

  /*! @brief Safety factor for neighbor list calculation */
  double d_safetyFactor;

  /*! @brief Flag which specifies if the partially inside nodes (elements if
   * the discretization is weak finite element) should be included in the
   * neighbor list
   */
  bool d_addPartialElems;

  /** @}*/

  /*!
   * @brief Constructor
   */
  NeighborDeck() : d_safetyFactor(1.0), d_addPartialElems(false){};
};

} // namespace inp
#endif // NEIGHBORDECK_H
