// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef INP_NEIGHBORDECK_H
#define INP_NEIGHBORDECK_H

namespace inp {

/**
 * \ingroup Input
 */
/**@{*/

/*! @brief Structure to read and store neighbor list related input data */
struct NeighborDeck {

  /**
   * @name Data members
   */
  /**@{*/

  /*! @brief Safety factor for neighbor list calculation */
  double d_safetyFactor;

  /*!
   * @brief Flag to include partially inside nodes in neighbor list
   *
   *
   * Flag which specifies if the partially inside nodes (elements if
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

/** @}*/

} // namespace inp

#endif // INP_NEIGHBORDECK_H
