// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef INP_OUTPUTDECK_H
#define INP_OUTPUTDECK_H

#include <string>
#include <vector>

namespace inp {

/**
 * \ingroup Input
 */
/**@{*/

/*! @brief Structure to read and store output data */
struct OutputDeck {

  /**
   * @name Data members
   */
  /**@{*/

  /*! @brief Output format: currently supports vtk output */
  std::string d_outFormat;

  /*! @brief Output Path where the files will be written */
  std::string d_path;

  /*! @brief List of tags of data to be dumped */
  std::vector<std::string> d_outTags;

  /*! @brief Size of time steps (or frequency) for output operation */
  size_t d_dtOut;

  /*! @brief Flag specifying debug level */
  size_t d_debug;

  /*!
   * @brief Flag specifying if element-node connectivity should not be
   * dumped
   *
   * for large mesh, writing element-node connectivity using vtk writer
   * creates problem.
   */
  bool d_performFEOut;

  /** @}*/

  /*!
   * @brief Constructor
   */
  OutputDeck() : d_dtOut(0), d_debug(0), d_performFEOut(true){};
};

/** @}*/

} // namespace inp

#endif // INP_OUTPUTDECK_H
