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

  /*! @brief Factor to increase dt out */
  size_t d_factorDtOut;

  /*! @brief Initial size of time steps (or frequency) for output operation */
  size_t d_dtOutInit;

  /*!
   * @brief Criteria to be used in deciding when to increase the output
   * frequency (more details will be in the model where criteria based
   * output frequency is implemented, e.g. model::FDModel)
   *
   * Valid values are
   * - ""
   * - "Z_max"
   */
  std::string d_criteriaOut;

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

  /*! @brief Compressor type for .vtu files */
  std::string d_compressType;

  /** @}*/

  /*!
   * @brief Constructor
   */
  OutputDeck()
      : d_dtOut(0), d_dtOutInit(0), d_factorDtOut(1), d_debug(0),
        d_performFEOut(true){};

  /*!
   * @brief Searches list of tags and returns true if the asked tag is in the
   * list
   * @param tag Tag to search
   * @return True/false If tag is find return true or else false
   */
  bool isTagInOutput(const std::string &tag) {

    // search for tag in output tag list
    for (const auto &type : d_outTags)
      if (tag == type)
        return true;

    return false;
  };
};

/** @}*/

} // namespace inp

#endif // INP_OUTPUTDECK_H
