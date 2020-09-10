////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef INP_POLICY_H
#define INP_POLICY_H

#include <string>
#include <vector>

// forward declaration of policy deck
namespace inp {
struct PolicyDeck;
}

namespace inp {

/*! @brief A class to enforce certain policies to reduce memory loads
 *
 * We implement simple method to control population of data in simulation.
 * Implementation can be made specific to particular model by simply
 * assigning a tag to d_modelTag and defining a new rule for the tag in the
 * inp::Policy::init.
 *
 * For a given memory control flag i (can be 0,1,2,3), we look at the list of
 * tags in inp::Policy::d_lTags to know whether we populate the data (given
 * by tag, e.g. tag for data d_u in model::Model is Model_g_u) in the
 * simulation. If it is in the list inp::Policy::d_lTags[i] then we do not
 * populate this data in the simulation.
 *
 * inp::Policy::d_lTags data is created in init() which is called when the
 * getInstance() is invoked for the first time.
 *
 * @note If inp::Policy::d_enablePostProcessing is set to false then no
 * postprocessing calculation will be carried out and no postprocessing data
 * will be populated. So in a way inp::Policy::d_enablePostProcessing = false
 * acts as  strictest control.
 *
 * @todo Modify constructor and getInstance() to set d_modelTag.
 */
class Policy {

public:

  /*!
   * @brief Returns the pointer to static class. Creates instance in its
   * first call
   * @param deck Input deck which contains user-specified information
   * @return Policy instance of static class Policy
   */
  static Policy *getInstance(inp::PolicyDeck *deck = nullptr);

  /*!
   * @brief Destroys the instance
   */
  static void destroyInstance();

  /**
   * @name Setter method
   */
  /**@{*/

  /*!
   * @brief Adds tag to specified level tag list
   * @param level Level of tag list
   * @param tag Tag to be appended to list
   */
  void addToTags(const size_t &level, const std::string &tag);

  /*!
   * @brief Looks for tag in the level d_memControlFlag and if present
   * removes it
   * @param tag Tag to be appended to list
   */
  void removeTag(const std::string &tag);

  /** @}*/

  /**
   * @name Getter method
   */
  /**@{*/

  /*!
   * @brief Returns true/false depending on whether tag is found
   * @param tag Tag to search for
   * @return bool True if it can be populated, false otherwise
   */
  bool populateData(const std::string &tag);

  /*!
   * @brief Returns memory control flag
   * @return Flag
   */
  int getMemoryControlFlag();

  /*!
   * @brief Returns true if post-processing computation is to be done
   * @return Flag
   */
  bool enablePostProcessing();

  /** @}*/

private:
  /*!
   * @brief Private constructor
   * @param deck Input deck which contains user-specified information
   */
  explicit Policy(inp::PolicyDeck *deck = nullptr);

  /*! @brief Private operator */
  Policy(Policy const &);

  /*!
   * @brief Private copy operator
   * @return Policy object after copying
   */
  const Policy &operator=(const Policy &);

  /*! @brief Private destructor */
  ~Policy();

  /*! @brief Initializes the data */
  void init();

  /*! @brief Static instance of Policy class */
  static Policy *d_instance_p;

  /*!
   * @brief Flag which indicates level of memory control to be enforced
   *
   * Default is 0 which means no control. Max at present is 2 which means as
   * much control as possible.
   */
  int d_memControlFlag;

  /*!
   * @brief Enable post-processing calculation
   *
   * Default is true.
   */
  bool d_enablePostProcessing;

  /*! @brief Specify model tag */
  std::string d_modelTag;

  /*! @brief Specify maximum level of memory control */
  size_t d_maxLevel;

  /*! @brief List of variable names in different levels to help enforce the
   * memory control
   */
  std::vector<std::vector<std::string>> d_lTags;
};

} // namespace inp

#endif // IO_POLICY_H
