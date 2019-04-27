// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

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
 * This class introduces policies which restrict population (declaration) of
 * post-processing data which are not important for running simulation and
 * are postprocessing data.
 *
 * For example, if the simulation is large, this class will restrict
 * population of data such as fracture energy, damage data, strain and stress
 * data, etc.
 *
 * This class also enforces lumping approximation of mass matrix if the level
 * of restriction is set to 2 or higher.
 *
 * The level of memory restriction can be set in the input file. Default
 * value is 0 which means no restriction.
 */
class Policy {

public:
  /**
   * @name Get and destroy instance
   */
  /**@{*/

  /*!
   * @brief Returns the pointer to static class. Creates instance in its
   * first call
   */
  static Policy *getInstance();

  /*!
   * @brief Returns the pointer to static class. Creates instance in its
   * first call
   * @param deck Input deck which contains user-specified information
   * @return Policy instance of static class Policy
   */
  static Policy *getInstance(inp::PolicyDeck *deck);

  /*!
   * @brief Destroys the instance
   */
  static void destroyInstance();

  /** @}*/

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

  /** @}*/

  /**
   * @name Getter method
   */
  /**@{*/

  /*!
   * @brief Returns true/false depending on whether tag is found
   * @param tag Tag to search for
   * @return True/False
   */
  bool populateData(const std::string &tag);

  /*!
   * @brief Returns memory control flag
   * @return Flag
   */
  int getMemoryControlFlag();

  /*!
   * @brief Returns true if postprocessing computation is to be done
   * @return Flag
   */
  bool enablePostProcessing();

  /** @}*/

private:
  /*!
   * @brief Private constructor
   */
  Policy();

  /*!
   * @brief Private constructor
   * @param deck Input deck which contains user-specified information
   */
  explicit Policy(inp::PolicyDeck *deck);

  /*! @brief Private operator */
  Policy(Policy const &);

  /*!
   * @brief Private copy operator
   * @return Policy object after copying
   */
  const Policy &operator=(const Policy &);

  /*! @brief Private destructor */
  ~Policy();

  /**
   * @name Private methods
   */
  /**@{*/

  /*! @brief Initializes the data */
  void init();

  /** @}*/

  /**
   * @name Internal data
   */
  /**@{*/

  /*! @brief Static instance of Policy class */
  static Policy *d_instance_p;

  /*! @brief Policy deck which contains input data */
  inp::PolicyDeck *d_policyDeck_p;

  /*! @brief List of variable names in different levels to help enforce the
   * memory control
   */
  std::vector<std::string> d_l0Tags;
  std::vector<std::string> d_l1Tags;
  std::vector<std::string> d_l2Tags;
  std::vector<std::string> d_otherTags;

  /** @}*/
};

} // namespace inp

#endif // IO_POLICY_H
