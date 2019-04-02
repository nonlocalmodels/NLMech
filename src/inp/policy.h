// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef IO_POLICY_H
#define IO_POLICY_H

#include <string>
#include <vector>

// forward declaration of policy deck
namespace inp {
struct PolicyDeck;
}

namespace inp {

/*! @brief Methods and database associated to the mesh */
class Policy {

public:
  /**
   * \defgroup Methods to get and destroy instance
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
   */
  static Policy *getInstance(inp::PolicyDeck *deck);

  /*!
   * @brief Destroys the instance
   */
  static void destroyInstance();

  /** @}*/

  /**
   * \defgroup Methods to get and set tags from the list
   */
  /**@{*/

  /*!
   * @brief Adds tag to specified level tag list
   * @param level Level of tag list
   * @param tag Tag to be appended to list
   */
  void addToTags(size_t level, const std::string &tag);

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
  Policy(inp::PolicyDeck *deck);

  /*! @brief Private operator */
  Policy(Policy const &);

  /*! @brief Private operator */
  const Policy &operator=(const Policy &);

  /*! @brief Private destructor */
  ~Policy();

  /**
   * \defgroup Private methods
   */
  /**@{*/

  /*! @brief Initializes the data */
  void init();

  /** @}*/

  /**
   * \defgroup Data members
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
