// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef GEOM_FRACTURE_H
#define GEOM_FRACTURE_H

#include "util/point.h" // definition of Point3
#include <inp/decks/fractureDeck.h>
#include <stdint.h> // uint8_t type
#include <string.h> // size_t type
#include <vector>

// forward declaration of fracture deck
namespace inp {
struct EdgeCrack;
struct FractureDeck;
} // namespace inp

/*!
 * @brief Collection of methods and database related to geometry
 *
 * This namespace provides methods and data members specific to geometry. It
 * consists of class Fracture, InteriorFlags, and Neighbor.
 *
 * @sa Fracture, InteriorFlags, Neighbor
 */
namespace geometry {

/*! @brief A struct to compute crack tip and crack tip velocity */
struct CrackOutData {

  /*! @brief Old update time for top (right) side of crack */
  double d_timet;

  /*! @brief Old update time for bottom (left) side of crack */
  double d_timeb;

  /*! @brief Total number of updates at current time */
  size_t d_updateCount;

  /*! @brief Total number of files created so far */
  size_t d_fileOutCount;

  /*! @brief Flag which specifies if we need to create new file */
  bool d_needNewFile;

  /*! @brief File to which crack data will be written */
  FILE *d_file;

  /*! @brief Default constructor */
  CrackOutData()
      : d_timet(0.), d_timeb(0.), d_updateCount(0), d_fileOutCount(0),
        d_needNewFile(false), d_file(nullptr){};
};

/*! @brief A class for fracture state of bonds
 *
 * In this class fracture state of each bonds (i.e. whether the bond is
 * broken or not) is stored. It also comes with the access to the state
 * of bond.
 */
class Fracture {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   * @param nodes Pointer to nodal coordinates
   * @param neighbor_list Pointer to neighbor list
   */
  Fracture(inp::FractureDeck *deck, const std::vector<util::Point3> *nodes,
           const std::vector<std::vector<size_t>> *neighbor_list);

  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   */
  explicit Fracture(inp::FractureDeck *deck);

  /*!
   * @brief Sets the bond state
   *
   * @param i Nodal id
   * @param j Local id of bond in neighbor list of i
   * @param state State which is applied to the bond
   */
  void setBondState(const size_t &i, const size_t &j, const bool &state);

  /*!
   * @brief Sets the output time step to given value
   *
   * @param dt Output time step interval
   */
  void setUpdateCrack(const size_t &dt);

  /*!
   * @brief Sets the bond state
   *
   * @param i Nodal id
   * @param j Local id of bond in neighbor list of i
   * @return state True if bond is fractured otherwise false
   */
  bool getBondState(const size_t &i, const size_t &j);

  /*!
   * @brief Returns the bonds of given node i
   *
   * @param i Nodal id
   * @return List Bonds of node i
   */
  const std::vector<uint8_t> getBonds(const size_t &i);

  /*!
   * @brief Updates crack tip and velocity and performs output
   *
   * @param n Time step
   * @param time Current time
   * @param output_path Path where crack data file will be created
   * @param horizon Horizon
   * @param nodes Reference coordinates of all nodes
   * @param u Vector of displacement
   * @param Z Vector of damage at nodes
   */
  void
  updateCrackAndOutput(const size_t &n, const double &time,
                       const std::string &output_path,
                       const double &horizon,
                       const std::vector<util::Point3> *nodes,
                       std::vector<util::Point3> *u, std::vector<float> *Z);

  /*!
   * @brief Updates crack tip and velocity and performs output
   * @return dt Time step interval for crack output
   */
  size_t getDtCrackOut();

private:
  /*!
   * @brief Sets flag of bonds of i as fractured which intersect the
   * pre-crack
   *
   * @param i Nodal id
   * @param crack Pointer to the pre-crack
   * @param nodes Pointer to nodal coordinates
   * @param neighbors Pointer to neighbors of node i
   */
  void computeFracturedBondFd(const size_t &i, inp::EdgeCrack *crack,
                              const std::vector<util::Point3> *nodes,
                              const std::vector<size_t> *neighbors);

  /*!
   * @brief Updates crack tip location and crack velocity
   * @param n Time step
   * @param time Current time
   * @param horizon Horizon
   * @param nodes Reference coordinates of all nodes
   * @param u Vector of displacement
   * @param Z Vector of damage at nodes
   */
  void updateCrack(const size_t &n, const double &time, const double &horizon,
                   const std::vector<util::Point3> *nodes,
                   std::vector<util::Point3> *u, std::vector<float> *Z);

  /*!
   * @brief Performs output
   * @param n Time step
   * @param time Current time
   * @param output_path Path where crack data file will be created
   * @param nodes Reference coordinates of all nodes
   * @param u Vector of displacement
   * @param Z Vector of damage at nodes
   */
  void output(const size_t &n, const double &time,
              const std::string &output_path,
              const std::vector<util::Point3> *nodes,
              std::vector<util::Point3> *u, std::vector<float> *Z);

  /*! @brief Interior flags deck */
  inp::FractureDeck *d_fractureDeck_p;

  /*! @brief Data for crack tip calculation */
  geometry::CrackOutData d_crackOutData;

  /*! @brief Vector which stores the state of bonds
   *
   * Given node i, vector d_fracture[i] is the list of state of bonds of node
   * i.
   *
   * This is the most memory efficient data where 1 bit represents the state
   * of bond.
   */
  std::vector<std::vector<uint8_t>> d_fracture;
};

} // namespace geometry

#endif // GEOM_FRACTURE_H
