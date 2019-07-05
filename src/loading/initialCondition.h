// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef LOADING_INITIALCONDITION_H
#define LOADING_INITIALCONDITION_H

#include "util/point.h"          // definition of Point3
#include <string>
#include <vector>

// forward declaration of initial condition deck
namespace inp {
struct InitialConditionDeck;
}

// forward declaration
namespace fe {
class Mesh;
}

namespace loading {

/*!
 * @brief A class to apply initial condition
 *
 * This class processes input data and provides method to apply initial
 * condition.
 */
class InitialCondition {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   */
  explicit InitialCondition(inp::InitialConditionDeck *deck);

  /*!
   * @brief Applies initial condition to displacement and velocity
   * @param u Vector nodal displacements
   * @param v Vector nodal velocities
   * @param mesh Mesh object
   */
  void apply(std::vector<util::Point3> *u, std::vector<util::Point3> *v,
             fe::Mesh *mesh);

private:
  /*!
   * @brief Computes the formula specified by input file
   * @param fn_type Type of function in formula
   * @param params List of required parameters
   * @param x Coordinate of point
   * @param dof Degree of freedom (0 for X, 1 for Y, 2 for Z)
   * @param dim Dimension
   * @return value Value of formula at the point
   */
  double getICFormula(const std::string &fn_type,
                      const std::vector<double> &params, const util::Point3 &x,
                      const size_t &dof, const size_t &dim);

  /*! @brief Initial condition deck */
  inp::InitialConditionDeck *d_deck_p;
};

} // namespace loading

#endif // LOADING_INITIALCONDITION_H
