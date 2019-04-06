// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef INP_FRACTUREDECK_H
#define INP_FRACTUREDECK_H

#include "../../util/point.h"           // definition of Point3
#include "../../util/compare.h"         // compare utility
#include "../../util/transfomation.h"   // rotation of point
#include <cstdlib>                      // definition of siz_t
#include <vector>

namespace inp {

/*! @brief A structure to edge crack of any orientation */
struct EdgeCrack {

  /**
   * @name Data members
   */
  /**@{*/

  /*!
   * @brief Orientation of crack
   *
   * - 1 for orientation along horizontal axis (x-axis),
   * - -1 for orientation along vertical axis (y-axis),
   * - 0 for any crack which makes nonzero angle with the horizontal axis
   */
  int d_o;

  /*!
   * @brief Angle the edge crack makes with the horizontal axis
   *
   * - 0 if orientation = 1,
   * - \f$ \pi/2 \f$ if orientation = -1
   * - Value for orientation = 0
   */
  double d_theta;

  /*! @brief Total current length of crack */
  double d_l;

  /*! @brief Current length of crack on top side (right side) */
  double d_lt;

  /*! @brief Current length of crack on bottom side (left side) */
  double d_lb;

  /*! @brief Velocity of top (right) crack tip */
  std::vector<double> d_vt;

  /*! @brief Velocity of bottom (left) crack tip */
  std::vector<double> d_vb;

  /*! @brief Closest id of node to top (right) crack tip */
  int d_it;

  /*! @brief Closest id of node to bottom (left) crack tip */
  int d_ib;

  /*! @brief Closest id of node to old top (right) crack tip */
  int d_iot;

  /*! @brief Closest id of node to old bottom (left) crack tip */
  int d_iob;

  /*! @brief Top (right) crack tip location */
  util::Point3 d_pt;

  /*! @brief Bottom (left) crack tip location */
  util::Point3 d_pb;

  /*! @brief Old top (right) crack tip location */
  util::Point3 d_pot;

  /*! @brief Old bottom (left) crack tip location */
  util::Point3 d_pob;

  /** @}*/

  /*!
   * @brief Constructor
   */
  EdgeCrack()
      : d_o(1), d_theta(0.), d_l(0.), d_lt(0.), d_lb(0.), d_it(-1), d_ib(-1),
        d_iot(-1), d_iob(-1), d_vt(std::vector<double>(3, 0.)),
        d_vb(std::vector<double>(3, 0.)), d_pt(util::Point3()),
        d_pb(util::Point3()), d_pot(util::Point3()), d_pob(util::Point3()){};

  /*!
   * @brief Checks if point lies outside the crack
   * @param p Point p to check
   * @param o Orientation of crack
   * @param pb Bottom (left) point of crack line
   * @param pt Top (right) point of crack line
   * @param theta Angle that crack line makes with horizontal axis
   * @return true if lies outside
   */
  bool ptOutside(util::Point3 p, int o, util::Point3 pb, util::Point3 pt,
                 double theta = 0.0) {

    if (o == -1) {

      // straight crack along y-axis
      return util::compare::definitelyLessThan(p.d_y, pb.d_y) or
             util::compare::definitelyGreaterThan(p.d_y, pt.d_y);

    } else if (o == 1) {

      // straight crack along x-axis ===> pb == pl, pt == pr
      return util::compare::definitelyLessThan(p.d_x, pb.d_x) or
             util::compare::definitelyGreaterThan(p.d_x, pt.d_x);

    } else if (o == 0) {

      // straight crack at theta angle with x-axis ===> pb == pl, pt == pr

      // 1. pb ==> (0,0), pt ==> pt - pb, p ==> p - pb
      // 2. Apply CW rotation to new pt and new p so that crack line
      // after rotation is simply along x-axis with left point at
      // origin and right point at transformed new pt

      util::Point3 pmap = util::transformation::rotateCW2D(p - pb, theta);
      util::Point3 ptmap = util::transformation::rotateCW2D(pt - pb, theta);

      return util::compare::definitelyLessThan(pmap.d_x, 0.0) or
             util::compare::definitelyGreaterThan(pmap.d_x, ptmap.d_x);
    }

    return true;
  }

  /*!
   * @brief Checks if point lies on left(top) or right(bottom) of crack
   * @param p Point p to check
   * @param pb Bottom (left) point of crack line
   * @param pt Top (right) point of crack line
   * @return true if lies on left(top)
   */
  bool ptLeftside(util::Point3 p, util::Point3 pb, util::Point3 pt) {

    // algorithm is same for any orientation of crack
    //            pt
    //           o
    //          /
    //  p      /
    //   o    /
    //       /
    //      /
    //     /
    //    o
    //    pb
    double a = 0.5 *( (pt.d_x - pb.d_x)*(p.d_y - pb.d_y)
                    - (p.d_x - pb.d_x)*(pt.d_y - pb.d_y));

    // crack is closer to left nodes
    return !util::compare::definitelyLessThan(a, 0.0);
  }

  /*!
   * @brief Checks if point lies on left(top) or right(bottom) of crack
   * @param p Point p to check
   * @param pb Bottom (left) point of crack line
   * @param pt Top (right) point of crack line
   * @return true if lies on right(bottom)
   */
  bool ptRightside(util::Point3 p, util::Point3 pb, util::Point3 pt) {

    return !this->ptLeftside(p, pb, pt);
  }
};

/**
 * \ingroup Input
 */
/**@{*/

/*! @brief Structure to read and store fracture related input data */
struct FractureDeck {

  /**
   * @name Data members
   */
  /**@{*/

  /*! @brief Vector of pre-crack data*/
  std::vector<inp::EdgeCrack> d_cracks;

  /*! @brief Flag which indicates if pre-crack is present */
  bool d_crackActive;

  /*!
   * @brief Output interval
   *
   * I.e. code perform crack data output every N number of steps, when N is
   * given by this data
   */
  size_t d_crackOutput;

  /** @}*/

  /*!
   * @brief Constructor
   */
  FractureDeck() : d_crackActive(false), d_crackOutput(0){};
};

/** @}*/

} // namespace inp

#endif // INP_FRACTUREDECK_H
