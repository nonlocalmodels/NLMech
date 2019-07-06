////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef MATERIAL_PD_INFLUENCEFN_H
#define MATERIAL_PD_INFLUENCEFN_H

#include <cstring>
#include <vector>

namespace material {

namespace pd {

/*! @brief A base class for computing influence function */
class BaseInfluenceFn {

public:
  /*! @brief Constructor */
  BaseInfluenceFn() = default;

  /*!
   * @brief Returns the value of influence function
   *
   * @param r Reference (initial) bond length
   * @return value Influence function at r
   */
  virtual double getInfFn(const double &r) = 0;

  /*!
   * @brief Returns the moment of influence function
   *
   * If \f$ J(r) \f$ is the influence function for \f$ r\in [0,1)\f$ then \f$
   * i^{th} \f$ moment is given by \f[ M_i = \int_0^1 J(r) r^i dr. \f]
   *
   * @param i ith moment
   * @return moment Moment
   */
  virtual double getMoment(const size_t &i) = 0;
};

/*! @brief A class to implement constant influence function */
class ConstInfluenceFn : public BaseInfluenceFn {

public:
  /*!
   * @brief Constructor
   *  @param params List of parameters
   *  @param dim Dimension
   */
  ConstInfluenceFn(const std::vector<double> &params, const size_t &dim);

  /*!
   * @brief Returns the value of influence function
   *
   * @param r Reference (initial) bond length
   * @return value Influence function at r
   */
  double getInfFn(const double &r) override;

  /*!
   * @brief Returns the moment of influence function
   *
   * If \f$ J(r) \f$ is the influence function for \f$ r\in [0,1)\f$ then \f$
   * i^{th}\f$ moment is given by \f[ M_i = \int_0^1 J(r) r^i dr. \f]
   *
   * @param i ith moment
   * @return moment Moment
   */
  double getMoment(const size_t &i) override;

private:
  /*! @brief Constant such that J(r) = Constant */
  double d_a0;
};

/*! @brief A class to implement linear influence function
 *
 * \f$ J(r) = a0 + a1 r \f$
 */
class LinearInfluenceFn : public BaseInfluenceFn {

public:
  /*!
   * @brief Constructor
   *  @param params List of parameters
   *  @param dim Dimension
   */
  LinearInfluenceFn(const std::vector<double> &params, const size_t &dim);

  /*!
   * @brief Returns the value of influence function
   *
   * @param r Reference (initial) bond length
   * @return value Influence function at r
   */
  double getInfFn(const double &r) override;

  /*!
   * @brief Returns the moment of influence function
   *
   * If \f$ J(r) \f$ is the influence function for \f$ r\in [0,1)\f$ then \f$
   * i^{th}\f$ moment is given by \f[ M_i = \int_0^1 J(r) r^i dr. \f]
   *
   * @param i ith moment
   * @return moment Moment
   */
  double getMoment(const size_t &i) override;

private:
  /*! @brief Constants such that J(r) = d_a0 + d_a1 * r */
  double d_a0;

  /*! @brief Constants such that J(r) = d_a0 + d_a1 * r */
  double d_a1;
};

/*! @brief A class to implement Gaussian influence function
 *
 * \f$ J(r) = \alpha \exp(-r^2/\beta) \f$
 */
class GaussianInfluenceFn : public BaseInfluenceFn {

public:
  /*!
   * @brief Constructor
   *  @param params List of parameters
   *  @param dim Dimension
   */
  GaussianInfluenceFn(const std::vector<double> &params, const size_t &dim);

  /*!
   * @brief Returns the value of influence function
   *
   * @param r Reference (initial) bond length
   * @return value Influence function at r
   */
  double getInfFn(const double &r) override;

  /*!
   * @brief Returns the moment of influence function
   *
   * If \f$ J(r) \f$ is the influence function for \f$ r\in [0,1)\f$ then \f$
   * i^{th}\f$ moment is given by \f[ M_i = \int_0^1 J(r) r^i dr. \f]
   *
   * @param i ith moment
   * @return moment Moment
   */
  double getMoment(const size_t &i) override;

private:
  /*! @brief Constants */
  double d_alpha;

  /*! @brief Constants */
  double d_beta;
};

} // namespace pd

} // namespace material

#endif // MATERIAL_PD_INFLUENCEFN_H
