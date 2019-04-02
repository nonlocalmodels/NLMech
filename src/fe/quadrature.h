// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef FE_QUADRATURE_H
#define FE_QUADRATURE_H

// forward declaration
namespace inp {
struct QuadratureDeck;
}

namespace fe {
class Quadrature {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   */
  explicit Quadrature(inp::QuadratureDeck *deck);

private:
  /**
   * \defgroup Data which do not belong to mesh directly
   */
  /**@{*/

  /*! @brief Quadrature deck */
  inp::QuadratureDeck *d_quadratureDeck_p;

  /** @}*/
};

} // namespace fe

#endif // FE_QUADRATURE_H
