// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef GEOM_INTERIORFLAGS_H
#define GEOM_INTERIORFLAGS_H

#include <vector>

// forward declaration of interior flags deck
namespace inp {
struct InteriorFlagsDeck;
}

namespace geometry {

/*! @brief A class to store interior/exterior flags of node
 *
 * In this class we store the the flag which indicates if the node is inside
 * the material domain or it is near the boundary. This is useful when we
 * implement \a no-fail \a region in \b Peridynamics.
 */
class InteriorFlags {

public:
  /*!
   * @brief Constructor
   * @param deck Input deck which contains user-specified information
   */
  explicit InteriorFlags(inp::InteriorFlagsDeck *deck);

private:
  /*! @brief Interior flags deck */
  inp::InteriorFlagsDeck *d_interiorFlagsDeck_p;

  /*! @brief Vector of flags
   *
   * Since char data provides 4 bits, we club flags of 4 nodes as one store the
   * flags in char data.
   *
   * Given node \a i, to find the interior flag, we proceed as follows:
   *
   * - Location in vector d_intFlags: \a j = \a i/4
   * - Bit location: \a b = \a i%4
   *
   * Now we check bit at location \a b of d_intFlags[ \a j ]. If that bit is
   * \a 0, then the node \a i is in the interior, otherwise it is in the
   * exterior.
   */
  std::vector<char> d_intFlags;

  /** @}*/
};

} // namespace geometry

#endif // GEOM_INTERIORFLAGS_H
