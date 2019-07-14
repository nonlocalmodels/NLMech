// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef TOOLS_PP_UTIL_H
#define TOOLS_PP_UTIL_H

#include "inp/decks/fractureDeck.h" // definition of EdgeCrack
#include "util/point.h"             // definition of Point3
#include <vector>
#include <string>

namespace tools {

namespace pp {

/*!
 * @brief Datatype used in sorting rectangles for crack tip search
 */
struct SortZ {

  /*! @brief Id of node */
  size_t d_i;

  /*! @brief Id of rectangle */
  size_t d_r;

  /*! @brief Damage associated to the node */
  double d_Z;

  SortZ() : d_i(0), d_r(0), d_Z(0.){};
};

/*!
 * @brief Datatype to hold crack tip data
 */
struct CrackTipData {

  /*! @brief Output time step */
  int d_n;

  /*! @brief Crack tip location */
  util::Point3 d_p;

  /*! @brief Crack tip velocity */
  util::Point3 d_v;

  /*!
   * @brief Constructor
   */
  CrackTipData() : d_n(0), d_p(util::Point3()), d_v(util::Point3()){};

  /*!
   * @brief Constructor
   *
   * @param n Output time step
   * @param p Crack tip location
   * @param v Crack tip velocity
   */
  CrackTipData(size_t n, util::Point3 p, util::Point3 v)
      : d_n(n), d_p(p), d_v(v){};
};

/*!
 * @brief Data for transforming displacement operation
 */
struct TransformU {
  /*! @brief Factor for scaling the displacement */
  double d_scale;

  /*!
   * @brief Constructor
   */
  TransformU() : d_scale(1.) {};
};

/*!
 * @brief Data for marking and symmetrizing of velocity operation
 */
struct TransformVelocity {
  /*! @brief Specify if velocity is to be marked */
  bool d_markVAsZero;

  /*! @brief True if rectangle region for marking is specified */
  bool d_markVInRectGiven;

  /*! @brief Rectangle region for marking */
  std::pair<util::Point3, util::Point3> d_markVRect;

  /*! @brief List of points at which velocities should be marked */
  std::vector<util::Point3> d_markVPts;

  /*! @brief Are points provided in current configuration */
  bool d_markVPtsAreInCurrentConfig;

  /*! @brief List of nodes at which velocities should be marked */
  std::vector<size_t> d_markVNodes;

  /*! @brief Specify if velocity is to be symmetrized */
  bool d_symmetrizeV;

  /*! @brief Specify if combine mark and symmetrize operation */
  bool d_combineMarkV;

  /*! @brief Axis of symmetry, e.g. "x", "y" */
  std::string d_symmAxis;

  /*! @brief Specify the axis location */
  double d_symmLine;

  /*!
   * @brief Constructor
   */
  TransformVelocity()
      : d_markVAsZero(false), d_markVInRectGiven(false),
        d_markVPtsAreInCurrentConfig(false), d_symmetrizeV(false),
        d_combineMarkV(false), d_symmLine(0.) {};
};

/*!
 * @brief Data for strain and its magnitude computation operation
 */
struct ComputeStrain {
  /*! @brief Factor for scaling the displacement */
  bool d_computeStrain;

  /*! @brief Compute magnitude of strain tensor */
  bool d_magStrainTensor;

  /*!
   * @brief Specify component (if any) of which absolute value is to be
   * computed
   */
  std::string d_magStrainComp;

  /*!
   * @brief Specify element ids at which magnitude of strain should be
   * marked
   */
  std::vector<std::pair<size_t, double>> d_markMagStrainCells;

  /*!
   * @brief Constructor
   */
  ComputeStrain() : d_computeStrain(false), d_magStrainTensor(false) {};
};

/*!
 * @brief Data for crack tip computation
 */
struct FindCrackTip {

  /*! @brief Old update time for top (right) side of crack */
  double d_timet;

  /*! @brief Old update time for bottom (left) side of crack */
  double d_timeb;

  /*! @brief Total number of updates at current time */
  size_t d_updateCount;

  /*! @brief Set of edge cracks */
  std::vector<inp::EdgeCrack> d_cracks;

  /*! @brief Set of edge cracks */
  double d_minZAllowed;

  /*! @brief Set of edge cracks */
  double d_maxZAllowed;

  /*! @brief File (for top tip) to which crack data will be written */
  FILE *d_filet;

  /*! @brief File (for bottom tip) to which crack data will be written */
  FILE *d_fileb;

  /*!
   * @brief Constructor
   */
  FindCrackTip()
      : d_timet(0.), d_timeb(0.), d_updateCount(0),
        d_minZAllowed(1.), d_maxZAllowed(50.), d_filet(nullptr),
        d_fileb(nullptr){};
};

/*!
 * @brief Data for J integral calculation computation
 */
struct ComputeJIntegral {
  /*! @brief Specify orientation of crack. 1 for horizontal, -1 for vertical */
  int d_crackOrient;

  /*! @brief Specify crack id for which data to be read from file */
  int d_crackId;

  /*! @brief Specify file from which crack tip information is to be read */
  std::string d_crackTipFile;

  /*! @brief Factor of horizon for defining contour */
  std::vector<double> d_contourFactor;

  /*! @brief Data to hold crack tip information */
  std::vector<tools::pp::CrackTipData> d_crackTipData;

  /*! @brief File to write J integral data */
  FILE *d_file;

  /*! @brief Set lateral component of velocity of crack tip as zero */
  bool d_setLateralCompVZero;

  /*! @brief Set lateral component of displacement of crack tip as zero */
  bool d_setLateralCompUZero;

  /*!
   * @brief Set lateral component of crack tip location to given value (it
   * should be initial value of lateral component of crack tip)
   */
  double d_setLateralCompX;

  /*!
   * @brief Constructor
   */
  ComputeJIntegral()
      : d_crackOrient(0), d_crackId(1), d_file(nullptr),
        d_setLateralCompVZero(false), d_setLateralCompUZero(false),
        d_setLateralCompX(0.){};
};

/*!
 * @brief Datatype to hold instructions for postprocessing operation
 */
struct InstructionData {

  /*!
   * @brief Tag assigned to this compute instruction (used for creating output
   * filename)
   */
  std::string d_tagFilename;

  /*! @brief Specify format of the output file, e.g. msh, vtu, legacy_vtk
   *
   * Default is vtu format
   */
  std::string d_outFormat;

  /*! @brief Start output step */
  int d_start;

  /*! @brief End output step */
  int d_end;

  /*!
   * @brief Number which indicates how many simulation file we skip after
   * processing 1 simulation file. Default is 1. Valid number >= 1.
   */
  int d_interval;

  /*!
   * @brief Flag for performing only nodal output and not full fem output
   * (for large mesh, vtk produces error when writing element-node
   * connectivity)
   */
  bool d_outOnlyNodes;

  /*!
   * @brief Flag for removing elements which contain damages quadrature
   * point (currently this is not implemented)
   */
  bool d_removeElements;

  /*! @brief Compression type for .vtu files */
  std::string d_compressType;

  /*! @brief Transforming of displacement operation */
  tools::pp::TransformU *d_transformU_p;

  /*! @brief Transform velocity operation */
  tools::pp::TransformVelocity *d_transformV_p;

  /*! @brief Compute strain and its magnitude operation */
  tools::pp::ComputeStrain *d_compStrain_p;

  /*! @brief Find crack tip and compute velocity operation */
  tools::pp::FindCrackTip *d_findCrackTip_p;

  /*! @brief Computing damage at nodes operation */
  bool d_damageAtNodes;

  /*! @brief J-integral calculation operation */
  tools::pp::ComputeJIntegral *d_computeJInt_p;

  InstructionData()
      : d_outFormat("vtu"), d_start(-1), d_end(-1), d_interval(1),
        d_outOnlyNodes(true), d_removeElements(false), d_transformU_p(nullptr),
        d_transformV_p(nullptr), d_compStrain_p(nullptr),
        d_findCrackTip_p(nullptr), d_damageAtNodes(false),
        d_computeJInt_p(nullptr){};
};

} // namespace pp

} // namespace tools

#endif //TOOLS_PP_UTIL_H
