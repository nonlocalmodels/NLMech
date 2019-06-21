// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef TOOLS_PP_COMPUTE_H
#define TOOLS_PP_COMPUTE_H

#include "fe/mesh.h"                // definition of Mesh
#include "geometry/fracture.h"      // definition of Fracture
#include "geometry/neighbor.h"      // definition of Neighbor
#include "inp/decks/materialDeck.h" // definition of MaterialDeck
#include "inp/decks/modelDeck.h"    // definition of ModelDeck
#include "inp/decks/outputDeck.h"   // definition of OutputDeck
#include "inp/decks/fractureDeck.h"   // definition of FractureDeck
#include "inp/input.h"              // definition of Input
#include "material/pdMaterial.h"    // definition of Material
#include "util.h"
#include "util/matrix.h" // definition of SymMatrix3
#include "rw/writer.h"              // definition of VtkWriterInterface

namespace tools {

/*!
 * @brief Namespace for postprocessing of simulation results
 *
 * In this namespace we define the methods to compute quantities of interest
 * from the simulation results. Examples are: strain and stress, scaling of
 * displacement for fracture plot, symmetrization of plots, etc.
 */
namespace pp {

/*!
 * @brief Processes simulation results and computes postprocessing
 * quantities
 */
class Compute {

public:
  /*!
   * @brief Constructor
   *
   * @param filename Name of input (yaml) file
   */
  explicit Compute(const std::string &filename);

  /*!
   * @brief Initialize the class based on input data
   */
  void init();

  /*!
   * @brief Read compute instruction
   *
   * @param set Tag for given compute set
   * @param data Pointer to instruction data
   */
  void readComputeInstruction(const std::string &set,
                              tools::pp::InstructionData *data);

  /*!
   * @brief Read crack tip data
   */
  void readCrackTipData(const std::string &filename,
                        std::vector<tools::pp::CrackTipData> *data);

  /*!
   * @brief Create new writer object if it is not already done so
   *
   * @param writer Pointer to vtk writer
   */
  void initWriter(rw::writer::VtkWriterInterface *writer,
                  std::vector<util::Point3> *u);

  /**
   * @name Postprocessing methods
   */
  /**@{*/

  /*!
   * @brief Transform displacement
   *
   * Currently following transformation is implemented
   * - Scale displacement by specified factor
   *
   * @param writer Pointer to vtk writer
   */
  void transformU(rw::writer::VtkWriterInterface *writer);

  /*!
   * @brief Transform velocity
   *
   * There are two transformation implemented
   * - Mark velocity of certain nodes as zero
   * - Symmetrize the velocity field
   *
   * @param writer Pointer to vtk writer
   */
  void transformV(rw::writer::VtkWriterInterface *writer);

  /*!
   * @brief Compute strain and stress
   *
   * We perform following postprocessing calculation
   * - compute strain and stress
   * - compute magnitude (or absolute value of specified component of tensor)
   * of strain tensor
   * - mark strain of specified cells to given value
   *
   * 1. Given an element \f$ T \f$ and its shape functions \f$ N_1, N_2,..., N_n
   * \f$, the displacement at quadrature point is given by
   * \f[ u(x_q,y_q) = \sum_{i=1}^n N_i(x_q, y_q) u^i, \f]
   * where \f$ u^i \f$ is displacement of node \f$ i\f$. Linear strain is
   * given by
   * \f[ E(x_q, y_q) = \frac{1}{2}( \nabla u(x_q, y_q) + \nabla u(x_q,y_q)^T ).
   * \f]
   * Components of tensor are given by
   * \f[ E_{xx} = \frac{\partial u_x(x_q, y_q)}{\partial x}, \quad E_{yy} =
   * \frac{\partial u_y(x_q, y_q)}{\partial y} \f]
   * \f[ E_{xy} = E_{yx} = \frac{1}{2}( \frac{\partial u_x(x_q, y_q)
   * }{\partial y} + \frac{\partial u_y(x_q, y_q)}{\partial x}). \f]
   * Therefore we have
   * \f[ E_{xx} = \sum_{i=1}^n\frac{\partial N_i(x_q, y_q)}{\partial x} u^i_x,
   * \quad E_{yy} =  \sum_{i=1}^n \frac{\partial N_i(x_q, y_q)}{\partial y}
   * u^i_y \f]
   * \f[ E_{xy} = E_{yx} = \frac{1}{2}( \sum_{i=1}^n \frac{\partial N_i(x_q, y_q)
   * }{\partial y} u^i_x + \sum_{i=1}^n \frac{\partial N_i(x_q, y_q)
   * }{\partial x} u^i_y). \f]
   *
   * 2. Out of plane component of strain is zero if it is \a plane-stress
   * simulation. If it is \a plane-strain simulation then we use formula
   * \f[
   * E_{zz} = -\frac{\nu}{1 - \nu} (E_{xx} + E_{yy}).\f]
   *
   * 3. We use linear elastic constitutive equation to compute stress
   * \f[ \sigma_{xx} = \lambda (E_{xx} + E_{yy} + E_{zz}) + 2\mu E_{xx},\f]
   * \f[\sigma_{yy} = \lambda (E_{xx} + E_{yy} + E_{zz}) + 2\mu E_{yy},\f]
   * \f[ \sigma_{xy} = 2\mu E_{xy},\quad \sigma_{xz} = 2\mu E_{xz}, \quad
   * \sigma_{yz} = 2\mu E_{yz}.\f]
   *
   * 4. If it is \a plane-strain, we have \f$ \sigma_{zz} = 0 \f$. Otherwise,
   * \f[ \sigma_{zz} = \nu (\sigma_{xx} + \sigma_{yy}) .\f]
   *
   * @param writer Pointer to vtk writer
   */
  void computeStrain(rw::writer::VtkWriterInterface *writer);

  /*!
   * @brief Compute damage at nodes
   *
   * We compute the damage function
   * \f[ Z(x) := \sup_{y\in B_\epsilon(x)} \frac{|S(y,x;u)|}{S_c(y,x)}, \f]
   * where \f$  B_\epsilon(x) \f$ is the ball centered at \f$ x\f$, \f$ S(y,
   * x;u) \f$ is the bond-strain between points \f$y,x\f$, and \f$S_c(y,x)
   * \f$ is the critical strain.
   *
   * @param Z Pointer to nodal damage
   * @param writer Pointer to vtk writer
   * @param perf_out Flag to perform output of damage data
   */
  void computeDamage(std::vector<float> *Z, rw::writer::VtkWriterInterface
  *writer, bool perf_out = false);

  /*!
   * @brief Find crack tip location and velocity
   *
   * To compute the crack tip location at current output time, we first compute
   * the damage at nodes. We then search for node which has smallest damage
   * among nodes with damage above 1.
   *
   * To compute crack tip velocity, we compare the tip at two different times
   * . If not specified, then the two different times are the previour output
   * time and current output time.
   *
   * Since we compute damage in this method, we store the damage in global
   * variable which can then be used in method computeDamage().
   *
   * @param Z Pointer to nodal damage
   * @param writer Pointer to vtk writer
   */
  void findCrackTip(std::vector<float> *Z, rw::writer::VtkWriterInterface
  *writer);

  /*!
   * @brief Compute J integral
   *
   * We consider rectangle near the crack tip and compute J-integral using
   * the contour of rectangle. Schematic for horizontal crack (similar for
   * vertical crack)
   *
   *                         D                    C
   *                         + + + + + + + + + + +
   *                         +                   +
   *                         +                   +
   *       ------------------+                   +
   *                         +                   +
   *                         +                   +
   *                         + + + + + + + + + + +
   *                        A                    B
   *
   * Note: contour is formed by lines A-B, B-C, C-D, D-A
   *
   * Let contour is denoted as \f$ \Gamma(t) \f$, where \f$ t \f$ indicates
   * contour moves with crack tip. Let domain inside contour is defined as
   * \f$ A(t) \f$. Let the outward normal to the domain \f$ A(t) \f$ is \f$ n
   * \f$ and the crack velocity is \f$ v \f$. Then the energy into crack is
   * given by
   * \f[ E(t) = \int_{\Gamma(t)} \left[ \frac{1}{2} \rho |\dot{u}(t)|^2 dt +
   * \bar{W}(x; u(t)) \right] v\cdot n dx - \frac{2}{\epsilon |B_\epsilon(0)|}
   * \int_{A^c(t)} \int_{A(t) \cap B_\epsilon(x)} \partial_S W(S(y,x;
   * u(t))) \frac{y-x}{|y-x|} \cdot (\dot{u}(x,t) + \dot{u}(y,t)) dy dx. \f]
   * Here \f$ \bar{W}(x;u(t)) \f$ is the energy density at point \f$ x\f$ given
   * by \f[ \bar{W}(x;u(t)) = \frac{1}{|B_\epsilon(0)|} \int_{B_\epsilon(x)}
   * |y-x| W(S(y,x;u(t))) dy.\f] \f$ W(S(y,x;u)) \f$ is the pairwise energy
   * density. For regularized bond based model (see material::pd::RNPBond), it
   * is given by
   * \f[ W(S(y,x;u)) = J^\epsilon(|y-x|) \frac{1}{\epsilon |y-x|}
   * \psi(|y-x| S(y,x;u)^2). \f]
   * From above we have
   * \f[\partial_S W(S(y,x;u)) = \frac{2J^\epsilon(|y-x|) S(y,x;u)}{\epsilon}
   * \psi'(|y-x| S(y,x;u)^2). \f]
   * See material::pd::RNPBond for complete details about the material model.
   *
   * Method:
   *
   * 1. We store the ids of nodes which are within horizon distance from the
   * contour. We also store the ids of elements which intersect the contour.
   * All operations below are performed using this list of nodes and elements
   * . This significantly reduces the computation as we no longer loop over
   * whole list of nodes and elements.
   *
   * 2. To compute \f[ E(t) = \int_{\Gamma(t)} \left[ \frac{1}{2} \rho
   * |\dot{u}(t)|^2 dt +
   * \bar{W}(x; u(t)) \right] v\cdot n dx \f]
   * we discretize each edge in contour using 1-d line element. We use second
   * order quadrature rule. We first find the displacement and velocity of
   * the quadrature point using interpolation, see interpolateUV.
   *
   * 3. At quadrature point, we compute peridynamic energy density using
   * pdEnergy.
   *
   * 4. To compute \f[ \frac{2}{\epsilon |B_\epsilon(0)|}
   * \int_{A^c(t)} \int_{A(t) \cap B_\epsilon(x)} \partial_S W(S(y,x;
   * u(t))) \frac{y-x}{|y-x|} \cdot (\dot{u}(x,t) + \dot{u}(y,t)) dy dx \f]
   * we loop over nodes which are in complement of domain \f$A(t) \f$, i.e.
   * \f$ A^c(t) \f$. We compute the contribution of each node in \f$ A^c(t)
   * \f$ using pdForceWork.
   *
   * @param writer Pointer to vtk writer
   */
  void computeJIntegral(rw::writer::VtkWriterInterface *writer);

  /** @}*/

private:
  /**
   * @name Utility methods
   */
  /**@{*/

  /*!
 * @brief Inserts element to list if not present
 *
 * @param i Element to be inserted
 * @param list Pointer to list
 */
  void addUniqueToList(size_t i, std::vector<size_t> *list);

  /*!
   * @brief Find node closest to given point
   *
   * If u is not null then given point is assumed to be in current
   * configuration and therefore closest node in current configuration is
   * searched.
   *
   * @param x Point in reference/current configuration
   * @param nodes Pointer to list of nodes
   * @param u Pointer to nodal displacements
   * @return i Closest node to point x
   */
  size_t findNode(const util::Point3 &x, const std::vector<util::Point3> *nodes,
                  const std::vector<util::Point3> *u = nullptr);

  /*!
   * @brief Find node within contour and elements intersecting contour
   *
   * For nodes list, we consider domain which envelopes contour. Thickness of
   * domain is typically taken as size of horizon + 2 * mesh size.
   *
   * For element list, we look for elements which intersect the contour.
   *
   * @param cd Rectangle defining the contour
   * @param tol Thickness for node list search
   * @param nodes Pointer to ids of nodes
   * @param elements Pointer to ids of elements
   */
  void
  listElemsAndNodesInDomain(const std::pair<util::Point3, util::Point3> &cd,
                            const double &tol, std::vector<size_t> *nodes,
                            std::vector<size_t> *elements);

  /*!
   * @brief Interpolates displacement and velocity at given point in triangle
   *
   * Triangle is specified by global ids of nodes. We first check if point is
   * within the triangle \f$ T \f$. If yes then we use TriElem::getShapes to
   * compute the shape functions \f$ N_1, N_2, N_3 \f$ at the given point.
   * Using shape function we
   * can interpolate displacement and velocity as follows:
   * \f[ u = \sum_{i=1}^3 N_i u^i, \quad v = \sum_{i=1}^3 N_i v^i \f]
   * where \f$ u^i, v^i \f$ are displacement and velocity of vertex \f$ i \f$
   * of triangle.
   *
   * @param p Point at which we want to interpolate
   * @param up Displacement at the point p
   * @param vp Velocity at the point p
   * @param ids Global ids of vertices of triangle
   * @param status True if point is found in the triangle. Otherwise false.
   */
  bool triCheckAndInterpolateUV(const util::Point3 &p, util::Point3 &up,
                                util::Point3 &vp,
                                const std::vector<size_t> &ids);

  /*!
   * @brief Interpolates displacement and velocity at given point
   *
   * 1. For search over nodes/elements, we use list of nodes and elements which
   * are created in listElemsAndNodesInDomain.
   *
   * 2. If element-node connectivity is not available, this method uses
   * piecewise constant interpolation by searching for node closest to the
   * point and assigning displacement and velocity of that node to the point.
   *
   * 3. If element-node connectivity is available, then we find the element
   * containing point. If element type is triangle then we use inverse
   * mapping from given triangle to reference triangle to compute shape
   * functions, see triCheckAndInterpolateUV.
   *
   * 4. If element type is quadrangle, then we split quadrangle in two
   * elements. Suppose \f$ v^1, v^2, v^3, v^4 \f$ are vertices of quadrangle,
   * then we consider triangle \f$ T_1 = (v^1, v^2, v^3) \f$ and \f$ T_2 =
   * (v^1, v^3, v^4) \f$ and find if point belongs to \f$ T_1 \f$ or \f$ T_2
   * \f$ and proceed similar to the case of triangle element (see
   * triCheckAndInterpolateUV).
   *
   * @param p Point at which we want to interpolate
   * @param up Displacement at the point p
   * @param vp Velocity at the point p
   * @param nodes Pointer to ids of nodes to perform search
   * @param elements Pointer to ids of elements to perform search
   */
  void interpolateUV(const util::Point3 &p, util::Point3 &up, util::Point3 &vp,
                     const std::vector<size_t> *nodes,
                     const std::vector<size_t> *elements);

  /*!
   * @brief Computes contribution to energy into crack from the quadrature
   * point on contour
   *
   * Computes \f[ \frac{1}{2} \rho |\dot{u}(t)|^2 dt + \bar{W}(x; u(t)) \f]
   * at given point, see computeJIntegral for more details.
   *
   * @param p Point
   * @param nodes Pointer to ids of nodes to perform search
   * @param elements Pointer to ids of elements to perform search
   * @return energy Energy contribution
   */
  double getContourContribJInt(const util::Point3 &p,
                     const std::vector<size_t> *nodes,
                     const std::vector<size_t> *elements);

  /** @}*/

  /*! @brief Input filename for compute instruction */
  std::string d_inpFilename;

  /*! @brief List of compute sets */
  std::vector<tools::pp::InstructionData> d_computeData;

  /*! @brief Current active simulation output step */
  int d_nOut;

  /*! @brief Current active compute set */
  size_t d_nC;

  /*! @brief Current active compute data */
  tools::pp::InstructionData *d_currentData;

  /*! @brief Current active output filename to save postprocessed data */
  std::string d_outFilename;

  /*! @brief Output path where postprocessed files are created */
  std::string d_outPath;

  /*! @brief Initial tag for each postprocessing output */
  std::string d_outPreTag;

  /*! @brief Simulation input filename
   *
   * This filename is same input file used in running simulation
   */
  std::string d_simInpFilename;

  /*!
   * @brief Mesh filename to get displacement and velocity (results of
   * simulation)
   */
  std::string d_simOutFilename;

  /*! @brief Displacement of nodes */
  std::vector<util::Point3> d_u;

  /*! @brief Velocity of nodes */
  std::vector<util::Point3> d_v;

  /*! @brief Model deck */
  inp::ModelDeck *d_modelDeck_p;

  /*! @brief Output deck */
  inp::OutputDeck *d_outputDeck_p;

  /*! @brief Fracture deck */
  inp::FractureDeck *d_fractureDeck_p;

  /*! @brief Material deck */
  inp::MaterialDeck *d_matDeck_p;

  /*! @brief Pointer to Mesh object */
  fe::Mesh *d_mesh_p;

  /*! @brief Pointer to Fracture object */
  geometry::Fracture *d_fracture_p;

  /*! @brief Pointer to Neighbor object */
  geometry::Neighbor *d_neighbor_p;

  /*! @brief Pointer to Input object */
  inp::Input *d_input_p;

  /*! @brief Pointer to Material object */
  material::pd::Material *d_material_p;

  /*! @brief List of fracture objects (required for crack tip calculation) */
  std::vector<geometry::Fracture *> d_fractureSet;
};

} // namespace pp

} // namespace tools

#endif // TOOLS_PP_COMPUTE_H