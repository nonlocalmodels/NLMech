// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef INPUT_H
#define INPUT_H

#include <string>

/*!
 * @brief Collection of methods and database related to input
 *
 * This namespace provides methods and data members specific to reading input
 * data. It provides struct definitions to read and store data. Each struct
 * is unique and is used in initialization of higher level objects such as
 * fe::Mesh, geometry::Fracture, material::Material, loading::Loading, etc.
 *
 * The namespace consists of Input and Policy member classes. Input class is
 * the main class responsible of creating struct data and reading input file
 * into those structs.
 *
 * @sa Input, Policy
 */
namespace inp {

// forward declarations of decks
struct FractureDeck;
struct InitialConditionDeck;
struct InteriorFlagsDeck;
struct LoadingDeck;
struct MassMatrixDeck;
struct MaterialDeck;
struct MeshDeck;
struct ModelDeck;
struct NeighborDeck;
struct OutputDeck;
struct PolicyDeck;
struct QuadratureDeck;
struct RestartDeck;
struct SolverDeck;

/**
 * \defgroup Input Input
 *
 * Group which reads and stores input data
 */
/**@{*/

/*!
 * @brief A class to read input file
 *
 * In this class we create struct data types, and read input file and store
 * in the respective structs. The class depends on the YAML library.
 */
class Input {

public:
  /*!
   * @brief Constructor
   * @param filename of input file
   */
  explicit Input(const std::string &filename);

  /**
   * @name Accessor methods
   */
  /**@{*/

  /*!
   * @brief Get the pointer to fracture deck
   * @return Pointer to FractureDeck
   */
  inp::FractureDeck *getFractureDeck();

  /*!
   * @brief Get the pointer to initial condition deck
   * @return Pointer to InitialConditionDeck
   */
  inp::InitialConditionDeck *getInitialConditionDeck();

  /*!
   * @brief Get the pointer to interior flags deck
   * @return Pointer to InteriorFlagsDeck
   */
  inp::InteriorFlagsDeck *getInteriorFlagsDeck();

  /*!
   * @brief Get the pointer to loading deck
   * @return Pointer to LoadingDeck
   */
  inp::LoadingDeck *getLoadingDeck();

  /*!
   * @brief Get the pointer to mass matrix deck
   * @return Pointer to MassMatrixDeck
   */
  inp::MassMatrixDeck *getMassMatrixDeck();

  /*!
   * @brief Get the pointer to material deck
   * @return Pointer to MaterialDeck
   */
  inp::MaterialDeck *getMaterialDeck();

  /*!
   * @brief Get the pointer to mesh deck
   * @return Pointer to GeometryDeck
   */
  inp::MeshDeck *getMeshDeck();

  /*!
   * @brief Get the pointer to model deck
   * @return Pointer to ModelDeck
   */
  inp::ModelDeck *getModelDeck();

  /*!
   * @brief Get the pointer to neighbor list deck
   * @return Pointer to NeighborDeck
   */
  inp::NeighborDeck *getNeighborDeck();

  /*!
   * @brief Get the pointer to output deck
   * @return Pointer to OutputDeck
   */
  inp::OutputDeck *getOutputDeck();

  /*!
   * @brief Get the pointer to policy info deck
   * @return Pointer to PolicyDeck
   */
  inp::PolicyDeck *getPolicyDeck();

  /*!
   * @brief Get the pointer to quadrature deck
   * @return Pointer to QuadratureDeck
   */
  inp::QuadratureDeck *getQuadratureDeck();

  /*!
   * @brief Get the pointer to restart deck
   * @return Pointer to RestartDeck
   */
  inp::RestartDeck *getRestartDeck();

  /*!
   * @brief Get the pointer to solver deck
   * @return Pointer to SolverDeck
   */
  inp::SolverDeck *getSolverDeck();


  /*!
   * @brief Get the name of spatial discretization
   *
   * Return value can be of four kind:
   * - "" (none)
   * - finite_difference
   * - weak_finite_element
   * - nodal_finite_element
   * - truss_finite_element
   *
   * @return Tag
   */
  const std::string getSpatialDiscretization();

  /** @}*/

private:
  /**
   * @name Setter methods
   *
   * Reads input file into the respective structs
   */
  /**@{*/

  /*!
   * @brief Read data into fracture deck and store its pointer
   */
  void setFractureDeck();

  /*!
   * @brief Read data into initial condition deck and store its pointer
   */
  void setInitialConditionDeck();

  /*!
   * @brief Read data into interior flags deck and store its pointer
   */
  void setInteriorFlagsDeck();

  /*!
   * @brief Read data into loading deck and store its pointer
   */
  void setLoadingDeck();

  /*!
   * @brief Read data into mass matrix deck and store its pointer
   */
  void setMassMatrixDeck();

  /*!
   * @brief Read data into material deck and store its pointer
   */
  void setMaterialDeck();

  /*!
   * @brief Read data into mesh deck and store its pointer
   */
  void setMeshDeck();

  /*!
   * @brief Read data into model deck and store its pointer
   */
  void setModelDeck();

  /*!
   * @brief Read data into neighbor list deck and store its pointer
   */
  void setNeighborDeck();

  /*!
   * @brief Read data into output deck and store its pointer
   */
  void setOutputDeck();

  /*!
   * @brief Read data into policy deck and store its pointer
   */
  void setPolicyDeck();

  /*!
   * @brief Read data into quadrature deck and store its pointer
   */
  void setQuadratureDeck();

  /*!
   * @brief Read data into restart deck and store its pointer
   */
  void setRestartDeck();

  /*!
   * @brief Read data into solver deck and store its pointer
   */
  void setSolverDeck();

  /** @}*/

  /**
   * @name Internal data
   */
  /**@{*/

  /*! @brief Name of input file */
  std::string d_inputFilename;

  /** @}*/

  /**
   * @name Struct data
   */
  /**@{*/

  /*!
   * @brief Pointer to deck holding fracture related data
   *
   * E.g. pre-crack location and orientation, crack-path update frequency,
   * etc
   */
  inp::FractureDeck *d_fractureDeck_p;

  /*!
   * @brief Pointer to deck holding initial condition related data
   *
   * E.g. initial condition function type for velocity and displacement,
   * parameters to compute initial condition function, projection method or
   * interpolation method, etc
   */
  inp::InitialConditionDeck *d_initialConditionDeck_p;

  /*!
   * @brief Pointer to deck holding interior flags information
   */
  inp::InteriorFlagsDeck *d_interiorFlagsDeck_p;

  /*!
   * @brief Pointer to deck holding loading related data.
   *
   * E.g. displacement loading information, force loading information, etc
   */
  inp::LoadingDeck *d_loadingDeck_p;

  /*!
   * @brief Pointer to deck holding mass matrix calculation related data
   *
   * E.g. mass matrix approximation type
   */
  inp::MassMatrixDeck *d_massMatrixDeck_p;

  /*!
   * @brief Pointer to deck holding material related data
   *
   * E.g. type of material, influence function information, parameters, etc
   */
  inp::MaterialDeck *d_materialDeck_p;

  /*!
   * @brief Pointer to deck holding geometry related data
   *
   * E.g. dimension, discretization type, mesh file, etc
   */
  inp::MeshDeck *d_meshDeck_p;

  /*!
   * @brief Pointer to deck holding problem related data
   *
   * E.g. type of simulation (central-difference, velocity-verlet, implicit)
   * etc
   */
  inp::ModelDeck *d_modelDeck_p;

  /*!
   * @brief Pointer to deck holding neighbor list related data
   *
   * E.g. factor of safety, volume correction, etc
   */
  inp::NeighborDeck *d_neighborDeck_p;

  /*!
   * @brief Pointer to deck holding output related data
   *
   * E.g. output frequency, output file format, output element-node
   * connectivity flag, etc
   */
  inp::OutputDeck *d_outputDeck_p;

  /*!
   * @brief Pointer to deck holding policy related data
   *
   * E.g. level of restriction in memory allocation, etc
   */
  inp::PolicyDeck *d_policyDeck_p;

  /*!
   * @brief Pointer to deck holding quadrature point approximation related
   * data
   *
   * E.g. order of approximation, order of approximation for mass matrix
   */
  inp::QuadratureDeck *d_quadratureDeck_p;

  /*!
   * @brief Pointer to deck holding restart related data such as restart
   * filename and restart time step
   */
  inp::RestartDeck *d_restartDeck_p;

  /*!
   * @brief Pointer to deck holding solver related data
   *
   * E.g. solver parameters like tolerance, maximum iterations, solver type,
   * etc
   */
  inp::SolverDeck *d_solverDeck_p;

  /** @}*/
};

/** @}*/

} // namespace inp

#endif // INPUT_H
