// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#ifndef INPUT_H
#define INPUT_H

#include <string>

//! Collection of methods performing read and write operations
namespace inp {

// forward declarations of decks
struct FractureDeck;
struct GeometryDeck;
struct InitialConditionDeck;
struct InteriorFlagsDeck;
struct LoadingDeck;
struct MaterialDeck;
struct NeighborDeck;
struct OutputDeck;
struct PolicyDeck;
struct ModelDeck;
struct SolverDeck;

/*! @brief Class to read and store input data */
class Input {

public:
  /*!
   * @brief Constructor
   * @param Filename of input file
   */
  explicit Input(const std::string& filename);

  /**
   * \defgroup Accessor methods
   */
  /**@{*/

  /*!
   * @brief Return the pointer to fracture deck
   * @return Pointer to FractureDeck
   */
  inp::FractureDeck *getFractureDeck();

  /*!
   * @brief Return the pointer to geometry deck
   * @return Pointer to GeometryDeck
   */
  inp::GeometryDeck *getGeometryDeck();

  /*!
   * @brief Return the pointer to initial condition deck
   * @return Pointer to InitialConditionDeck
   */
  inp::InitialConditionDeck *getInitialConditionDeck();

  /*!
   * @brief Return the pointer to interior flags deck
   * @return Pointer to InteriorFlagsDeck
   */
  inp::InteriorFlagsDeck *getInteriorFlagsDeck();

  /*!
   * @brief Return the pointer to loading deck
   * @return Pointer to LoadingDeck
   */
  inp::LoadingDeck *getLoadingDeck();

  /*!
   * @brief Return the pointer to material deck
   * @return Pointer to MaterialDeck
   */
  inp::MaterialDeck *getMaterialDeck();

  /*!
   * @brief Return the pointer to neighbor list deck
   * @return Pointer to NeighborDeck
   */
  inp::NeighborDeck *getNeighborDeck();

  /*!
   * @brief Return the pointer to output deck
   * @return Pointer to OutputDeck
   */
  inp::OutputDeck *getOutputDeck();

  /*!
   * @brief Return the pointer to policy info deck
   * @return Pointer to PolicyDeck
   */
  inp::PolicyDeck *getPolicyDeck();

  /*!
   * @brief Return the pointer to problem deck
   * @return Pointer to ProblemDeck
   */
  inp::ModelDeck *getModelDeck();

  /*!
   * @brief Return the pointer to solver deck
   * @return Pointer to SolverDeck
   */
  inp::SolverDeck *getSolverDeck();

  /** @}*/

  /**
   * \defgroup Accessor methods which give a peak into the data of higher
   * level data members
   */
  /**@{*/

  /*!
   * @brief Return the name of spatial discretization.
   * Return value can be of four kind: "finite_difference",
   * "weak_finite_element", "nodal_finite_element", "truss_finite_element"
   * @return Tag of spatial discretization
   */
   const std::string getSpatialDiscretization();

  /** @}*/

private:
  /**
   * \defgroup Setter methods. Reads input file and fills the higher level
   * objects
   */
  /**@{*/

  /*!
   * @brief Read data into fracture deck and store its pointer
   */
  void setFractureDeck();

  /*!
   * @brief Read data into geometry deck and store its pointer
   */
  void setGeometryDeck();

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
   * @brief Read data into material deck and store its pointer
   */
  void setMaterialDeck();

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
   * @brief Read data into problem deck and store its pointer
   */
  void setModelDeck();

  /*!
   * @brief Read data into solver deck and store its pointer
   */
  void setSolverDeck();

  /** @}*/

  /**
   * \defgroup Internal data
   */
  /**@{*/

  /*! @brief Name of input file */
  std::string d_inputFilename;

  /** @}*/

  /**
   * \defgroup Modular struct data containing input data
   */
  /**@{*/

  /*! @brief Pointer to deck holding fracture related data.
   * E.g. pre-crack location and orientation, crack-path update frequency,
   * etc
   */
  inp::FractureDeck *d_fractureDeck_p;

  /*! @brief Pointer to deck holding geometry related data.
   * E.g. dimension, discretization type, mesh file, etc
   */
  inp::GeometryDeck *d_geometryDeck_p;

  /*! @brief Pointer to deck holding initial condition related data.
   * E.g. initial condition function type for velocity and displacement,
   * parameters to compute initial condition function, projection method or
   * interpolation method, etc
   */
  inp::InitialConditionDeck *d_initialConditionDeck_p;

  /*! @brief Pointer to deck holding interior flags information */
  inp::InteriorFlagsDeck *d_interiorFlagsDeck_p;

  /*! @brief Pointer to deck holding loading related data.
   * E.g. displacement loading information, force loading information, etc
   */
  inp::LoadingDeck *d_loadingDeck_p;

  /*! @brief Pointer to deck holding material related data.
   * E.g. type of material, influence function information, parameters, etc
   */
  inp::MaterialDeck *d_materialDeck_p;

  /*! @brief Pointer to deck holding neighbor list related data.
   * E.g. factor of safety, volume correction, etc
   */
  inp::NeighborDeck *d_neighborDeck_p;

  /*! @brief Pointer to deck holding output related data.
   * E.g. output frequency, output file format, output element-node
   * connectivity flag, etc
   */
  inp::OutputDeck *d_outputDeck_p;

  /*! @brief Pointer to deck holding policy related data.
   * E.g. level of restriction in memory allocation, etc
   */
  inp::PolicyDeck *d_policyDeck_p;

  /*! @brief Pointer to deck holding problem related data.
   * E.g. type of simulation (central-difference, velocity-verlet, implicit)
   * etc
   */
  inp::ModelDeck *d_modelDeck_p;

  /*! @brief Pointer to deck holding solver related data.
   * E.g. solver parameters like tolerance, maximum iterations, solver type,
   * etc
   */
  inp::SolverDeck *d_solverDeck_p;

  /** @}*/
};

} // namespace inp

#endif // INPUT_H
