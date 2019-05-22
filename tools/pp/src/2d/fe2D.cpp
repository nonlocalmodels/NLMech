// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "fe2D.h"
#include "inp/input.h"          // definition of Input
#include "inp/decks/materialDeck.h"       // definition of MaterialDeck
#include "inp/decks/modelDeck.h"          // definition of ModelDeck
#include "inp/decks/outputDeck.h"         // definition of OutputDeck
#include "inp/policy.h"         // definition of Policy
#include "rw/reader.h"          // definition of readVtuFileRestart
#include "rw/writer.h"          // definition of VtkWriterInterface
#include "util/compare.h"       // compare real numbers
#include "util/feElementDefs.h" // definition of fe element type
#include "util/point.h"         // definition of Point3
#include "util/matrix.h"        // definition of SymMatrix3
#include "fe/triElem.h"         // definition of TriElem
#include "fe/quadElem.h"        // definition of QuadElem
#include "fe/mesh.h"            // definition of Mesh
#include "material/pdMaterial.h"          // definition of Material
#include "util/utilGeom.h"      // definition of isPointInsideRectangle
#include <cmath>
#include <iostream>
#include <yaml-cpp/yaml.h>      // YAML reader

namespace {

struct InstructionData {
  std::string d_tagFilename;

  bool d_scaleUOut;
  double d_scaleU;

  bool d_damageAtNodes;

  bool d_outOnlyNodes;

  bool d_removeElements;

  bool d_markVAsZero;
  bool d_markVInRectGiven;
  std::pair<util::Point3, util::Point3> d_markVRect;
  std::vector<util::Point3> d_markVPts;
  bool d_markVPtsAreInCurrentConfig;
  std::vector<size_t> d_markVNodes;

  bool d_computeStrain;
  bool d_outStrainAtElements;

  bool d_magStrainTensor;
  std::string d_magStrainComp;
  std::vector<std::pair<size_t, double>> d_markMagStrainCells;

  bool d_symmetrizeV;
  bool d_combineMarkV;
  std::string d_symmAxis;
  double d_symmLine;

  bool d_crackTip;

  InstructionData()
      : d_scaleUOut(false), d_scaleU(1.), d_damageAtNodes(false),
        d_outOnlyNodes(true), d_removeElements(false), d_markVAsZero(false),
        d_markVInRectGiven(false),
        d_markVRect(std::make_pair(util::Point3(), util::Point3())),
        d_markVPtsAreInCurrentConfig(false), d_computeStrain(false),
        d_outStrainAtElements(false), d_magStrainTensor(false),
        d_symmetrizeV(false), d_combineMarkV(false), d_symmLine(0.),
        d_crackTip(false){};
};

void readInputFile(YAML::Node config, const std::string &set,
                   InstructionData *data) {

  if (config["Compute"][set]["Tag_Filename"]) {
    data->d_tagFilename = "_" +
        config["Compute"][set]["Tag_Filename"].as<std::string>();
  } else {
    // we work with default filename of form
    // pp_set_1_time_step_0.vtu
    data->d_tagFilename = "_" + set;
  }

  //
  // Output only nodes
  //
  if (config["Compute"][set]["Output_Only_Nodes"])
    data->d_outOnlyNodes = config["Compute"][set]["Output_Only_Nodes"]
        .as<bool>();

  //
  // Scale displacement
  //
  if (config["Compute"][set]["Scale_U_Ouptut"]) {
    data->d_scaleUOut = true;
    data->d_scaleU = config["Compute"][set]["Scale_U_Ouptut"].as<double>();
  }

  //
  // Compute damage at nodes
  //
  if (config["Compute"][set]["Damage_Z"])
    data->d_damageAtNodes = config["Compute"][set]["Damage_Z"].as<bool>();

  //
  // Mark velocity as zero
  //
  if (config["Compute"][set]["Mark_V_0"]) {
    data->d_markVAsZero = true;

    // read rectangle information
    data->d_markVInRectGiven = false;
    if (config["Compute"][set]["Mark_V_0"]["Rectangle"]) {

      data->d_markVInRectGiven = true;
      std::vector<double> locs;
      for (auto j : config["Compute"][set]["Mark_V_0"]["Rectangle"])
        locs.push_back(j.as<double>());

      if (locs.size() != 4) {
        std::cerr << "Error: Check Rectangle data for mark V as zero task.\n";
        exit(1);
      }

      data->d_markVRect.first = util::Point3(locs[0], locs[1], 0.);
      data->d_markVRect.second = util::Point3(locs[2], locs[3], 0.);
    }

    // read points
    if (config["Compute"][set]["Mark_V_0"]["Points"])
      for (auto &&pt : config["Compute"][set]["Mark_V_0"]["Points"]) {
        std::vector<double> locs;
        for (auto j : pt)
          locs.push_back(j.as<double>());

        if (locs.size() == 2)
          locs.push_back(0.0);

        data->d_markVPts.emplace_back(locs[0], locs[1], locs[2]);
      }

    if (!data->d_markVPts.empty() &&
        config["Compute"][set]["Mark_V_0"]["Points_Current_Config"])
      data->d_markVPtsAreInCurrentConfig =
          config["Compute"][set]["Mark_V_0"]["Points_Current_Config"]
              .as<bool>();

    if (config["Compute"][set]["Mark_V_0"]["Nodes"])
      for (auto &&pt : config["Compute"][set]["Mark_V_0"]["Nodes"])
        data->d_markVNodes.push_back(pt.as<size_t>());

    if (data->d_markVPts.empty() && data->d_markVNodes.empty() &&
        !data->d_markVInRectGiven)
      data->d_markVAsZero = false;
  }

  //
  // Compute strain and stress
  //
  if (config["Compute"][set]["Strain_Stress"]) {
    data->d_computeStrain = true;
    if (config["Compute"][set]["Strain_Stress"]["Output_At_Elements"])
      data->d_outStrainAtElements =
          config["Compute"][set]["Strain_Stress"]["Output_At_Elements"]
              .as<bool>();
  }

  //
  // Compute magnitude of strain and stress
  //
  if (config["Compute"][set]["Magnitude_Strain_Tensor"]) {
    data->d_magStrainTensor = true;

    // if d_computeStrain is set to false, then set it to true
    data->d_computeStrain = true;

    if (config["Compute"][set]["Magnitude_Strain_Tensor"]["Component"])
      data->d_magStrainComp =
          (config["Compute"][set]["Magnitude_Strain_Tensor"]["Component"])
              .as<std::string>();

    if (config["Compute"][set]["Magnitude_Strain_Tensor"]["Cells"])
      for (auto &&i :
           config["Compute"][set]["Magnitude_Strain_Tensor"]["Cells"])
        data->d_markMagStrainCells.emplace_back(i[0].as<size_t>(),
                                                i[1].as<double>());
  }

  //
  // Symmetrize velocity field
  //
  if (config["Compute"][set]["Symmetrize_V"]) {
    data->d_symmetrizeV = true;

    // check if we need to combine this with mark_v operation
    if (config["Compute"][set]["Symmetrize_V"]["Combine_Mark_V_0"])
      data->d_combineMarkV =
          config["Compute"][set]["Symmetrize_V"]["Combine_Mark_V_0"].as<bool>();

    if (config["Compute"][set]["Symmetrize_V"]["Axis"])
      data->d_symmAxis =
          config["Compute"][set]["Symmetrize_V"]["Axis"].as<std::string>();
    else {
      std::cerr << "Error: Need Axis of symmetry for symmetrization of "
                   "velocity.\n";
      exit(1);
    }

    if (config["Compute"][set]["Symmetrize_V"]["Axis_Line"])
      data->d_symmLine =
          config["Compute"][set]["Symmetrize_V"]["Axis_Line"].as<double>();
    else {
      std::cerr << "Error: Need location of symmetry axis for symmetrization "
                   "of velocity.\n";
      exit(1);
    }
  }

  //
  // Crack tip calculation
  //
  if (config["Compute"][set]["Crack_Tip"]) {
    data->d_crackTip = config["Compute"][set]["Crack_Tip"].as<bool>();
    if (data->d_crackTip)
      data->d_damageAtNodes = true;
  }
}

size_t findNode(const util::Point3 &x, const std::vector<util::Point3>
  *nodes, const std::vector<util::Point3> *u = nullptr) {

  long int i_found = -1;
  double dist = 1000.0;
  for (size_t i = 0; i < nodes->size(); i++)
    if (u == nullptr) {
      if (util::compare::definitelyLessThan(x.dist((*nodes)[i]), dist)) {
        dist = x.dist((*nodes)[i]);
        i_found = i;
      }
    } else {
      if (util::compare::definitelyLessThan(x.dist((*nodes)[i] + (*u)[i]),
        dist)) {
        dist = x.dist((*nodes)[i] + (*u)[i]);
        i_found = i;
      }
    }

  if (i_found < 0) {
    std::cerr << "Error: Could not find the node.\n";
    exit(1);
  }

  return size_t(i_found);
}

} // namespace

void tools::pp::fe2D(const std::string &filename) {
  // read YAML file
  auto config = YAML::LoadFile(filename);

  // get simulation filename
  auto source_path = config["Source_Path"].as<std::string>();
  auto sim_filename =
      source_path + "/" + config["Simulation_Input_File"].as<std::string>();

  // read input data
  std::cout << "PP_fe2D: Reading simulation input file.\n";
  auto *deck = new inp::Input(filename);
  auto model_deck = deck->getModelDeck();
  auto output_deck = deck->getOutputDeck();

  // get policy deck
  auto policy = inp::Policy::getInstance(deck->getPolicyDeck());

  size_t process_single_file = 0; // default
  if (config["Override_Simulation_Input"]["Process_Single_File"])
    process_single_file =
        config["Override_Simulation_Input"]["Process_Single_File"].as<size_t>();

  std::string filename_to_read;
  if (config["Filename_To_Read"])
    filename_to_read = config["Filename_To_Read"].as<std::string>();
  else {
    std::cerr << "Error: Simulation results filename is not provided.\n";
    exit(1);
  }

  // get output path directory
  std::string out_path = "./"; // default
  if (config["Output"]["Path"])
    out_path = config["Output"]["Path"].as<std::string>();
  std::string out_filename = "pp"; // default
  if (config["Output"]["Filename"])
    out_filename = config["Output"]["Filename"].as<std::string>();

  // over ride final time step
  if (config["Override_Simulation_Input"]["Time_Steps"])
    model_deck->d_Nt = config["Override_Simulation_Input"]["Time_Steps"]
        .as<size_t>();
  size_t start_file = 1;
  if (process_single_file > 0)
    start_file = process_single_file;
  size_t end_file = model_deck->d_Nt / output_deck->d_dtOut;
  if (process_single_file > 0)
    end_file = process_single_file;

  // create mesh
  std::cout << "PP_fe2D: Creating mesh.\n";
  auto *mesh = new fe::Mesh(deck->getMeshDeck());

  // loop over output files
  for (size_t out=start_file; out<=end_file; out++) {
    // get number of compute set
    auto num_compute = config["Compute"]["Sets"].as<size_t>();

    // append path to filename
    auto sim_results_file =
        source_path + "/" + filename_to_read + "_" + std::to_string(out);

    // just read one file and do operations and move to the next compute set
    std::vector<util::Point3> u;
    std::vector<util::Point3> v;

    // get displacement and velocity from the simulation
    rw::reader::readVtuFileRestart(sim_results_file, &u, &v);

    // loop over compute sets and do as instructed in input file
    for (size_t c = 1; c <= num_compute; c++) {

      // get string to read instruction
      auto set = "Set_" + std::to_string(c);

      // read file for instructions
      auto data = InstructionData();
      readInputFile(config, set, &data);

      // keep track of whether mesh has been written to output file
      auto mesh_appended = false;

      // open a output vtu file
      auto writer = rw::writer::VtkWriterInterface(out_path + out_filename +
          data.d_tagFilename);

      //
      // operation : Scale displacement and write to output mesh
      //
      // Steps: Create displacement vector and scale it and write to mesh
      //
      if (data.d_scaleUOut) {

        std::vector<util::Point3> u_temp(mesh->getNumNodes(), util::Point3());
        for (size_t i = 0; i < mesh->getNumNodes(); i++)
          u_temp[i] =
              util::Point3(data.d_scaleU * u[i].d_x, data.d_scaleU * u[i].d_y,
                           data.d_scaleU * u[i].d_z);

        // append mesh (check if only nodes need to be written)
        if (data.d_outOnlyNodes)
          writer.appendNodes(mesh->getNodesP(), &u_temp);
        else
          writer.appendMesh(mesh->getNodesP(), mesh->getElementType(),
                            mesh->getElementConnectivitiesP(), &u_temp);

        mesh_appended = true;

        // append original displacement
        writer.appendPointData("Displacement", &u);

        // append velocity
        writer.appendPointData("Velocity", &v);
      }

      //
      // operation : Mark velocity in region as zero
      //
      // Steps: Read initial, current, and velocity, and modify
      // the velocity of nodes in the rectangle
      //
      // vector to hold new velocity
      std::vector<util::Point3> v_mark;

      if (data.d_markVAsZero) {
        v_mark = v;
        if (data.d_markVInRectGiven)
          for (size_t i = 0; i < mesh->getNumNodes(); i++)
            if (util::geometry::isPointInsideRectangle(
                    mesh->getNode(i), data.d_markVRect.first[0],
                    data.d_markVRect.second[0], data.d_markVRect.first[1],
                    data.d_markVRect.second[1]))
              v_mark[i] = util::Point3();

        if (!data.d_markVPts.empty())
          for (auto x : data.d_markVPts) {
            size_t i_found;
            if (data.d_markVPtsAreInCurrentConfig)
              i_found = findNode(x, mesh->getNodesP(), &u);
            else
              i_found = findNode(x, mesh->getNodesP());

            // modify v_new
            v_mark[i_found] = util::Point3();
          }

        if (!data.d_markVNodes.empty())
          for (auto i : data.d_markVNodes)
            v_mark[i] = util::Point3();

        if (!mesh_appended) {
          if (data.d_outOnlyNodes)
            writer.appendNodes(mesh->getNodesP(), &u);
          else
            writer.appendMesh(mesh->getNodesP(), mesh->getElementType(),
                              mesh->getElementConnectivitiesP(), &u);

          mesh_appended = true;
        }

        // append velocity
        writer.appendPointData("Mark_Velocity", &v_mark);
      }

      //
      // operation : Symmetrize velocity
      //
      // Steps: Create new velocity data by symmetrizing it along the
      // specified line of symmetry
      //
      if (data.d_symmetrizeV) {

        // use v_mark for modification
        // Reset it to current velocity if we are not combining this
        // with mark_v operation

        if (!data.d_combineMarkV)
          v_mark = v;

        for (size_t i = 0; i < mesh->getNumNodes(); i++) {

          auto x = mesh->getNode(i);

          if (data.d_symmAxis == "y" &&
              util::compare::definitelyLessThan(x.d_x, data.d_symmLine + 1.0E-8))
            continue;

          if (data.d_symmAxis == "x" &&
              util::compare::definitelyLessThan(x.d_y, data.d_symmLine + 1.0E-8))
            continue;

          // find the coordinate of point from where we want to copy
          // the velocity at this node.
          // Mirror image of point
          auto search_x = x;
          if (data.d_symmAxis == "y")
            search_x.d_x = data.d_symmLine - (x.d_x - data.d_symmLine);
          if (data.d_symmAxis == "x")
            search_x.d_y = data.d_symmLine - (x.d_y - data.d_symmLine);

          // search for node at search_x and obtain velocity
          size_t i_found = findNode(search_x, mesh->getNodesP());

          // write velocity
          v_mark[i] = v_mark[i_found];
          if (data.d_symmAxis == "y")
            v_mark[i].d_x *= -1.;
          if (data.d_symmAxis == "x")
            v_mark[i].d_y *= -1.;
        }

        if (!mesh_appended) {
          if (data.d_outOnlyNodes)
            writer.appendNodes(mesh->getNodesP(), &u);
          else
            writer.appendMesh(mesh->getNodesP(), mesh->getElementType(),
                              mesh->getElementConnectivitiesP(), &u);

          mesh_appended = true;
        }

        // append velocity
        writer.appendPointData("Symm_Velocity", &v_mark);
      }

      //
      // operation : Compute strain
      //
      // Steps: Call CCMFE class and compute strain and stress
      //
      std::vector<util::SymMatrix3> strain;
      std::vector<util::SymMatrix3> stress;
      if (data.d_computeStrain) {

        auto *material = new material::pd::Material(
            deck->getMaterialDeck(), model_deck->d_dim, model_deck->d_horizon);
        auto mat_deck = material->getMaterialDeck();

        //
        // compute strain and stress
        //
        strain.resize(mesh->getNumElements());
        stress.resize(mesh->getNumElements());

        // get Quadrature
        fe::BaseElem *quad;
        if (mesh->getElementType() == util::vtk_type_triangle)
          quad = new fe::TriElem(1);
        else if (mesh->getElementType() == util::vtk_type_quad)
          quad = new fe::QuadElem(1);
        else {

          std::cerr << "Error: Can not compute strain/stress as the element "
                       "type is not implemented.\n";
          exit(1);
        }

        for (size_t e = 0; e < mesh->getNumElements(); e++) {
          auto ssn = util::SymMatrix3();
          auto sss = util::SymMatrix3();

          // get ids of nodes of element, coordinate of nodes, 1st order quad
          // data, and first quad data
          auto id_nds = mesh->getElementConnectivity(e);
          auto nds = mesh->getElementConnectivityNodes(e);
          auto qds = quad->getQuadPoints(nds);
          auto qd0 = qds[0];

          // compute strain in xy plane
          for (size_t i = 0; i < id_nds.size(); i++) {
            auto id = id_nds[i];
            auto ui = u[id];

            ssn.d_xx += ui.d_x * qd0.d_derShapes[i][0];
            ssn.d_yy += ui.d_y * qd0.d_derShapes[i][1];
            ssn.d_xy += 0.5 * ui.d_x * qd0.d_derShapes[i][1];
            ssn.d_xy += 0.5 * ui.d_y * qd0.d_derShapes[i][0];
          }

          // if material deck does not have valid material properties, search
          // the properties in the pp input file
          if (mat_deck->d_matData.d_nu < 0. || mat_deck->d_matData.d_E < 0.) {
            if (config["Material"]["Poisson_Ratio"])
              mat_deck->d_matData.d_nu =
                  config["Material"]["Poisson_Ratio"].as<double>();
            else {
              std::cerr << "Error: Need Poisson ratio for strain and stress "
                           "computation.\n";
              exit(1);
            }
            if (config["Material"]["E"])
              mat_deck->d_matData.d_E = config["Material"]["E"].as<double>();
            else {
              std::cerr << "Error: Need Young's modulus for strain and stress "
                           "computation.\n";
              exit(1);
            }

            // compute lambda and mu
            mat_deck->d_matData.d_lambda = mat_deck->d_matData.toLambdaE(
                mat_deck->d_matData.d_E, mat_deck->d_matData.d_nu);
            mat_deck->d_matData.d_mu = mat_deck->d_matData.d_lambda;
          }

          if (mat_deck->d_isPlaneStrain)
            ssn.d_zz = -mat_deck->d_matData.d_nu * (ssn.d_xx + ssn.d_yy) /
                       (1. - mat_deck->d_matData.d_nu);

          // compute stress
          auto trace = ssn.d_xx + ssn.d_yy + ssn.d_zz;
          sss.d_xx = mat_deck->d_matData.d_lambda * trace * 1. +
                     2. * mat_deck->d_matData.d_mu * ssn.d_xx;
          sss.d_xy = 2. * mat_deck->d_matData.d_mu * ssn.d_xy;
          sss.d_xz = 2. * mat_deck->d_matData.d_mu * ssn.d_xz;

          sss.d_yy = mat_deck->d_matData.d_lambda * trace * 1. +
                     2. * mat_deck->d_matData.d_mu * ssn.d_yy;
          sss.d_yz = 2. * mat_deck->d_matData.d_mu * ssn.d_yz;
          if (!mat_deck->d_isPlaneStrain)
            sss.d_zz = mat_deck->d_matData.d_nu * (sss.d_xx + sss.d_yy);
        } // loop over elements

        //
        // output strain/stress data
        //
        if (data.d_outStrainAtElements) {

          // append mesh
          if (!mesh_appended) {
            writer.appendMesh(mesh->getNodesP(), mesh->getElementType(),
                              mesh->getElementConnectivitiesP(), &u);

            mesh_appended = true;
          }

          writer.appendCellData("Strain_Tensor", &strain);
          writer.appendCellData("Stress_Tensor", &stress);
        } else {

          // compute 1st order quad points and store them (quad points in
          // current configuration)
          std::vector<util::Point3> elem_quads =
              std::vector<util::Point3>(mesh->getNumElements(), util::Point3());

          for (size_t e = 0; e < mesh->getNumElements(); e++) {
            auto nds = mesh->getElementConnectivity(e);
            std::vector<util::Point3> nds_current(nds.size(), util::Point3());
            for (size_t j = 0; j < nds.size(); j++)
              nds_current[j] = mesh->getNode(nds[j]) + u[nds[j]];

            std::vector<fe::QuadData> qds = quad->getQuadPoints(nds_current);
            // store first quad point
            elem_quads[e] = qds[0].d_p;
          }

          // create unstructured vtk output
          std::string fname = data.d_filename + "_quads.vtu";
          auto writer1 = rw::writer::VtkWriterInterface(fname);
          writer1.appendNodes(&elem_quads);
          writer1.appendPointData("Strain_Tensor", &strain);
          writer1.appendPointData("Stress_Tensor", &stress);
          writer1.close();
        }
      }

      //
      // operation : Compute magnitude of strain
      //
      // Steps: Create new data and store magnitude of strain. Either
      // compute sum of absolute value of all components or the component
      // specified by the input file.
      //
      if (data.d_magStrainTensor) {
        std::vector<float> magS = std::vector<float>(strain.size(), 0.);
        for (size_t i = 0; i < strain.size(); i++) {
          if (data.d_magStrainComp.empty()) {
            magS[i] = std::abs(strain[i].d_xx);
            magS[i] += std::abs(strain[i].d_yy);
            magS[i] += std::abs(strain[i].d_zz);
            magS[i] += std::abs(strain[i].d_xy);
            magS[i] += std::abs(strain[i].d_xz);
            magS[i] += std::abs(strain[i].d_yz);
          } else if (data.d_magStrainComp == "xx") {
            magS[i] = std::abs(strain[i].d_xx);
          } else if (data.d_magStrainComp == "yy") {
            magS[i] = std::abs(strain[i].d_yy);
          }
        }

        // append mesh
        if (!mesh_appended) {
          writer.appendMesh(mesh->getNodesP(), mesh->getElementType(),
                            mesh->getElementConnectivitiesP(), &u);

          mesh_appended = true;
        }

        // append new cell data
        writer.appendCellData("Mag_Strain", &magS);

        // mark strains of certain cells prescribed value if asked
        for (auto cell : data.d_markMagStrainCells)
          magS[cell.first] = cell.second;

        // append new cell data
        writer.appendCellData("Mark_Mag_Strain", &magS);
      }

      //
      // operation : Compute damage at node
      //
      if (data.d_damageAtNodes) {
        auto *material = new material::pd::Material(
            deck->getMaterialDeck(), model_deck->d_dim, model_deck->d_horizon);

        std::vector<float> damage_Z(mesh->getNumNodes(), 0.);
        for (size_t i=0; i<mesh->getNumNodes(); i++) {
          auto xi = mesh->getNode(i);
          for (size_t j = 0; j < mesh->getNumNodes(); j++) {
            if (util::compare::definitelyGreaterThan(xi.dist(mesh->getNode(j)
            ), model_deck->d_horizon) || j == i)
              continue;

            auto xj = mesh->getNode(j);
            if (util::compare::definitelyGreaterThan(xj.dist(xi), 1.0E-10)) {
              auto Sr = std::abs(material->getS(xj - xi, u[j] - u[i])) /
                  material->getSc(xj.dist(xi));

              if(util::compare::definitelyLessThan(damage_Z[i], Sr))
                damage_Z[i] = Sr;
            }
          } // loop over neighbors
        } // loop over nodes

        if (!mesh_appended) {
          if (data.d_outOnlyNodes)
            writer.appendNodes(mesh->getNodesP(), &u);
          else
            writer.appendMesh(mesh->getNodesP(), mesh->getElementType(),
                              mesh->getElementConnectivitiesP(), &u);

          mesh_appended = true;
        }
        // append data to file
        writer.appendPointData("Damage_Z", &damage_Z);
      }

      //
      // operation : Compute crack tip location and velocity
      //
      if (data.d_crackTip) {

      }

      //
      // close file
      //
      if (mesh_appended)
        writer.close();
    } // loop over compute sets
  }// processing output files
}