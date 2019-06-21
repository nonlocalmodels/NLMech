// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "compute.h"
#include "fe/lineElem.h" // definition of LineElem
#include "fe/quadElem.h" // definition of QuadElem
#include "fe/triElem.h"  // definition of TriElem
#include "rw/reader.h"   // definition of readVtuFileRestart
#include "inp/policy.h"  // definition of Policy
#include "util/feElementDefs.h"     // definition of fe element type
#include "util/utilGeom.h"          // definition of isPointInsideRectangle
#include "external/csv.h"           // csv reader
#include <hpx/include/parallel_algorithm.hpp>
#include <cmath>
#include <iostream>
#include <yaml-cpp/yaml.h> // YAML reader

tools::pp::Compute::Compute(const std::string &filename)
    : d_inpFilename(filename), d_nC(0),
      d_currentData(nullptr), d_modelDeck_p(nullptr),
      d_outputDeck_p(nullptr), d_fractureDeck_p(nullptr), d_matDeck_p(nullptr),
      d_mesh_p(nullptr), d_fracture_p(nullptr), d_neighbor_p(nullptr),
      d_input_p(nullptr), d_material_p(nullptr) {

  init();

  // loop over output steps
  int N = d_modelDeck_p->d_Nt / d_outputDeck_p->d_dtOut;
  for (d_nOut = 1; d_nOut <= N; d_nOut++) {
    std::cout << "PP_fe2D: Processing output file = " << d_nOut << "\n";

    // mesh filename to read displacement and velocity
    d_simOutFilename.append(std::to_string(d_nOut) + ".vtu");

    // see if we need to read the data
    bool read_file = false;
    for (const auto &data : d_computeData)
      read_file = (d_nOut >= data.d_start && d_nOut <= data.d_end);
    if (!read_file)
      continue;

    // get displacement and velocity
    rw::reader::readVtuFileRestart(d_simOutFilename, &d_u, &d_v);

    // loop over compute sets and do as instructed in input file
    for (d_nC = 0; d_nC < d_computeData.size(); d_nC++) {
      d_currentData = &(d_computeData[d_nC]);
      // check if in the specified bound
      if (d_nOut < d_currentData->d_start || d_nOut > d_currentData->d_end)
        continue;

      std::cout << "  PP_fe2D: Processing compute set = " << d_nC + 1 << "\n";

      // filename for writing postprocessed data
      d_outFilename = d_outPreTag + "_" + d_currentData->d_tagFilename + "_" +
                      std::to_string(d_nOut);
      // writer
      rw::writer::VtkWriterInterface *writer = nullptr;

      // apply postprocessing
      transformU(writer);
      transformV(writer);
      computeStrain(writer);

      // nodal damage
      std::vector<float> Z;
      findCrackTip(&Z, writer);
      if (d_currentData->d_damageAtNodes)
        computeDamage(&Z, writer, true);
      if (!Z.empty())
        Z.shrink_to_fit();

      computeJIntegral(writer);

      // close file
      if (writer)
        writer->close();
    } // loop compute set
  } // loop simulation output
}

void tools::pp::Compute::init() {
  auto config = YAML::LoadFile(d_inpFilename);

  // get simulation filename
  auto source_path = config["Source_Path"].as<std::string>();
  d_simInpFilename =
      source_path + "/" + config["Simulation_Input_File"].as<std::string>();

  // get simulation results filename
  if (!config["Source_Path"] || !config["Filename_To_Read"]) {
    std::cerr << "Error: Either Source_Path or Filename_To_Read is not "
                 "provided in postprocessing input (yaml) file.\n";
    exit(1);
  }
  d_simOutFilename = config["Source_Path"].as<std::string>() + "/" +
                     config["Filename_To_Read"].as<std::string>() + "_";

  // read input data
  std::cout << "PP_fe2D: Reading simulation input file.\n";
  d_input_p = new inp::Input(d_simInpFilename);
  d_modelDeck_p = d_input_p->getModelDeck();
  d_outputDeck_p = d_input_p->getOutputDeck();
  d_fractureDeck_p = d_input_p->getFractureDeck();

  // get policy deck (Policy deck must be initialized)
  auto policy = inp::Policy::getInstance(d_input_p->getPolicyDeck());

  // get output path directory
  if (config["Output"]["Path"])
    d_outPath = config["Output"]["Path"].as<std::string>();
  else
    d_outPath = "./"; // default

  if (config["Output"]["Filename"])
    d_outPreTag =
        d_outPath + "/" + config["Output"]["Filename"].as<std::string>();
  else
    d_outPreTag = d_outPath + "/pp";

  // create mesh
  std::cout << "PP_fe2D: Creating mesh.\n";
  d_mesh_p = new fe::Mesh(d_input_p->getMeshDeck());

  // material deck and material
  d_material_p = new material::pd::Material(
      d_input_p->getMaterialDeck(), d_modelDeck_p->d_dim, d_modelDeck_p->d_horizon);
  d_matDeck_p = d_material_p->getMaterialDeck();
  // if material deck does not have valid material properties,
  // search the properties in the pp input file
  if (d_matDeck_p->d_matData.d_nu < 0. || d_matDeck_p->d_matData.d_E < 0.) {
    if (config["Material"]["Poisson_Ratio"])
      d_matDeck_p->d_matData.d_nu =
          config["Material"]["Poisson_Ratio"].as<double>();
    else {
      std::cerr << "Error: Need Poisson ratio for strain and stress "
                   "computation.\n";
      exit(1);
    }
    if (config["Material"]["E"])
      d_matDeck_p->d_matData.d_E = config["Material"]["E"].as<double>();
    else {
      std::cerr << "Error: Need Young's modulus for strain and stress "
                   "computation.\n";
      exit(1);
    }

    // compute lambda and mu
    d_matDeck_p->d_matData.d_lambda = d_matDeck_p->d_matData.toLambdaE(
        d_matDeck_p->d_matData.d_E, d_matDeck_p->d_matData.d_nu);
    d_matDeck_p->d_matData.d_mu = d_matDeck_p->d_matData.d_lambda;
  }

  // read compute instruction
  auto num_compute = config["Compute"]["Sets"].as<size_t>();
  for (size_t c = 0; c < num_compute; c++) {
    std::string set = "Set_" + std::to_string(c + 1);
    auto data = tools::pp::InstructionData();
    readComputeInstruction(set, &data);
    if (data.d_start == -1)
      data.d_start = 1;
    if (data.d_end == -1)
      data.d_end = d_modelDeck_p->d_Nt / d_outputDeck_p->d_dtOut;

    d_computeData.emplace_back(data);
  }

  // vector of null pointer Fracture class
  d_fractureSet = std::vector<geometry::Fracture *>(num_compute, nullptr);

  // read crack tip data for J integral calculation
  for (auto &d : d_computeData) {
    auto data = d.d_computeJInt_p;
    if (data) {
      readCrackTipData(data->d_crackTipFile, &(data->d_crackTipData));
      if (!data->d_crackTipData.empty()) {
        data->d_start = data->d_crackTipData[0].d_n;
        data->d_end =
            data->d_crackTipData[data->d_crackTipData.size() - 1].d_n;
      }
    }
  }

  // compute neighbor list (if required)
  bool compute_neighbor_list = false;
  for (const auto& d : d_computeData)
    compute_neighbor_list = (d.d_computeJInt_p != nullptr);

  if (compute_neighbor_list)
    d_neighbor_p = new geometry::Neighbor(d_modelDeck_p->d_horizon,
                                          d_input_p->getNeighborDeck(),
                                          d_mesh_p->getNodesP());
}

void tools::pp::Compute::readComputeInstruction(
    const std::string &set, tools::pp::InstructionData *data) {

  auto config = YAML::LoadFile(d_inpFilename);

  // read tag for output file
  if (config["Compute"][set]["Tag_Filename"]) {
    data->d_tagFilename =
        config["Compute"][set]["Tag_Filename"].as<std::string>();
  } else {
    // we work with default filename of form pp_set_1_0.vtu
    data->d_tagFilename = set;
  }

  // Output only nodes
  if (config["Compute"][set]["Output_Only_Nodes"])
    data->d_outOnlyNodes =
        config["Compute"][set]["Output_Only_Nodes"].as<bool>();

  // check if start and end time step are specified
  if (config["Compute"][set]["Dt_Start"])
    data->d_start = config["Compute"][set]["Dt_Start"].as<int>();
  if (config["Compute"][set]["Dt_End"])
    data->d_end = config["Compute"][set]["Dt_End"].as<int>();

  // Scale displacement
  if (config["Compute"][set]["Scale_U_Ouptut"]) {
    if (!data->d_transformU_p)
      data->d_transformU_p = new tools::pp::TransformU();
    data->d_transformU_p->d_scale =
        config["Compute"][set]["Scale_U_Ouptut"].as<double>();
  }

  // Compute damage at nodes
  if (config["Compute"][set]["Damage_Z"])
    data->d_damageAtNodes = config["Compute"][set]["Damage_Z"].as<bool>();

  // Mark velocity as zero
  if (config["Compute"][set]["Mark_V_0"]) {
    if (!data->d_transformV_p)
      data->d_transformV_p = new tools::pp::TransformVelocity();
    data->d_transformV_p->d_markVAsZero = true;

    // read rectangle information
    if (config["Compute"][set]["Mark_V_0"]["Rectangle"]) {

      data->d_transformV_p->d_markVInRectGiven = true;
      std::vector<double> locs;
      for (auto j : config["Compute"][set]["Mark_V_0"]["Rectangle"])
        locs.push_back(j.as<double>());

      if (locs.size() != 4) {
        std::cerr << "Error: Check Rectangle data for mark V as zero task.\n";
        exit(1);
      }

      data->d_transformV_p->d_markVRect.first =
          util::Point3(locs[0], locs[1], 0.);
      data->d_transformV_p->d_markVRect.second =
          util::Point3(locs[2], locs[3], 0.);
    }

    // read points
    if (config["Compute"][set]["Mark_V_0"]["Points"])
      for (auto &&pt : config["Compute"][set]["Mark_V_0"]["Points"]) {
        std::vector<double> locs;
        for (auto j : pt)
          locs.push_back(j.as<double>());

        if (locs.size() == 2)
          locs.push_back(0.0);

        data->d_transformV_p->d_markVPts.emplace_back(locs[0], locs[1],
                                                      locs[2]);
      }

    if (!data->d_transformV_p->d_markVPts.empty() &&
        config["Compute"][set]["Mark_V_0"]["Points_Current_Config"])
      data->d_transformV_p->d_markVPtsAreInCurrentConfig =
          config["Compute"][set]["Mark_V_0"]["Points_Current_Config"]
              .as<bool>();

    if (config["Compute"][set]["Mark_V_0"]["Nodes"])
      for (auto &&pt : config["Compute"][set]["Mark_V_0"]["Nodes"])
        data->d_transformV_p->d_markVNodes.push_back(pt.as<size_t>());

    if (data->d_transformV_p->d_markVPts.empty() &&
        data->d_transformV_p->d_markVNodes.empty() &&
        !data->d_transformV_p->d_markVInRectGiven)
      data->d_transformV_p->d_markVAsZero = false;
  }

  // Symmetrize velocity field
  if (config["Compute"][set]["Symmetrize_V"]) {
    if (!data->d_transformV_p)
      data->d_transformV_p = new tools::pp::TransformVelocity();
    data->d_transformV_p->d_symmetrizeV = true;

    // check if we need to combine this with mark_v operation
    if (config["Compute"][set]["Symmetrize_V"]["Combine_Mark_V_0"])
      data->d_transformV_p->d_combineMarkV =
          config["Compute"][set]["Symmetrize_V"]["Combine_Mark_V_0"].as<bool>();

    if (config["Compute"][set]["Symmetrize_V"]["Axis"])
      data->d_transformV_p->d_symmAxis =
          config["Compute"][set]["Symmetrize_V"]["Axis"].as<std::string>();
    else {
      std::cerr << "Error: Need Axis of symmetry for symmetrization of "
                   "velocity.\n";
      exit(1);
    }

    if (config["Compute"][set]["Symmetrize_V"]["Axis_Line"])
      data->d_transformV_p->d_symmLine =
          config["Compute"][set]["Symmetrize_V"]["Axis_Line"].as<double>();
    else {
      std::cerr << "Error: Need location of symmetry axis for symmetrization "
                   "of velocity.\n";
      exit(1);
    }
  }

  // Compute strain and stress
  if (config["Compute"][set]["Strain_Stress"]) {
    if (!data->d_compStrain_p)
      data->d_compStrain_p = new tools::pp::ComputeStrain();

    data->d_compStrain_p->d_computeStrain =
        config["Compute"][set]["Strain_Stress"].as<bool>();
  }

  // Compute magnitude of strain and stress
  if (config["Compute"][set]["Magnitude_Strain_Tensor"]) {
    if (!data->d_compStrain_p)
      data->d_compStrain_p = new tools::pp::ComputeStrain();

    data->d_compStrain_p->d_magStrainTensor = true;
    data->d_compStrain_p->d_computeStrain = true;

    if (config["Compute"][set]["Magnitude_Strain_Tensor"]["Component"])
      data->d_compStrain_p->d_magStrainComp =
          (config["Compute"][set]["Magnitude_Strain_Tensor"]["Component"])
              .as<std::string>();

    if (config["Compute"][set]["Magnitude_Strain_Tensor"]["Cells"])
      for (auto &&i :
           config["Compute"][set]["Magnitude_Strain_Tensor"]["Cells"])
        data->d_compStrain_p->d_markMagStrainCells.emplace_back(
            i[0].as<size_t>(), i[1].as<double>());
  }

  // Crack tip calculation
  if (config["Compute"][set]["Crack_Tip"]) {
    if (!data->d_findCrackTip_p)
      data->d_findCrackTip_p = new tools::pp::FindCrackTip();

    if (config["Compute"][set]["Crack_Tip"]["Same_Dt_Out"])
      data->d_findCrackTip_p->d_crackSameDtOut =
          config["Compute"][set]["Crack_Tip"]["Same_Dt_Out"].as<bool>();
  }

  // J integral
  if (config["Compute"][set]["J_Integral"]) {
    if (!data->d_computeJInt_p)
      data->d_computeJInt_p = new tools::pp::ComputeJIntegral();

    auto e = config["Compute"][set]["J_Integral"];
    if (e["Crack_Orient"])
      data->d_computeJInt_p->d_crackOrient = e["Crack_Orient"].as<int>();
    else {
      std::cerr << "Error: Crack orientation is not provided.\n";
      exit(1);
    }

    if (e["Crack_Tip_File"])
      data->d_computeJInt_p->d_crackTipFile =
          e["Crack_Tip_File"].as<std::string>();
    else {
      std::cerr << "Error: Crack tip information filename is not provided.\n";
      exit(1);
    }

    if (e["Contour_Size"]) {
      for (auto f : e["Contour_Size"])
        data->d_computeJInt_p->d_contourFactor.push_back(f.as<double>());

      if (data->d_computeJInt_p->d_contourFactor.size() == 1)
        data->d_computeJInt_p->d_contourFactor.push_back(
            data->d_computeJInt_p->d_contourFactor[0]);
    } else {
      std::cerr << "Error: Factors to create contour for J integral not "
                   "provided.\n";
      exit(1);
    }
  }
}

void tools::pp::Compute::readCrackTipData(const std::string &filename,
                      std::vector<tools::pp::CrackTipData> *data) {
  // expected format of file:
  // <output step>, <tip x>, <tip y>, <tip vx>, <tip vy>
  io::CSVReader<5> in(filename);
  in.read_header(io::ignore_extra_column, "id", "x", "y", "z", "volume");

  double px, py, vx, vy;
  int n;
  while (in.read_row(n, px, py, vx, vy)) {
    data->emplace_back(tools::pp::CrackTipData(
        size_t(n), util::Point3(px, py, 0.), util::Point3(vx, vy, 0.)));
  }
}

void tools::pp::Compute::initWriter(rw::writer::VtkWriterInterface *writer,
    std::vector<util::Point3> *u) {
  if (writer)
    return;

  writer = new rw::writer::VtkWriterInterface(d_outFilename);
  // append mesh (check if only nodes need to be written)
  if (d_currentData->d_outOnlyNodes)
    writer->appendNodes(d_mesh_p->getNodesP(), u);
  else
    writer->appendMesh(d_mesh_p->getNodesP(), d_mesh_p->getElementType(),
                       d_mesh_p->getElementConnectivitiesP(), u);
}

//
// compute methods
//
void tools::pp::Compute::transformU(rw::writer::VtkWriterInterface *writer) {

  if (!d_currentData->d_transformU_p)
    return;

  std::vector<util::Point3> u_temp(mesh->getNumNodes(), util::Point3());
  auto scale = d_currentData->d_transformU_p->d_scale;
  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_mesh_p->getNumNodes(), [&u_temp, scale, this](boost::uint64_t i) {
        u_temp[i] = util::Point3(scale * d_u[i].d_x, scale * d_u[i].d_y,
                                 scale * d_u[i].d_z);
      });
  f.get();

  if (!writer) {
    initWriter(writer, &u_temp);

    // append original displacement
    writer->appendPointData("Displacement", &d_u);

    // append velocity
    writer->appendPointData("Velocity", &d_v);
  } else {
    std::cerr << "Error: Writer object must be empty when calling scaleU().\n";
    exit(1);
  }
}

void tools::pp::Compute::transformV(rw::writer::VtkWriterInterface *writer) {

  std::vector<util::Point3> v_mark;

  auto data = d_currentData->d_transformV_p;

  // mark v operation (this preceeds symmetrizing of v)
  if (data->d_markVAsZero) {
    v_mark = d_v;

    if (data->d_markVInRectGiven) {
      auto f = hpx::parallel::for_loop(
          hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
          d_mesh_p->getNumNodes(), [&v_mark, data, this](boost::uint64_t i) {
            if (util::geometry::isPointInsideRectangle(
                    d_mesh_p->getNode(i), data->d_markVRect.first.d_x,
                    data->d_markVRect.second.d_x, data->d_markVRect.first.d_y,
                    data->d_markVRect.second.d_y))
              v_mark[i] = util::Point3();
          });
      f.get();
    }

    if (!data->d_markVPts.empty())
      for (auto x : data->d_markVPts) {
        size_t i_found;
        if (data->d_markVPtsAreInCurrentConfig)
          i_found = findNode(x, d_mesh_p->getNodesP(), &d_u);
        else
          i_found = findNode(x, d_mesh_p->getNodesP());

        // modify v_new
        v_mark[i_found] = util::Point3();
      }

    if (!data->d_markVNodes.empty())
      for (auto i : data->d_markVNodes)
        v_mark[i] = util::Point3();

    if (!writer)
      initWriter(writer, &d_u);

    // append velocity
    writer->appendPointData("Mark_Velocity", &v_mark);
  } // mark v operation

  // symmetrize v operation
  if (data->d_symmetrizeV) {

    // use v_mark for modification
    // Reset it to current velocity if we are not combining this
    // with mark_v operation

    if (!data->d_combineMarkV)
      v_mark = d_v;

    auto f = hpx::parallel::for_loop(
        hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
        d_mesh_p->getNumNodes(), [&v_mark, data, this](boost::uint64_t i) {
          auto x = d_mesh_p->getNode(i);

          bool proceed = true;
          if (data->d_symmAxis == "y" && util::compare::definitelyLessThan(
                                             x.d_x, data->d_symmLine + 1.0E-8))
            proceed = false;

          if (data->d_symmAxis == "x" && util::compare::definitelyLessThan(
                                             x.d_y, data->d_symmLine + 1.0E-8))
            proceed = false;

          if (proceed) {
            // find the coordinate of point from where we want to copy
            // the velocity at this node.
            // Mirror image of point
            auto search_x = x;
            if (data->d_symmAxis == "y")
              search_x.d_x = data->d_symmLine - (x.d_x - data->d_symmLine);
            if (data->d_symmAxis == "x")
              search_x.d_y = data->d_symmLine - (x.d_y - data->d_symmLine);

            // search for node at search_x and obtain velocity
            size_t i_found = findNode(search_x, d_mesh_p->getNodesP());

            // write velocity
            v_mark[i] = v_mark[i_found];
            if (data->d_symmAxis == "y")
              v_mark[i].d_x *= -1.;
            if (data->d_symmAxis == "x")
              v_mark[i].d_y *= -1.;
          }
        }); // parallel for loop
    f.get();

    if (!writer)
      initWriter(writer, &d_u);

    // append velocity
    writer->appendPointData("Symm_Velocity", &v_mark);
  }
}

void tools::pp::Compute::computeStrain(rw::writer::VtkWriterInterface *writer) {

  auto data = d_currentData->d_compStrain_p;
  if (!data->d_computeStrain)
    return;

  std::vector<util::SymMatrix3> strain(d_mesh_p->getNumElements(),
                                       util::SymMatrix3());
  std::vector<util::SymMatrix3> stress(d_mesh_p->getNumElements(),
                                       util::SymMatrix3());
  std::vector<float> magS;

  // get Quadrature
  fe::BaseElem *quad;
  if (d_mesh_p->getElementType() == util::vtk_type_triangle)
    quad = new fe::TriElem(1);
  else if (mesh->getElementType() == util::vtk_type_quad)
    quad = new fe::QuadElem(1);
  else {
    std::cerr << "Error: Can not compute strain/stress as the element "
                 "type is not implemented.\n";
    exit(1);
  }

  // compute strain and stress
  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_mesh_p->getNumElements(),
      [&strain, &stress, data, quad, this](boost::uint64_t e) {
        auto ssn = util::SymMatrix3();
        auto sss = util::SymMatrix3();

        // get ids of nodes of element, coordinate of nodes, 1st order
        // quad data, and first quad data
        auto id_nds = d_mesh_p->getElementConnectivity(e);
        auto nds = d_mesh_p->getElementConnectivityNodes(e);
        auto qds = quad->getQuadDatas(nds);
        auto qd0 = qds[0];

        // compute strain in xy plane
        for (size_t i = 0; i < id_nds.size(); i++) {
          auto id = id_nds[i];
          auto ui = d_u[id];

          ssn.d_xx += ui.d_x * qd0.d_derShapes[i][0];
          ssn.d_yy += ui.d_y * qd0.d_derShapes[i][1];
          ssn.d_xy += 0.5 * ui.d_x * qd0.d_derShapes[i][1];
          ssn.d_xy += 0.5 * ui.d_y * qd0.d_derShapes[i][0];
        }

        if (d_matDeck_p->d_isPlaneStrain)
          ssn.d_zz = -d_matDeck_p->d_matData.d_nu * (ssn.d_xx + ssn.d_yy) /
                     (1. - d_matDeck_p->d_matData.d_nu);

        // compute stress
        auto trace = ssn.d_xx + ssn.d_yy + ssn.d_zz;
        sss.d_xx = d_matDeck_p->d_matData.d_lambda * trace * 1. +
                   2. * d_matDeck_p->d_matData.d_mu * ssn.d_xx;
        sss.d_xy = 2. * d_matDeck_p->d_matData.d_mu * ssn.d_xy;
        sss.d_xz = 2. * d_matDeck_p->d_matData.d_mu * ssn.d_xz;

        sss.d_yy = d_matDeck_p->d_matData.d_lambda * trace * 1. +
                   2. * d_matDeck_p->d_matData.d_mu * ssn.d_yy;
        sss.d_yz = 2. * d_matDeck_p->d_matData.d_mu * ssn.d_yz;
        if (!d_matDeck_p->d_isPlaneStrain)
          sss.d_zz = d_matDeck_p->d_matData.d_nu * (sss.d_xx + sss.d_yy);

        strain[e] = ssn;
        stress[e] = sss;
      }); // parallel loop over elements
  f.get();

  // compute magnitude of strain
  if (data->d_magStrainTensor) {
    magS = std::vector<float>(strain.size(), 0.);
    auto f2 = hpx::parallel::for_loop(
        hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
        d_mesh_p->getNumElements(), [&magS, strain, data](boost::uint64_t e) {
          if (data->d_magStrainComp.empty()) {
            magS[e] = std::abs(strain[e].d_xx);
            magS[e] += std::abs(strain[e].d_yy);
            magS[e] += std::abs(strain[e].d_zz);
            magS[e] += std::abs(strain[e].d_xy);
            magS[e] += std::abs(strain[e].d_xz);
            magS[e] += std::abs(strain[e].d_yz);
          } else if (data->d_magStrainComp == "xx") {
            magS[e] = std::abs(strain[e].d_xx);
          } else if (data->d_magStrainComp == "yy") {
            magS[e] = std::abs(strain[e].d_yy);
          }
        });
    f2.get();
  }

  // output strain/stress data
  if (!d_currentData->d_outOnlyNodes) {

    // append mesh
    if (!writer)
      initWriter(writer, &d_u);

    writer->appendCellData("Strain_Tensor", &strain);
    writer->appendCellData("Stress_Tensor", &stress);
    if (data->d_magStrainTensor)
      writer->appendCellData("Mag_Strain", &magS);

    // mark magnitude of strain if asked
    if (!data->d_markMagStrainCells.empty() && data->d_magStrainTensor) {
      for (auto cell : data->d_markMagStrainCells)
        magS[cell.first] = cell.second;

      writer->appendCellData("Mark_Mag_Strain", &magS);
    }
  } else {

    // compute 1st order quad points and store them (quad points in
    // current configuration)
    std::vector<util::Point3> elem_quads =
        std::vector<util::Point3>(d_mesh_p->getNumElements(), util::Point3());

    auto f2 = hpx::parallel::for_loop(
        hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
        d_mesh_p->getNumElements(),
        [&elem_quads, quad, this](boost::uint64_t e) {
          auto nds = d_mesh_p->getElementConnectivity(e);
          std::vector<util::Point3> nds_current(nds.size(),
                                                util::Point3());
          for (size_t j = 0; j < nds.size(); j++)
            nds_current[j] = d_mesh_p->getNode(nds[j]) + d_u[nds[j]];

          auto qds = quad->getQuadPoints(nds_current);
          // store first quad point
          elem_quads[e] = qds[0].d_p;
        });
    f2.get();

    // create unstructured vtk output
    std::string fname = d_outPreTag + "_" + d_currentData->d_tagFilename +
                        "_quads_" + std::to_string(d_nOut);
    auto writer1 = rw::writer::VtkWriterInterface(fname);
    writer1.appendNodes(&elem_quads);
    writer1.appendPointData("Strain_Tensor", &strain);
    writer1.appendPointData("Stress_Tensor", &stress);
    if (data->d_magStrainTensor)
      writer1.appendPointData("Mag_Strain", &magS);

    // mark magnitude of strain if asked
    if (!data->d_markMagStrainCells.empty() && data->d_magStrainTensor) {
      for (auto cell : data->d_markMagStrainCells)
        magS[cell.first] = cell.second;

      writer1.appendPointData("Mark_Mag_Strain", &magS);
    }
    writer1.close();
  }
}

void tools::pp::Compute::computeDamage(std::vector<float> *Z,
                                       rw::writer::VtkWriterInterface *writer,
                                       bool perf_out) {
  if (Z->empty())
    Z->resize(d_mesh_p->getNumNodes());

  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_mesh_p->getNumNodes(),
      [Z, this](boost::uint64_t i) {
        auto xi = d_mesh_p->getNode(i);

        float locz = 0.;
        for (size_t j = 0; j < d_mesh_p->getNumNodes(); j++) {
          if (util::compare::definitelyGreaterThan(xi.dist(d_mesh_p->getNode(j)),
                                                   d_modelDeck_p->d_horizon) ||
              j == i)
            continue;

          auto xj = d_mesh_p->getNode(j);
          if (util::compare::definitelyGreaterThan(xj.dist(xi), 1.0E-10)) {
            auto Sr = std::abs(d_material_p->getS(xj - xi, d_u[j] - d_u[i])) /
                d_material_p->getSc(xj.dist(xi));

            if (util::compare::definitelyLessThan(locz, Sr))
              locz = Sr;
          }
        } // loop over neighbors

        (*Z)[i] = locz;
      }); // parallel loop over nodes
  f.get();

  if (!perf_out)
    return;

  if (!writer)
    initWriter(writer, &d_u);

  writer->appendPointData("Damage_Z", Z);
}

void tools::pp::Compute::findCrackTip(std::vector<float> *Z,
                                      rw::writer::VtkWriterInterface *writer) {

  auto data = d_currentData->d_findCrackTip_p;
  if (d_fractureSet[d_nC] == nullptr) {
    // modify fracture deck for crack tip calculation
    auto fract_deck = *d_fractureDeck_p;
    if (!data->d_crackSameDtOut) {
      fract_deck.d_dtCrackOut = d_outputDeck_p->d_dtOut;
      fract_deck.d_dtCrackVelocity = 1;
    } else {
      fract_deck.d_dtCrackOut = d_outputDeck_p->d_dtOut;
      fract_deck.d_dtCrackVelocity = d_outputDeck_p->d_dtOut;
    }
    fract_deck.d_crackOutFilename = d_currentData->d_tagFilename;

    // get fracture class
    d_fractureSet[d_nC] = new geometry::Fracture(&fract_deck);
  }

  auto n = d_nOut * d_outputDeck_p->d_dtOut;
  auto time = n * d_modelDeck_p->d_dt;

  // compute displacement, damage, and crack tip at k-1 step if
  // crack update is not same as simulation output interval
  if (!data->d_crackSameDtOut) {
    n -= 1;
    time -= d_modelDeck_p->d_dt;
    auto f = hpx::parallel::for_loop(
        hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
        d_mesh_p->getNumNodes(), [this](boost::uint64_t i) {
          d_u[i] -= util::Point3(d_modelDeck_p->d_dt * d_v[i].d_x,
                                 d_modelDeck_p->d_dt * d_v[i].d_y,
                                 d_modelDeck_p->d_dt * d_v[i].d_z);
        });
    f.get();

    // compute damage at n-1
    if (Z->empty())
      Z->resize(d_mesh_p->getNumNodes());
    computeDamage(Z, writer, false);

    // compute crack tip location and crack tip velocity
    d_fractureSet[d_nC]->updateCrackAndOutput(n, time, d_outPath,
                                       d_modelDeck_p->d_horizon,
                                       d_mesh_p->getNodesP(), &d_u, Z);
  }

  // get current displacement
  n = d_nOut * d_outputDeck_p->d_dtOut;
  time = n * d_modelDeck_p->d_dt;
  auto f2 = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_mesh_p->getNumNodes(), [this](boost::uint64_t i) {
        d_u[i] += util::Point3(d_modelDeck_p->d_dt * d_v[i].d_x,
                               d_modelDeck_p->d_dt * d_v[i].d_y,
                               d_modelDeck_p->d_dt * d_v[i].d_z);
      });
  f2.get();

  // compute damage at current displacement
  if (Z->empty())
    Z->resize(d_mesh_p->getNumNodes());
  computeDamage(Z, writer, false);
  d_fractureSet[d_nC]->updateCrackAndOutput(n, time, d_outPath,
                                            d_modelDeck_p->d_horizon,
                                            d_mesh_p->getNodesP(), &d_u, Z);
}

void tools::pp::Compute::computeJIntegral(
    rw::writer::VtkWriterInterface *writer) {

  auto data = d_currentData->d_computeJInt_p;
  if (d_nOut < data->d_start || d_nOut > data->d_end)
    return;

  double energy = 0.;

  auto ctip = data->d_crackTipData[d_nOut - data->d_start];

  // Schematic for horizontal crack (similar for vertical crack)
  //
  //                         D                    C
  //                         + + + + + + + + + + +
  //                         +                   +
  //                         +                   +
  //       ------------------+                   +
  //                         +                   +
  //                         +                   +
  //                         + + + + + + + + + + +
  //                        A                    B
  //
  // Contour is formed by lines A-B, B-C, C-D, D-A
  std::pair<util::Point3, util::Point3> cd(
      std::make_pair(util::Point3(), util::Point3()));
  if (data->d_crackOrient == -1) {
    // vertical crack
    cd.first = util::Point3(ctip.d_p.d_x - 0.5 * data->d_contourFactor[0] *
                                               d_modelDeck_p->d_horizon,
                            ctip.d_p.d_y, 0.);
    cd.second = util::Point3(
        ctip.d_p.d_x +
            0.5 * data->d_contourFactor[0] * d_modelDeck_p->d_horizon,
        ctip.d_p.d_y + data->d_contourFactor[1] * d_modelDeck_p->d_horizon, 0.);
  } else if (data->d_crackOrient == 1) {
    // horizontal crack
    cd.first = util::Point3(ctip.d_p.d_x,
                            ctip.d_p.d_y - 0.5 * data->d_contourFactor[1] *
                                d_modelDeck_p->d_horizon,
                            0.);
    cd.second = util::Point3(
        ctip.d_p.d_x + data->d_contourFactor[0] * d_modelDeck_p->d_horizon,
        ctip.d_p.d_y + 0.5 * data->d_contourFactor[1] * d_modelDeck_p->d_horizon,
        0.);
  }

  // compute nodes and elements list for search
  std::vector<size_t> search_nodes;
  std::vector<size_t> search_elems;
  listElemsAndNodesInDomain(
      cd, d_modelDeck_p->d_horizon + 2. * d_mesh_p->getMeshSize(),
      &search_nodes, &search_nodes);

  // dot product with normal to edge A-B
  auto nv_AB = ctip.d_v.dot(util::Point3(0., -1., 0.));
  // dot product with normal to edge C-D
  auto nv_CD = ctip.d_v.dot(util::Point3(0., 1., 0.));
  // dot product with normal to edge D-A
  auto nv_DA = ctip.d_v.dot(util::Point3(-1., 0., 0.));
  // dot product with normal to B-C
  auto nv_BC = ctip.d_v.dot(util::Point3(1., 0., 0.));

  //
  // compute contribution from edge A-B and C-D
  //
  auto h = d_mesh_p->getMeshSize();
  size_t N = (cd.second.d_x - cd.first.d_x) / h;
  if (util::compare::definitelyLessThan(cd.first.d_x + double(N) * h,
                                        cd.second.d_x))
    N++;

  // create second order quadrature class for 1-d line element
  auto line_quad = fe::LineElem(2);

  for (size_t I = 0; I < N; I++) {
    // line element
    auto x1 = cd.first.d_x + double(I) * h;
    auto x2 = cd.first.d_x + double(I + 1) * h;
    if (I == N - 1)
      x2 = cd.second.d_x;

    // get quadrature points
    auto qds = line_quad.getQuadPoints(std::vector<util::Point3>{
        util::Point3(x1, 0., 0.), util::Point3(x2, 0., 0.)});

    // loop over quad points
    for (auto qd : qds) {
      // process edge A-B
      qd.d_p.d_y = cd.first.d_y;

      // get contribution
      energy += getContourContribJInt(qd.d_p, &search_nodes, &search_elems) *
                nv_AB * qd.d_w;

      // process edge C-D
      qd.d_p.d_y = cd.second.d_y;

      // get contribution
      energy += getContourContribJInt(qd.d_p, &search_nodes, &search_elems) *
          nv_CD * qd.d_w;
    }
  }

  //
  // compute contribution from edge B-C and D-A
  //
  N = (cd.second.d_y - cd.first.d_y) / h;
  if (util::compare::definitelyLessThan(cd.first.d_y + double(N) * h,
                                        cd.second.d_y))
    N++;

  for (size_t I = 0; I < N; I++) {
    // line element
    auto y1 = cd.first.d_y + double(I) * h;
    auto y2 = cd.first.d_y + double(I + 1) * h;
    if (I == N - 1)
      y2 = cd.second.d_y;

    // get quadrature points
    auto qds = line_quad.getQuadPoints(std::vector<util::Point3>{
        util::Point3(0., y1, 0.), util::Point3(0., y2, 0.)});

    // loop over quad points
    for (auto qd : qds) {
      // process edge B-C
      qd.d_p.d_x = cd.second.d_x;

      // get contribution
      energy += getContourContribJInt(qd.d_p, &search_nodes, &search_elems) *
          nv_BC * qd.d_w;

      // process edge D-A
      qd.d_p.d_x = cd.first.d_x;

      // get contribution
      energy += getContourContribJInt(qd.d_p, &search_nodes, &search_elems) *
          nv_DA * qd.d_w;
    }
  }

  //
  // Contribution from work done by peridynamic force
  //
  for (auto i : search_nodes) {

    auto xi = d_mesh_p->getNode(i);
  }
}

//
// utility methods
//
void tools::pp::Compute::addUniqueToList(size_t i, std::vector<size_t> *list) {

  for (auto a : *list) {
    if (a == i)
      return;
  }
  list->emplace_back(i);
}

size_t tools::pp::Compute::findNode(const util::Point3 &x,
                                    const std::vector<util::Point3> *nodes,
                                    const std::vector<util::Point3> *u) {
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

void tools::pp::Compute::listElemsAndNodesInDomain(
    const std::pair<util::Point3, util::Point3> &cd, const double &tol,
    std::vector<size_t> *nodes, std::vector<size_t> *elements) {

  // nodes list
  nodes->clear();
  for (size_t i = 0; i < d_mesh_p->getNumNodes(); i++) {

    auto x = d_mesh_p->getNode(i);

    // check if node is in the bigger domain and not in smaller domain
    if (util::geometry::isPointInsideRectangle(
            x, cd.first.d_x - tol, cd.second.d_x + tol, cd.first.d_y - tol,
            cd.second.d_y + tol) &&
        !util::geometry::isPointInsideRectangle(
            x, cd.first.d_x + tol, cd.second.d_x - tol, cd.first.d_y + tol,
            cd.second.d_y - tol))
      addUniqueToList(i, nodes);
  }

  // element list
  elements->clear();
  for (size_t e = 0; e < d_mesh_p->getNumElements(); e++) {
    auto ids = d_mesh_p->getElementConnectivity(e);
    bool add_e = false;
    // one of the node has to be inside the domain and at least one node
    // has to be outside
    for (auto i : ids) {
      if (add_e)
        continue;
      // check if this node is inside the domain
      if (util::geometry::isPointInsideRectangle(d_mesh_p->getNode(i),
                                                 cd.first.d_x, cd.second.d_x,
                                                 cd.first.d_y, cd.second.d_y)) {

        // now check if there is at least one node outside the domain
        for (auto j : ids) {
          if (j == i)
            continue;
          if (!util::geometry::isPointInsideRectangle(
                  d_mesh_p->getNode(j), cd.first.d_x, cd.second.d_x,
                  cd.first.d_y, cd.second.d_y)) {
            add_e = true;
          }
        }
      }
    } // loop over nodes of element e

    // add e to list
    if (add_e)
      addUniqueToList(e, elements);
  }
}

bool tools::pp::Compute::triCheckAndInterpolateUV(
    const util::Point3 &p, util::Point3 &up, util::Point3 &vp,
    const std::vector<size_t> &ids) {

  // get triangle element object
  auto tri_quad = fe::TriElem(0);
  auto nodes = std::vector<util::Point3>{
      d_mesh_p->getNode(ids[0]), d_mesh_p->getNode(ids[1]), d_mesh_p->getNode(ids[2])};
  double area = std::abs(tri_quad.elemSize(nodes));
  double sum_area = 0.;
  for (size_t a = 0; a < nodes.size(); a++)
    sum_area += std::abs(tri_quad.elemSize(std::vector<util::Point3>{
        p, nodes[a], nodes[(a == nodes.size() - 1 ? 0 : a + 1)]}));

  if (util::compare::definitelyGreaterThan(sum_area, area))
    return false;

  // point p is inside triangle and so we compute shape functions at this
  // point and perform interpolation
  auto shapes = tri_quad.getShapes(p, nodes);

  // interpolate
  up = util::Point3();
  vp = util::Point3();
  for (size_t i=0; i<3; i++) {
    up.d_x += shapes[i] * d_u[ids[i]].d_x;
    up.d_y += shapes[i] * d_u[ids[i]].d_y;

    vp.d_x += shapes[i] * d_v[ids[i]].d_x;
    vp.d_y += shapes[i] * d_v[ids[i]].d_y;
  }

  return true;
}

void tools::pp::Compute::interpolateUV(const util::Point3 &p, util::Point3 &up,
                                       util::Point3 &vp,
                                       const std::vector<size_t> *nodes,
                                       const std::vector<size_t> *elements) {
  // check if element data is available
  if (elements->empty()) {
    // use piecewise constant interpolation
    double dist = 1000.;
    long int loc_i = -1;
    // search for closest node
    for (auto i : *nodes) {
      if (util::compare::definitelyLessThan(p.dist(d_mesh_p->getNode(i)), dist)) {
        dist = p.dist(d_mesh_p->getNode(i));
        loc_i = i;
      }
    }
    up = d_u[loc_i];
    vp = d_v[loc_i];

    return;
  }
  else {
    for (auto e : *elements) {
      auto ids = d_mesh_p->getElementConnectivity(e);

      // cases
      if (d_mesh_p->getElementType() == util::vtk_type_triangle) {
        if (triCheckAndInterpolateUV(p, up, vp, ids))
          return;
      } else if (d_mesh_p->getElementType() == util::vtk_type_quad) {
        // check in triangle {v1, v2, v3}
        if (triCheckAndInterpolateUV(
            p, up, vp, std::vector<size_t>{ids[0], ids[1], ids[2]}))
          return;

        // check in triangle {v1, v3, v4}
        if (triCheckAndInterpolateUV(
            p, up, vp, std::vector<size_t>{ids[0], ids[2], ids[3]}))
          return;
      }
    }

    // issue error since element containing the point is not found
    std::cerr << "Error: Can not find element for point p = (" << p.d_x
              << ", " << p.d_y << ") for interpolation.\n";
    exit(1);
  }
}

double tools::pp::Compute::getContourContribJInt(
    const util::Point3 &p, const std::vector<size_t> *nodes,
    const std::vector<size_t> *elements) {

  // get displacement and velocity at point
  auto uq = util::Point3();
  auto vq = util::Point3();
  interpolateUV(p, uq, vq, nodes, elements);

  // add kinetic energy
  double energy = 0.5 * d_material_p->getDensity() * vq.dist(vq);

  // add peridynamic energy
  for (auto j : *nodes) {
    auto xj = d_mesh_p->getNode(j);
    if (util::compare::definitelyGreaterThan(xj.dist(p),
                                             d_modelDeck_p->d_horizon) ||
        util::compare::definitelyLessThan(xj.dist(p), 1.0E-10))
      continue;

    // get displacement and strain
    auto uj = d_u[j];
    auto rjq = xj.dist(p);
    auto Sjq = d_material_p->getS(xj - p, uj - uq);

    // get volume correction
    auto volj = d_mesh_p->getNodalVolume(j);
    if (util::compare::definitelyGreaterThan(
            rjq, d_modelDeck_p->d_horizon - 0.5 * d_mesh_p->getMeshSize()))
      volj *= (d_modelDeck_p->d_horizon + 0.5 * d_mesh_p->getMeshSize() - rjq) /
              d_mesh_p->getMeshSize();

    bool fracture_state = false;
    auto ef = d_material_p->getBondEF(rjq, Sjq, fracture_state, true);

    // add contribution to energy
    energy += ef.first * volj;
  } // loop over nodes for pd energy density

  return energy;
}
