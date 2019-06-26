// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "compute.h"
#include "external/csv.h"       // csv reader
#include "fe/lineElem.h"        // definition of LineElem
#include "fe/quadElem.h"        // definition of QuadElem
#include "fe/triElem.h"         // definition of TriElem
#include "inp/policy.h"         // definition of Policy
#include "rw/reader.h"          // definition of readVtuFileRestart
#include "util/feElementDefs.h" // definition of fe element type
#include "util/utilGeom.h"      // definition of isPointInsideRectangle
#include <cmath>
#include <hpx/include/parallel_algorithm.hpp>
#include <iostream>
#include <util/fastMethods.h>
#include <yaml-cpp/yaml.h> // YAML reader

bool minSortZ(tools::pp::SortZ a, tools::pp::SortZ b) {
  return util::compare::definitelyLessThan(a.d_Z, b.d_Z);
}

tools::pp::Compute::Compute(const std::string &filename)
    : d_inpFilename(filename), d_nOut(0), d_nC(0), d_currentData(nullptr),
      d_writerReady(false), d_modelDeck_p(nullptr), d_outputDeck_p(nullptr),
      d_fractureDeck_p(nullptr), d_matDeck_p(nullptr), d_mesh_p(nullptr),
      d_fracture_p(nullptr), d_neighbor_p(nullptr), d_input_p(nullptr),
      d_material_p(nullptr) {

  init();

  // loop over output steps
  int N = d_modelDeck_p->d_Nt / d_outputDeck_p->d_dtOut;
  for (d_nOut = 1; d_nOut <= N; d_nOut++) {
    std::cout << "PP_fe2D: Processing output file = " << d_nOut << "\n";

    // mesh filename to read displacement and velocity
    std::string sim_out_filename =
        d_simOutFilename + std::to_string(d_nOut) + ".vtu";

    // see if we need to read the data
    bool read_file = false;
    for (const auto &data : d_computeData)
      read_file = (d_nOut >= data.d_start && d_nOut <= data.d_end);
    if (!read_file)
      continue;

    // get displacement and velocity
    rw::reader::readVtuFileRestart(sim_out_filename, &d_u, &d_v,
      d_mesh_p->getNodesP());

    // loop over compute sets and do as instructed in input file
    for (d_nC = 0; d_nC < d_computeData.size(); d_nC++) {
      d_currentData = &(d_computeData[d_nC]);
      // check if in the specified bound
      if (d_nOut < d_currentData->d_start || d_nOut > d_currentData->d_end)
        continue;

      // skip d_interval number of simulation file after processing 1 file
      if (d_nOut % d_currentData->d_interval != 0)
        continue;

      std::cout << "  PP_fe2D: Processing compute set = " << d_nC + 1 << "\n";

      // filename for writing postprocessed data
      d_outFilename = d_outPreTag + d_currentData->d_tagFilename + "_" +
                      std::to_string(d_nOut);

      // writer
      rw::writer::VtkWriterInterface writer;
      d_writerReady = false;

      // apply postprocessing
      transformU(&writer);
      transformV(&writer);
      computeStrain(&writer);

      // nodal damage
      std::vector<double> Z;
      findCrackTip(&Z, &writer);
      if (d_currentData->d_damageAtNodes)
        computeDamage(&writer, &Z, true);
      if (!Z.empty())
        Z.shrink_to_fit();

      computeJIntegral();

      // close file
      if (d_writerReady)
        writer.close();
    } // loop compute set
  }   // loop simulation output

  finalize();
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

//  if (config["Output"]["Filename"])
//    d_outPreTag =
//        d_outPath + "/" + config["Output"]["Filename"].as<std::string>();
//  else
    d_outPreTag = d_outPath + "/";

  // create mesh
  std::cout << "PP_fe2D: Creating mesh.\n";
  d_mesh_p = new fe::Mesh(d_input_p->getMeshDeck());

  // material deck and material
  d_material_p = new material::pd::Material(d_input_p->getMaterialDeck(),
                                            d_modelDeck_p->d_dim,
                                            d_modelDeck_p->d_horizon);
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

  // copy edge crack set for findCrackTip() function
  for (auto &d : d_computeData) {
    if (d.d_findCrackTip_p) {
      d.d_findCrackTip_p->d_cracks = d_fractureDeck_p->d_cracks;

      // check which point, top/bottom, we need to track
      auto bbox = d_mesh_p->getBoundingBox();
      for (auto &ck : d.d_findCrackTip_p->d_cracks) {

        // point should be in smaller box inside bounding box
        ck.d_trackt = util::geometry::isPointInsideRectangle(
            ck.d_pt, bbox.first[0] + d_modelDeck_p->d_horizon,
            bbox.second[0] - d_modelDeck_p->d_horizon,
            bbox.first[1] + d_modelDeck_p->d_horizon,
            bbox.second[1] - d_modelDeck_p->d_horizon);

        ck.d_trackb = util::geometry::isPointInsideRectangle(
            ck.d_pb, bbox.first[0] + d_modelDeck_p->d_horizon,
            bbox.second[0] - d_modelDeck_p->d_horizon,
            bbox.first[1] + d_modelDeck_p->d_horizon,
            bbox.second[1] - d_modelDeck_p->d_horizon);
      } // modify crack track status
    }
  }

  // read crack tip data for J integral calculation
  for (auto &d : d_computeData) {
    auto data = d.d_computeJInt_p;
    if (data) {
      readCrackTipData(data->d_crackTipFile, data->d_crackId,
                       &(data->d_crackTipData));
      if (!data->d_crackTipData.empty()) {
        d.d_start = data->d_crackTipData[0].d_n;
        d.d_end = data->d_crackTipData[data->d_crackTipData.size() - 1].d_n;
        data->d_start = data->d_crackTipData[0].d_n;
        data->d_end = data->d_crackTipData[data->d_crackTipData.size() - 1].d_n;
      }
    }
  }
}

void tools::pp::Compute::finalize() {

  // close open file
  for (const auto &c : d_computeData) {
    if (c.d_computeJInt_p) {
      if (c.d_computeJInt_p->d_file)
        fclose(c.d_computeJInt_p->d_file);
    }

    if (c.d_findCrackTip_p) {
      if (c.d_findCrackTip_p->d_filet)
        fclose(c.d_findCrackTip_p->d_filet);
      if (c.d_findCrackTip_p->d_fileb)
        fclose(c.d_findCrackTip_p->d_fileb);
    }
  }
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

  // compress type
  if (config["Compute"][set]["Compress_Type"])
    data->d_compressType =
      config["Compute"][set]["Compress_Type"].as<std::string>();

  // check if start and end time step are specified
  if (config["Compute"][set]["Dt_Start"])
    data->d_start = config["Compute"][set]["Dt_Start"].as<int>();
  if (config["Compute"][set]["Dt_End"])
    data->d_end = config["Compute"][set]["Dt_End"].as<int>();
  if (config["Compute"][set]["Dt_Interval"])
    data->d_interval = config["Compute"][set]["Dt_Interval"].as<int>();
  if (data->d_interval < 1) {
    std::cerr << "Error: Specify valid number (greater than or equal to 1) in"
                 " Dt_Interval of compute " << set << ".\n";
    exit(1);
  }

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

    if (e["Crack_Id"])
      data->d_computeJInt_p->d_crackId = e["Crack_Orient"].as<int>();

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

    if (e["Set_V_Lateral_Comp_Zero"])
      data->d_computeJInt_p->d_setLateralCompVZero =
          e["Set_V_Lateral_Comp_Zero"].as<bool>();
  }
}

void tools::pp::Compute::readCrackTipData(
    const std::string &filename, int crack_id,
    std::vector<tools::pp::CrackTipData> *data) {

  // expected format of file:
  // <crack id>, <output step>, <tip x>, <tip y>, <tip vx>, <tip vy>
  io::CSVReader<6> in(filename);
  in.read_header(io::ignore_extra_column, "'crack_id'", "'dt_out'", "'x'",
                 "'y'", "'vx'", "'vy'");

  double px, py, vx, vy;
  int n, ci;
  while (in.read_row(ci, n, px, py, vx, vy)) {
    if (ci == crack_id)
      data->emplace_back(tools::pp::CrackTipData(n, util::Point3(px, py, 0.),
                                                 util::Point3(vx, vy, 0.)));
  }
}

void tools::pp::Compute::initWriter(rw::writer::VtkWriterInterface *writer,
                                    std::vector<util::Point3> *u) {
  if (d_writerReady)
    return;

  writer->open(d_outFilename, d_currentData->d_compressType);
  // append mesh (check if only nodes need to be written)
  if (d_currentData->d_outOnlyNodes)
    writer->appendNodes(d_mesh_p->getNodesP(), u);
  else
    writer->appendMesh(d_mesh_p->getNodesP(), d_mesh_p->getElementType(),
                       d_mesh_p->getElementConnectivitiesP(), u);

  d_writerReady = true;
}

//
// compute methods
//
void tools::pp::Compute::transformU(rw::writer::VtkWriterInterface *writer) {

  if (!d_currentData->d_transformU_p)
    return;

  std::vector<util::Point3> u_temp(d_mesh_p->getNumNodes(), util::Point3());
  auto scale = d_currentData->d_transformU_p->d_scale;
  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_mesh_p->getNumNodes(), [&u_temp, scale, this](boost::uint64_t i) {
        u_temp[i] = util::Point3(scale * d_u[i].d_x, scale * d_u[i].d_y,
                                 scale * d_u[i].d_z);
      });
  f.get();

  if (!d_writerReady) {
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

  if (!d_currentData->d_transformV_p)
    return;

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

    if (!d_writerReady)
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

    if (!d_writerReady)
      initWriter(writer, &d_u);

    // append velocity
    writer->appendPointData("Symm_Velocity", &v_mark);
  }
}

void tools::pp::Compute::computeStrain(rw::writer::VtkWriterInterface *writer) {

  if (!d_currentData->d_compStrain_p)
    return;

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
  else if (d_mesh_p->getElementType() == util::vtk_type_quad)
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
    if (!d_writerReady)
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
          std::vector<util::Point3> nds_current(nds.size(), util::Point3());
          for (size_t j = 0; j < nds.size(); j++)
            nds_current[j] = d_mesh_p->getNode(nds[j]) + d_u[nds[j]];

          auto qds = quad->getQuadPoints(nds_current);
          // store first quad point
          elem_quads[e] = qds[0].d_p;
        });
    f2.get();

    // create unstructured vtk output
    std::string fname = d_outPreTag + d_currentData->d_tagFilename +
                        "_quads_" + std::to_string(d_nOut);
    auto writer1 =
        rw::writer::VtkWriterInterface(fname, d_currentData->d_compressType);
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

void tools::pp::Compute::computeDamage(rw::writer::VtkWriterInterface *writer,
                                       std::vector<double> *Z, bool perf_out) {
  //  if (!d_currentData->d_damageAtNodes)
  //    return;

  if (Z->size() != d_mesh_p->getNumNodes())
    Z->resize(d_mesh_p->getNumNodes());

  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_mesh_p->getNumNodes(), [&Z, this](boost::uint64_t i) {
        auto xi = d_mesh_p->getNode(i);

        double locz = 0.;
        for (size_t j = 0; j < d_mesh_p->getNumNodes(); j++) {
          if (util::compare::definitelyGreaterThan(
                  xi.dist(d_mesh_p->getNode(j)), d_modelDeck_p->d_horizon) ||
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

  if (!d_writerReady)
    initWriter(writer, &d_u);

  writer->appendPointData("Damage_Z", Z);
}

void tools::pp::Compute::findCrackTip(std::vector<double> *Z,
                                      rw::writer::VtkWriterInterface *writer) {

  auto data = d_currentData->d_findCrackTip_p;
  if (!data)
    return;

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
    if (Z->size() != d_mesh_p->getNumNodes())
      Z->resize(d_mesh_p->getNumNodes());
    computeDamage(writer, Z, false);

    // compute crack tip location and crack tip velocity
    updateCrack(time, Z);
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
  if (Z->size() != d_mesh_p->getNumNodes())
    Z->resize(d_mesh_p->getNumNodes());
  computeDamage(writer, Z, false);

  // compute crack tip location and crack tip velocity
  updateCrack(time, Z);

  // perform output
  crackOutput();
}

void tools::pp::Compute::computeJIntegral() {

  auto data = d_currentData->d_computeJInt_p;
  if (!data)
    return;

  // to hold energy into crack
  double energy = 0.;

  // get crack tip data
  auto ctip = data->d_crackTipData[d_nOut - data->d_start];

  // set lateral component of crack velocity zero
  if (data->d_setLateralCompVZero) {
    if (data->d_crackOrient == -1)
      ctip.d_v.d_x = 0.;
    else if (data->d_crackOrient == 1)
      ctip.d_v.d_y = 0.;
  }

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
  double tolx = 0.5 * data->d_contourFactor[0] * d_modelDeck_p->d_horizon;
  double toly = 0.5 * data->d_contourFactor[1] * d_modelDeck_p->d_horizon;
  auto cd = std::make_pair(
      util::Point3(ctip.d_p.d_x - tolx, ctip.d_p.d_y - toly, 0.),
      util::Point3(ctip.d_p.d_x + tolx, ctip.d_p.d_y + toly, 0.));

  // compute nodes and elements list for search
  std::vector<size_t> search_nodes;
  std::vector<size_t> search_elems;
  listElemsAndNodesInDomain(
      cd, d_modelDeck_p->d_horizon + 2. * d_mesh_p->getMeshSize(),
      &search_nodes, &search_elems);

  //
  // Compute contour integral
  //
  // create second order quadrature class for 1-d line element
  auto line_quad = fe::LineElem(2);
  auto h = d_mesh_p->getMeshSize();
  for (size_t E=0; E<2; E++) {
    long N = 0;
    if (E==0) {
      // number of elements for horizontal edge
      N = (cd.second.d_x - cd.first.d_x) / h;
      if (util::compare::definitelyLessThan(cd.first.d_x + double(N) * h,
                                            cd.second.d_x))
        N++;
    } else {
      // number of elements for vertical edge
      N = (cd.second.d_y - cd.first.d_y) / h;
      if (util::compare::definitelyLessThan(cd.first.d_y + double(N) * h,
                                            cd.second.d_y))
        N++;
    }

    auto energies = std::vector<double>(N, 0.);
    auto f = hpx::parallel::for_loop(
        hpx::parallel::execution::par(hpx::parallel::execution::task), 0, N,
        [&energies, N, h, cd, ctip, &line_quad, search_nodes, search_elems, E,
         this](boost::uint64_t I) {

          double loc_energy = 0.;

          // line element
          auto x1 = 0.;
          auto x2 = 0.;
          if (E==0) {
            // discretization of horizontal line
            x1 = cd.first.d_x + double(I) * h;
            x2 = cd.first.d_x + double(I + 1) * h;
            if (I == N - 1)
              x2 = cd.second.d_x;
          } else {
            // discretization of vertical line
            x1 = cd.first.d_y + double(I) * h;
            x2 = cd.first.d_y + double(I + 1) * h;
            if (I == N - 1)
              x2 = cd.second.d_y;
          }

          // get quadrature points
          auto qds = line_quad.getQuadPoints(std::vector<util::Point3>{
              util::Point3(x1, 0., 0.), util::Point3(x2, 0., 0.)});

          // loop over quad points
          for (auto qd : qds) {
            if (E==0) {
              // process edge A-B
              qd.d_p.d_y = cd.first.d_y;

              // get contribution
              // n dot v for edge A-B = - (y component of v)
              loc_energy +=
                  getContourContribJInt(qd.d_p, &search_nodes, &search_elems) *
                  (-ctip.d_v.d_y) * qd.d_w;

              // process edge C-D
              qd.d_p.d_y = cd.second.d_y;

              // get contribution
              // n dot v for edge C-D = y component of v
              loc_energy +=
                  getContourContribJInt(qd.d_p, &search_nodes, &search_elems) *
                  ctip.d_v.d_y * qd.d_w;
            } else {
              // process edge B-C
              // transform quad point along vertical line to correct coordinate
              auto p_temp = qd.d_p;
              qd.d_p = util::Point3(cd.second.d_x, p_temp.d_x, 0.);

              // get contribution
              // n dot v for edge B-C = (x component of v)
              loc_energy +=
                  getContourContribJInt(qd.d_p, &search_nodes, &search_elems) *
                      ctip.d_v.d_x * qd.d_w;

              // process edge D-A
              // transform quad point along vertical line to correct coordinate
              qd.d_p = util::Point3(cd.first.d_x, p_temp.d_x, 0.);

              // get contribution
              // n dot v for edge D-A = - (x component of v)
              loc_energy +=
                  getContourContribJInt(qd.d_p, &search_nodes, &search_elems) *
                      (-ctip.d_v.d_x) * qd.d_w;
            }
          } // loop over quad points

          energies[I] = loc_energy;
        });
    f.get();

    // add energies
    energy += util::methods::add(energies);
  }

  //
  // Contribution from work done by peridynamic force
  //
  // decompose search_nodes list in two parts: one list for nodes outside
  // domain A formed by contour and other for nodes on contour and inside
  // domain A.
  std::vector<size_t> search_node_comp;
  decomposeSearchNodes(cd, &search_nodes, &search_node_comp);
  auto energies = std::vector<double>(search_node_comp.size(), 0.);

  // loop over nodes in compliment of domain A
  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      search_node_comp.size(),
      [&energies, h, cd, &line_quad, search_nodes, search_node_comp,
       this](boost::uint64_t i) {
        auto id = search_node_comp[i];
        auto xi = d_mesh_p->getNode(id);
        auto ui = d_u[id];
        auto vi = d_v[id];
        auto voli = d_mesh_p->getNodalVolume(id);

        double energy_loc = 0.;

        // loop over nodes in domain A
        for (auto j : search_nodes) {

          auto xj = d_mesh_p->getNode(j);
          auto rji = xj.dist(xi);
          if (util::compare::definitelyGreaterThan(rji,
                                                   d_modelDeck_p->d_horizon) ||
              j == id)
            continue;

          auto v_sum = d_v[j] + vi;
          auto xji = xj - xi;
          auto Sji = d_material_p->getS(xji, d_u[j] - ui);

          // get corrected volume of node j
          auto volj = d_mesh_p->getNodalVolume(j);
          if (util::compare::definitelyGreaterThan(
                  rji, d_modelDeck_p->d_horizon - 0.5 * h))
            volj *= (d_modelDeck_p->d_horizon + 0.5 * h - rji) / h;

          // get bond force
          bool fracture_state = false;
          auto ef = d_material_p->getBondEF(rji, Sji, fracture_state, true);

          // compute the contribution
          energy_loc -= ef.second * volj *
                        (xji.d_x * v_sum.d_x + xji.d_y * v_sum.d_y +
                         xji.d_z * v_sum.d_z) /
                        rji;
        } // loop over neighboring nodes

        energies[i] = energy_loc * voli;
      });
  f.get();

  // add energies
  energy += util::methods::add(energies);

  // create file in first call
  if (!data->d_file) {
    std::string filename =
        d_outPreTag + d_currentData->d_tagFilename + ".csv";
    data->d_file = fopen(filename.c_str(), "w");

    // write header
    fprintf(data->d_file, "dt_out, vmag, E, v*Gc, Gc_experiment, Gc_theory\n");
  }

  // write data
  double Gcompute = 0.;
  auto vmag = ctip.d_v.length();
  if (util::compare::definitelyGreaterThan(vmag, 1.E-10))
    Gcompute = energy / vmag;
  fprintf(data->d_file, "%u, %4.6e, %4.6e, %4.6e, %4.6e, %4.6e\n", d_nOut, vmag,
          energy, vmag * d_matDeck_p->d_matData.d_Gc, Gcompute,
          d_matDeck_p->d_matData.d_Gc);
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

    // special care is required for corner points
    if (!add_e) {
      // search element which contains any of the corner point
      for (size_t corner = 0; corner < 4; corner++) {

        if (add_e)
          continue;

        auto p = util::Point3();
        if (corner == 0) {
          p = cd.first;
        } else if (corner == 1) {
          p = cd.second;
        } else if (corner == 2) {
          p.d_x = cd.first.d_x;
          p.d_y = cd.second.d_y;
        } else if (corner == 3) {
          p.d_x = cd.second.d_x;
          p.d_y = cd.first.d_y;
        }

        // search
        // cases
        auto up = util::Point3();
        auto vp = util::Point3();
        if (d_mesh_p->getElementType() == util::vtk_type_triangle) {
          add_e = triCheckAndInterpolateUV(p, up, vp, ids, true);
        } else if (d_mesh_p->getElementType() == util::vtk_type_quad) {
          // check in triangle {v1, v2, v3}
          add_e = triCheckAndInterpolateUV(
              p, up, vp, std::vector<size_t>{ids[0], ids[1], ids[2]}, true);

          // check in triangle {v1, v3, v4}
          add_e = triCheckAndInterpolateUV(
              p, up, vp, std::vector<size_t>{ids[0], ids[2], ids[3]}, true);
        }
      } // loop over corners
    }

    // add e to list
    if (add_e)
      addUniqueToList(e, elements);
  }
}

void tools::pp::Compute::decomposeSearchNodes(
    const std::pair<util::Point3, util::Point3> &cd, std::vector<size_t> *nodes,
    std::vector<size_t> *nodes_new) {

  std::vector<size_t> nodes_temp = *nodes;
  nodes_new->clear();
  nodes->clear();
  for (auto i : nodes_temp) {
    auto xi = d_mesh_p->getNode(i);
    if (util::geometry::isPointInsideRectangle(d_mesh_p->getNode(i),
                                               cd.first.d_x, cd.second.d_x,
                                               cd.first.d_y, cd.second.d_y))
      nodes->emplace_back(i);
    else
      nodes_new->emplace_back(i);
  }
}

bool tools::pp::Compute::triCheckAndInterpolateUV(
    const util::Point3 &p, util::Point3 &up, util::Point3 &vp,
    const std::vector<size_t> &ids, bool check_only) {

  // get triangle element object
  auto tri_quad = fe::TriElem(0);
  auto nodes = std::vector<util::Point3>{d_mesh_p->getNode(ids[0]),
                                         d_mesh_p->getNode(ids[1]),
                                         d_mesh_p->getNode(ids[2])};
  double area = std::abs(tri_quad.elemSize(nodes));
  double sum_area = 0.;
  for (size_t a = 0; a < nodes.size(); a++)
    sum_area += std::abs(tri_quad.elemSize(std::vector<util::Point3>{
        p, nodes[a], nodes[(a == nodes.size() - 1 ? 0 : a + 1)]}));

  if (util::compare::definitelyGreaterThan(sum_area, area))
    return false;

  if (check_only)
    return true;

  // point p is inside triangle and so we compute shape functions at this
  // point and perform interpolation
  auto shapes = tri_quad.getShapes(p, nodes);

  // interpolate
  up = util::Point3();
  vp = util::Point3();
  for (size_t i = 0; i < 3; i++) {
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
      if (util::compare::definitelyLessThan(p.dist(d_mesh_p->getNode(i)),
                                            dist)) {
        dist = p.dist(d_mesh_p->getNode(i));
        loc_i = i;
      }
    }
    if (loc_i == -1) {
      std::cerr << "Error: Can not find node closer to point p = (" << p.d_x
                << ", " << p.d_y << ").\n";
      exit(1);
    }
    up = d_u[loc_i];
    vp = d_v[loc_i];

    return;
  } else {
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
    std::cerr << "Error: Can not find element for point p = (" << p.d_x << ", "
              << p.d_y << ") for interpolation.\n";

    // for debug
    std::vector<util::Point3> enodes;
    for (auto e : *elements) {
      for (auto n : d_mesh_p->getElementConnectivity(e))
        enodes.emplace_back(d_mesh_p->getNode(n));
    }
    enodes.emplace_back(p);
    std::vector<size_t> etags(enodes.size(), 1);
    etags[etags.size() - 1] = 100;
    auto writer1 = rw::writer::VtkWriterInterface(d_outPreTag +
        d_currentData->d_tagFilename + "_debug_nodes_" + std::to_string(d_nOut));
    writer1.appendNodes(&enodes);
    writer1.appendPointData("Tag", &etags);
    writer1.close();
    exit(1);
  }
}

double
tools::pp::Compute::getContourContribJInt(const util::Point3 &p,
                                          const std::vector<size_t> *nodes,
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

void tools::pp::Compute::updateCrack(const double &time,
                                     const std::vector<double> *Z) {
  auto compute_data = d_currentData->d_findCrackTip_p;

  // loop over crack lines
  for (auto &crack : compute_data->d_cracks) {

    if (crack.d_o == 0)
      continue;

    auto pb = crack.d_pb;
    auto pt = crack.d_pt;

    // find maximum of damage
    double max_Z = util::methods::max(*Z);

    if (util::compare::definitelyLessThan(max_Z, compute_data->d_minZAllowed)) {
      addNewCrackTip(crack, pt, time, true);
      addNewCrackTip(crack, pb, time, false);
      return;
    }

    //
    // Step 1
    //
    // create vector of rectangle domain and data for nodes in each
    // rectangle and its damage
    std::vector<std::pair<util::Point3, util::Point3>> rects_t;
    std::vector<std::pair<util::Point3, util::Point3>> rects_b;
    std::vector<std::vector<size_t>> nodes_t;
    std::vector<std::vector<size_t>> nodes_b;
    std::vector<std::vector<double>> Z_t;
    std::vector<std::vector<double>> Z_b;
    getRectsAndNodesForCrackTip(crack, rects_t, rects_b, nodes_t, nodes_b, Z_t,
                                Z_b, Z);

    // process top point
    if (crack.d_trackt) {
      auto pnew = findTipInRects(crack, max_Z, rects_t, nodes_t, Z_t, Z, true);
      addNewCrackTip(crack, pnew, time, true);
    }

    if (crack.d_trackb) {
      auto pnew = findTipInRects(crack, max_Z, rects_b, nodes_b, Z_b, Z,
          false);
      addNewCrackTip(crack, pnew, time, false);
    }
  } // loop over cracks
}

void tools::pp::Compute::crackOutput() {

  auto compute_data = d_currentData->d_findCrackTip_p;

  // stop when update exceeds upper bound
  int up_bound = 10000;
  if (compute_data->d_updateCount >= up_bound) {
    std::cout << "Warning: Number of times crack data output requested "
                 "exceeds the upper limit 10000.\n";
    if (compute_data->d_filet)
      fclose(compute_data->d_filet);
    if (compute_data->d_fileb)
      fclose(compute_data->d_fileb);
    return;
  }

  // create file in first call
  if (compute_data->d_updateCount == 0) {
    std::string filename =
        d_outPreTag + d_currentData->d_tagFilename + "_t.csv";
    compute_data->d_filet = fopen(filename.c_str(), "w");

    filename = d_outPreTag + d_currentData->d_tagFilename + "_b.csv";
    compute_data->d_fileb = fopen(filename.c_str(), "w");

    // write header
    fprintf(compute_data->d_filet, "'crack_id', 'dt_out', 'x', 'y', 'vx', "
                                   "'vy'\n");
    fprintf(compute_data->d_fileb, "'crack_id', 'dt_out', 'x', 'y', 'vx', "
                                   "'vy'\n");
  }

  // write to file
  for (size_t i = 0; i < compute_data->d_cracks.size(); i++) {
    auto crack = compute_data->d_cracks[i];
    fprintf(compute_data->d_filet, "%lu, %u, %6.8e, %6.8e, %6.8e, %6.8e\n",
            i + 1, d_nOut, crack.d_pt.d_x, crack.d_pt.d_y, crack.d_vt.d_x,
            crack.d_vt.d_y);
    fprintf(compute_data->d_fileb, "%lu, %u, %6.8e, %6.8e, %6.8e, %6.8e\n",
            i + 1, d_nOut, crack.d_pb.d_x, crack.d_pb.d_y, crack.d_vb.d_x,
            crack.d_vb.d_y);
  }

  compute_data->d_updateCount++;
}

void tools::pp::Compute::getRectsAndNodesForCrackTip(
    inp::EdgeCrack &crack,
    std::vector<std::pair<util::Point3, util::Point3>> &rects_t,
    std::vector<std::pair<util::Point3, util::Point3>> &rects_b,
    std::vector<std::vector<size_t>> &nodes_t,
    std::vector<std::vector<size_t>> &nodes_b,
    std::vector<std::vector<double>> &Z_t,
    std::vector<std::vector<double>> &Z_b, const std::vector<double> *Z) {
  auto compute_data = d_currentData->d_findCrackTip_p;

  //
  // Method: For both pt and pb, create sequence of rectangles with
  // increasing distance from pt and pb, and find nodes with suitable
  // damage within each rectangle
  //
  auto pt = crack.d_pt;
  auto pb = crack.d_pb;

  auto h = d_mesh_p->getMeshSize();
  auto horizon = d_modelDeck_p->d_horizon;
  auto bbox = d_mesh_p->getBoundingBox();
  auto bbox_small = bbox;
  bbox_small.first[0] += horizon;
  bbox_small.first[1] += horizon;
  bbox_small.second[0] -= horizon;
  bbox_small.second[1] -= horizon;
  if (!util::geometry::isPointInsideRectangle(
          pt, bbox_small.first[0], bbox_small.second[0], bbox_small.first[1],
          bbox_small.second[1]))
    crack.d_trackt = false;
  if (!util::geometry::isPointInsideRectangle(
          pb, bbox_small.first[0], bbox_small.second[0], bbox_small.first[1],
          bbox_small.second[1]))
    crack.d_trackb = false;

  // size of rectangles for search
  double seq_size = 2. * horizon;
  int Nt = 0;
  int Nb = 0;
  if (crack.d_o == 1) {

    if (crack.d_trackt) {
      Nt = (bbox.second[0] - pt.d_x) / seq_size;
      if (util::compare::definitelyLessThan(Nt * seq_size, bbox.second[0]))
        Nt++;
    }

    if (crack.d_trackb) {
      Nb = (pb.d_x - bbox.first[0]) / seq_size;
      if (util::compare::definitelyLessThan(bbox.first[0] + Nb * seq_size,
                                            pb.d_x))
        Nb++;
    }
  } else if (crack.d_o == -1) {
    if (crack.d_trackt) {
      Nt = (bbox.second[1] - pt.d_y) / seq_size;
      if (util::compare::definitelyLessThan(Nt * seq_size, bbox.second[1]))
        Nt++;
    }

    if (crack.d_trackb) {
      Nb = (pb.d_y - bbox.first[1]) / seq_size;
      if (util::compare::definitelyLessThan(bbox.first[1] + Nb * seq_size,
                                            pb.d_y))
        Nb++;
    }
  }

  // create rectangles
  if (crack.d_trackt) {
    // last rectangle may overstep bounding box but this will not create
    // any problem
    for (int i = 1; i <= Nt; i++) {
      if (crack.d_o == 1)
        rects_t.emplace_back(std::make_pair(
            util::Point3(pt.d_x + (i - 1) * seq_size, pt.d_y - 2. * horizon,
                         0.),
            util::Point3(pt.d_x + i * seq_size, pt.d_y + 2. * horizon, 0.)));
      else if (crack.d_o == -1)
        rects_b.emplace_back(std::make_pair(
            util::Point3(pt.d_x - 2. * horizon, pt.d_y + (i - 1) * seq_size,
                         0.),
            util::Point3(pt.d_x + 2. * horizon, pt.d_y + i * seq_size, 0.)));
    }
  }
  if (crack.d_trackb) {
    // last rectangle may overstep bounding box but this will not create
    // any problem
    for (int i = 1; i <= Nb; i++) {
      if (crack.d_o == 1)
        rects_t.emplace_back(std::make_pair(
            util::Point3(pb.d_x - i * seq_size, pb.d_y - 2. * horizon, 0.),
            util::Point3(pb.d_x - (i - 1) * seq_size, pb.d_y + 2. * horizon,
                         0.)));
      else if (crack.d_o == -1)
        rects_b.emplace_back(std::make_pair(
            util::Point3(pt.d_x - 2. * horizon, pt.d_y - i * seq_size, 0.),
            util::Point3(pt.d_x + 2. * horizon, pt.d_y - (i - 1) * seq_size,
                         0.)));
    }
  }

  // resize nodes list and damage list
  nodes_t.resize(rects_t.size());
  nodes_b.resize(rects_b.size());
  Z_t.resize(rects_t.size());
  Z_b.resize(rects_b.size());

  for (size_t i = 0; i < d_mesh_p->getNumNodes(); i++) {
    auto xi = d_mesh_p->getNode(i);
    auto damage = (*Z)[i];
    if (util::compare::definitelyLessThan(damage,
                                          compute_data->d_minZAllowed) ||
        util::compare::definitelyGreaterThan(damage,
                                             compute_data->d_maxZAllowed))
      continue;

    for (size_t r = 0; r < rects_t.size(); r++) {
      if (util::geometry::isPointInsideRectangle(
              xi, rects_t[r].first.d_x, rects_t[r].second.d_x,
              rects_t[r].first.d_y, rects_t[r].second.d_y)) {
        if (r == 0) {
          nodes_t[r].emplace_back(i);
          Z_t[r].emplace_back(damage);
        } else {
          // see if the node i is in the previous vector
          bool found = false;
          for (auto j : nodes_t[r - 1]) {
            if (j == i) {
              found = true;
              continue;
            }
          }
          if (!found) {
            nodes_t[r].emplace_back(i);
            Z_t[r].emplace_back(damage);
          }
        }
      }
    } // top point
    for (size_t r = 0; r < rects_b.size(); r++) {
      if (util::geometry::isPointInsideRectangle(
              xi, rects_b[r].first.d_x, rects_b[r].second.d_x,
              rects_b[r].first.d_y, rects_b[r].second.d_y)) {
        if (r == 0) {
          nodes_b[r].emplace_back(i);
          Z_b[r].emplace_back(damage);
        } else {
          // see if the node i is in the previous vector
          bool found = false;
          for (auto j : nodes_b[r - 1]) {
            if (j == i) {
              found = true;
              continue;
            }
          }
          if (!found) {
            nodes_b[r].emplace_back(i);
            Z_b[r].emplace_back(damage);
          }
        }
      }
    } // bottom point
  }   // loop over nodes
}

void tools::pp::Compute::addNewCrackTip(inp::EdgeCrack &crack, util::Point3 pnew,
                                      double time, bool is_top) {
  auto compute_data = d_currentData->d_findCrackTip_p;

  if (is_top) {
    crack.d_oldPt = crack.d_pt;
    crack.d_pt = pnew;
    auto diff = crack.d_pt - crack.d_oldPt;
    auto delta_t = time - compute_data->d_timet;
    crack.d_lt += diff.length();
    crack.d_l += diff.length();
    crack.d_vt = util::Point3(diff.d_x / delta_t, diff.d_y / delta_t,
                              diff.d_z / delta_t);
    compute_data->d_timet = time;
  } else {
    crack.d_oldPb = crack.d_pb;
    crack.d_pb = pnew;
    auto diff = crack.d_pb - crack.d_oldPb;
    auto delta_t = time - compute_data->d_timet;
    crack.d_lb += diff.length();
    crack.d_l += diff.length();
    crack.d_vb =
        util::Point3(diff.d_x / delta_t, diff.d_y / delta_t, diff.d_z / delta_t);
    compute_data->d_timeb = time;
  }
}

util::Point3 tools::pp::Compute::findTipInRects(inp::EdgeCrack &crack,
    const double &max_Z,
    const std::vector<std::pair<util::Point3, util::Point3>> &rects,
    const std::vector<std::vector<size_t>> &nodes,
    const std::vector<std::vector<double>> &Zs,
    const std::vector<double> *Z,
    bool is_top) {

  if (is_top)
    std::cout << "    Tip information: Top     \n";
  else
    std::cout << "    Tip information: Bottom     \n";
  std::cout << "    max_Z = " << max_Z << "\n";
  std::cout << "    size rects = " << rects.size() << "\n";
  // find size of rectangle with atleast one node
  size_t rect_size_with_node = 0;
  for (size_t r=0; r<rects.size(); r++) {
    if (nodes[r].size() > 0)
      rect_size_with_node++;
  }
  std::cout << "    size valid rects = " << rect_size_with_node << "\n";

  // get old tip
  auto pold = crack.d_pt;
  if (!is_top)
    pold = crack.d_pb;

  // if there are no rectangles with desired nodes, return old tip
  if (rect_size_with_node == 0)
    return pold;

  //
  // Step 2
  //
  // order rectangles in increasing value of damage
  // for each rectangle find minimum damage and then order rectangles based
  // on the minimum value of damage
  std::vector<tools::pp::SortZ> sortZ;
  for (size_t r = 0; r < rects.size(); r++) {
    auto a = tools::pp::SortZ();
    size_t i = 0;
    a.d_r = r;
    if (Zs[r].size() > 0) {
      a.d_Z = util::methods::min(Zs[r], &i);
      a.d_i = nodes[r][i];
      sortZ.emplace_back(a);
    }
  }
  // sort the rectangles in increasing value of damage
  std::sort(sortZ.begin(), sortZ.end(), minSortZ);

  //
  // Step 3
  //
  //
  // for rectangle sortZ[0].r and sortZ[1].r we see if we can find a node
  // more closer to the crack line
  //
  for (size_t r = 0; r < sortZ.size(); r++) {
    if (r > 2)
      continue;

    auto mz = sortZ[r];
    auto y0 = d_mesh_p->getNode(mz.d_i) + d_u[mz.d_i];
    auto Z0 = mz.d_Z;
    double diff_Z_for_alternate_point = 0.1;
    if (sortZ.size() > r + 1) {
      auto diff_p1 = sortZ[r + 1].d_Z - Z0;
      if (util::compare::definitelyLessThan(diff_p1,
                                            diff_Z_for_alternate_point))
        diff_Z_for_alternate_point = diff_p1;
    }
    double dist_crack_line = 0.;
    if (crack.d_o == -1)
      dist_crack_line = std::abs(y0.d_x - pold.d_x);
    else if (crack.d_o == 1)
      dist_crack_line = std::abs(y0.d_y - pold.d_y);

    long i_new = -1;
    for (auto x : nodes[mz.d_r]) {
      if (x == mz.d_i)
        continue;

      auto y = d_mesh_p->getNode(x) + d_u[x];
      double dist = 0.;
      if (crack.d_o == -1)
        dist = std::abs(y.d_x - pold.d_x);
      else if (crack.d_o == 1)
        dist = std::abs(y.d_y - pold.d_y);

      if (util::compare::definitelyLessThan(dist, dist_crack_line) &&
          util::compare::definitelyLessThan(std::abs(Z0 - (*Z)[x]),
                                            diff_Z_for_alternate_point)) {

        i_new = x;
        dist_crack_line = dist;
      }
    }

    // update sortZ[r]
    if (i_new >= 0) {
      sortZ[r].d_i = i_new;
      sortZ[r].d_Z = (*Z)[i_new];
    }
  }

  //
  // Step 4
  //
  //
  // Our first choice is rectangle sortZ[0].r
  // However, we see if rectangle sortZ[0].r and sortZ[1].r have symmetry
  //
  // Here by symmetry we mean if there are two nodes in opposite sides of
  // crack line with same damage
  //
  std::vector<long> sym_rect(2, -1);

  // Three choices of crack line in decreasing order of preference
  // Choice 1: crack line defined by initial crack tip, i.e. initial crack line
  // Choice 2: crack line defined by old crack tip
  // Choice 3: crack line defined by current crack tip
  auto p_choices = std::vector<util::Point3>{crack.d_initPt, crack.d_oldPt,
                                             crack.d_pt};
  if (!is_top)
    p_choices = std::vector<util::Point3>{crack.d_initPb, crack.d_oldPb,
                                          crack.d_pb};
  for (size_t r = 0; r < sortZ.size(); r++) {

    // only consider first and second rectangle in sorted rectangle list
    if (r > 2)
      continue;

    // get data for this rectangle
    auto mz = sortZ[r];

    // do not consider node with large damage difference to the minimum
    // damage in mz.d_Z
    // diff_z will be further reduced to find node with closest damage to mz.d_Z
    double diff_z = 1.0E-02;

    // loop over crack line choices
    for (size_t f = 0; f < p_choices.size(); f++) {
      auto p_search = p_choices[f];

      // check if we already found symmetrically opposite node for crack line
      // before this crack line
      if (sym_rect[r] >= 0)
        continue;

      for (auto x : nodes[mz.d_r]) {

        auto y = d_mesh_p->getNode(x) + d_u[x];
        auto y0 = d_mesh_p->getNode(mz.d_i) + d_u[mz.d_i];
        if (util::compare::definitelyLessThan(std::abs((*Z)[x] - mz.d_Z),
                                              diff_z) &&
            x != mz.d_i) {

          // find if it is symmetrically opposite to the node i
          if (crack.d_o == -1) {
            // compare x coordinate
            auto d0 = p_search.d_x - y0.d_x;
            auto d = p_search.d_x - y.d_x;

            if (util::compare::definitelyLessThan(d * d0, 0.)) {
              sym_rect[r] = x;
              diff_z = std::abs((*Z)[x] - mz.d_Z);
            }
          } else if (crack.d_o == 1) {
            // compare y coordinate
            auto d0 = p_search.d_y - y0.d_y;
            auto d = p_search.d_y - y.d_y;

            if (util::compare::definitelyLessThan(d * d0, 0.)) {
              sym_rect[r] = x;
              diff_z = std::abs((*Z)[x] - mz.d_Z);
            }
          }
        }
      } // loop over nodes within rectangle
    }
  } // loop over rectangle

  // output data for debugging
  std::cout << "    sortZ (r, i, Z) = {";
  for (size_t i = 0; i < sortZ.size(); i++) {
    std::cout << "(" << sortZ[i].d_r << "," << sortZ[i].d_i << ","
              << sortZ[i].d_Z << ")";
    if (i < sortZ.size() - 1)
      std::cout << ", ";
  }
  std::cout << "}\n";
  std::cout << "    sym_rect = {" << sym_rect[0] << ", " << sym_rect[1]
            << "}\n";

  //
  // Step 5
  //
  // there are three cases
  // 1. rect 0 is symmetric
  // 2. rect 0 is not symmetric and rect 1 is symmetric
  // 3. rect 0 is not symmetric and rect 0 is not symmetric
  auto pnew = pold;
  if (sym_rect[0] >= 0) {
    // crack tip is given by average of symmetrically opposite points
    size_t i1 = sortZ[0].d_i;
    auto y1 = d_mesh_p->getNode(i1) + d_u[i1];
    size_t i2 = sym_rect[0];
    auto y2 = d_mesh_p->getNode(i2) + d_u[i2];
    if (crack.d_o == -1) {
      pnew.d_y = y1.d_y;
      pnew.d_x = 0.5 * (y1.d_x + y2.d_x);
    } else if (crack.d_o == 1) {
      pnew.d_x = y1.d_x;
      pnew.d_y = 0.5 * (y1.d_y + y2.d_y);
    }

    std::cout << "    Candidate 1 info: id = " << sortZ[0].d_i
              << ", Z = " << sortZ[0].d_Z << ", r = " << sortZ[0].d_r
              << ", p = (" << y1.d_x << ", " << y1.d_y << ")\n";
    std::cout << "    Candidate 2 info: id = " << sym_rect[0]
              << ", Z = " << (*Z)[sym_rect[0]] << ", r = " << sortZ[0].d_r
              << ", p = (" << y2.d_x << ", " << y2.d_y << ")\n";
  } // rect 0 is symmetric
  else {
    if (sym_rect[1] >= 0) {
      size_t i1 = sortZ[1].d_i;
      auto y1 = d_mesh_p->getNode(i1) + d_u[i1];
      size_t i2 = sym_rect[1];
      auto y2 = d_mesh_p->getNode(i2) + d_u[i2];
      if (crack.d_o == -1) {
        pnew.d_y = y1.d_y;
        pnew.d_x = 0.5 * (y1.d_x + y2.d_x);
      } else if (crack.d_o == 1) {
        pnew.d_x = y1.d_x;
        pnew.d_y = 0.5 * (y1.d_y + y2.d_y);
      }
      std::cout << "    Candidate 1 info: id = " << sortZ[1].d_i
                << ", Z = " << sortZ[1].d_Z << ", r = " << sortZ[1].d_r
                << ", p = (" << y1.d_x << ", " << y1.d_y << ")\n";
      std::cout << "    Candidate 2 info: id = " << sym_rect[1]
                << ", Z = " << (*Z)[sym_rect[1]] << ", r = " << sortZ[1].d_r
                << ", p = (" << y2.d_x << ", " << y2.d_y << ")\n";
    } // rect 1 is symmetric
    else {
      // use average between best node in rectangle and old tip
      size_t i1 = sortZ[0].d_i;
      auto y1 = d_mesh_p->getNode(i1) + d_u[i1];
      auto y2 = pold;
      if (crack.d_o == -1) {
        pnew.d_y = y1.d_y;
        pnew.d_x = 0.5 * (y1.d_x + y2.d_x);
      } else if (crack.d_o == 1) {
        pnew.d_x = y1.d_x;
        pnew.d_y = 0.5 * (y1.d_y + y2.d_y);
      }

      std::cout << "    Candidate 1 info: id = " << sortZ[0].d_i
                << ", Z = " << sortZ[0].d_Z << ", r = " << sortZ[0].d_r
                << ", p = (" << y1.d_x << ", " << y1.d_y << ")\n";
    } // rect 1 is not symmetric
  }   // rect 0 is not symmetric

  // output data for debugging
  std::cout << "    Old tip = (" << pold.d_x << ", " << pold.d_y << "), "
            << "New tip = (" << pnew.d_x << ", " << pnew.d_y << ").\n";

  return pnew;
}