////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "compute.h"
#include "external/csv.h"            // csv reader
#include "fe/lineElem.h"             // definition of LineElem
#include "fe/mesh.h"                 // definition of Mesh
#include "fe/quadElem.h"             // definition of QuadElem
#include "fe/triElem.h"              // definition of TriElem
#include "geometry/fracture.h"       // definition of Fracture
#include "geometry/neighbor.h"       // definition of Neighbor
#include "inp/decks/materialDeck.h"  // definition of MaterialDeck
#include "inp/decks/modelDeck.h"     // definition of ModelDeck
#include "inp/decks/outputDeck.h"    // definition of OutputDeck
#include "inp/input.h"               // definition of Input
#include "inp/policy.h"              // definition of Policy
#include "material/materials.h"      // definition of Material
#include "rw/reader.h"               // definition of readVtuFileRestart
#include "rw/writer.h"               // definition of WriterInterface
#include "util/fastMethods.h"        // max and min operation
#include "util/feElementDefs.h"      // definition of fe element type
#include "util/utilGeom.h"           // definition of isPointInsideRectangle
#include <cmath>
#include <hpx/include/parallel_algorithm.hpp>
#include <iostream>
#include <yaml-cpp/yaml.h>  // YAML reader
#include <cfloat>

namespace {

std::ostringstream oss;

// struct to send to lambda function in hpx
struct ContourDataLambdaFn {
  std::vector<double> d_contourPdStrainEnergies;
  std::vector<double> d_contourPdStrainEnergiesRate;
  std::vector<double> d_contourKineticEnergiesRate;
  std::vector<double> d_contourElasticInternalWorksRate;

  ContourDataLambdaFn(){};
};

bool minSortZ(tools::pp::SortZ a, tools::pp::SortZ b) {
  return util::compare::definitelyLessThan(a.d_Z, b.d_Z);
}

void processQuadPointForContour(size_t E, int top_side,
                                const std::pair<util::Point3, util::Point3> &cd,
                                const util::Point3 &ctip_v,
                                const util::Point3 &ctip_d, util::Point3 &p,
                                util::Point3 &edge_normal, double &n_dot_n_c,
                                double &n_dot_v_c) {
  if (E == 0) {
    // process horizontal edges
    if (top_side == 0) {
      // process edge A-B
      // translate quadrature point's y-coordinate to y-coordinate
      // of edge A-B
      p.d_y = cd.first.d_y;
      // normal
      edge_normal = util::Point3(0., -1., 0.);
    } else {
      // process edge C-D
      // translate quadrature point's y-coordinate to y-coordinate
      // of edge C-D
      p.d_y = cd.second.d_y;
      // normal
      edge_normal = util::Point3(0., 1., 0.);
    }

    // n dot n_c where n_c is direction along crack
    n_dot_n_c = edge_normal.dot(ctip_d);
    // n dot v_c where v_c is velocity direction
    n_dot_v_c = edge_normal.dot(ctip_v);
  } else {
    // process vertical edges

    // store original quad point
    double py = p.d_x;

    if (top_side == 0) {
      // process edge B-C
      // transform quad point along vertical line to correct
      // coordinate
      p.d_x = cd.second.d_x;
      p.d_y = py;
      // normal
      edge_normal = util::Point3(1., 0., 0.);
    } else {
      // process edge D-A
      // transform quad point along vertical line to correct
      // coordinate
      p.d_x = cd.first.d_x;
      p.d_y = py;
      // normal
      edge_normal = util::Point3(-1., 0., 0.);
    }

    // n dot n_c where n_c is direction along crack
    n_dot_n_c = edge_normal.dot(ctip_d);
    // n dot v_c where v_c is velocity direction
    n_dot_v_c = edge_normal.dot(ctip_v);
  }  // process vertical edges
}

bool doesBondIntersectCrack(const util::Point3 &p1, const util::Point3 &p2,
                            const tools::pp::CrackTipData &ctip,
                            const int &crack_orient) {
  if (crack_orient == -1) {
    // vertical crack (use x coordinate to know the veritcal line)
    double p1_check = p1.d_x - ctip.d_p.d_x;
    double p2_check = p2.d_x - ctip.d_p.d_x;
    return util::compare::definitelyLessThan(p1_check * p2_check, 0.);

  } else if (crack_orient == 1) {
    // horizontal crack (use y coordinate to know the horizontal line)
    double p1_check = p1.d_y - ctip.d_p.d_y;
    double p2_check = p2.d_y - ctip.d_p.d_y;
    return util::compare::definitelyLessThan(p1_check * p2_check, 0.);
  } else {
    std::cerr << "Error: Extend doesBondIntersectCrack() to arbitrarily "
                 "oriented crack.\n";
    exit(1);
  }
}
}  // namespace

tools::pp::Compute::Compute(const std::string &filename)
    : d_inpFilename(filename),
      d_nOut(0),
      d_nC(0),
      d_time(0.),
      d_currentData(nullptr),
      d_dtOutChange(0),
      d_writerReady(false),
      d_uPlus(false),
      d_dtN(0),
      d_dtStart(0),
      d_dtEnd(0),
      d_fnSuccess(true),
      d_fnErrMsg(""),
      d_needDamageZ(false),
      d_needNeighborList(false),
      d_outputDeck_p(nullptr),
      d_fractureDeck_p(nullptr),
      d_matDeck_p(nullptr),
      d_input_p(nullptr),
      d_material_p(nullptr) {
  d_dataManager_p = new data::DataManager();
  init();

  // build neighbor list if we need it
  if (d_needNeighborList) {
    std::cout << "PP: Computing neighbor list\n";
    d_dataManager_p->setNeighborP(new geometry::Neighbor(
        d_dataManager_p->getModelDeckP()->d_horizon, d_input_p->getNeighborDeck(),
        d_dataManager_p->getMeshP()->getNodesP()));
  }

  // take the smaller output interval as dt_out
  size_t dt_out = d_outputDeck_p->d_dtOutCriteria;
  size_t dt_out_old = d_outputDeck_p->d_dtOutOld;
  size_t dt_interval_factor = dt_out_old / dt_out;

  for (d_nOut = 1; d_nOut <= d_dtN; d_nOut++) {
    size_t current_step = d_nOut * dt_out;
    // proceed every dt_out_old interval if the current time step is
    // less than the time step when the change in output will happen
    bool process_iN = true;
    if (d_nOut < d_dtOutChange)
      if (current_step % dt_out_old != 0) process_iN = false;

    if (!process_iN) continue;

    // get correct factor for dt interval
    if (d_nOut >= d_dtOutChange) dt_interval_factor = 1;

    // current time
    d_time = current_step * d_dataManager_p->getModelDeckP()->d_dt;

    std::cout << "PP_fe2D: Processing output file = " << d_nOut << "\n";

    // mesh filename to read displacement and velocity
    std::string sim_out_filename =
        d_simOutFilename + std::to_string(d_nOut) + ".vtu";

    // see if we need to read the data
    bool read_file = false;
    for (const auto &data : d_computeData) {
      // modify read file flag only if it is false
      if (!read_file)
        read_file = (d_nOut >= data.d_start && d_nOut <= data.d_end);
    }

    if (!read_file) continue;

    // get displacement and velocity
    if (d_outputDeck_p->d_outFormat == "vtu")
      rw::reader::readVtuFileRestart(sim_out_filename, &d_u, &d_v,
                                     d_dataManager_p->getMeshP()->getNodesP());
    else if (d_outputDeck_p->d_outFormat == "msh")
      rw::reader::readMshFileRestart(sim_out_filename, &d_u, &d_v,
                                     d_dataManager_p->getMeshP()->getNodesP());

    if (d_uPlus) {
      auto f2 = hpx::parallel::for_loop(
          hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
          d_dataManager_p->getMeshP()->getNumNodes(),
          [this](boost::uint64_t i) {
            d_u[i] += util::Point3(d_dataManager_p->getModelDeckP()->d_dt * d_v[i].d_x,
                                   d_dataManager_p->getModelDeckP()->d_dt * d_v[i].d_y,
                                   d_dataManager_p->getModelDeckP()->d_dt * d_v[i].d_z);
          });
      f2.get();
    }

    bool damage_computed = false;
    // loop over compute sets and do as instructed in input file
    for (d_nC = 0; d_nC < d_computeData.size(); d_nC++) {
      d_currentData = &(d_computeData[d_nC]);
      // check if in the specified bound
      if (d_nOut < d_currentData->d_start || d_nOut > d_currentData->d_end)
        continue;

      // skip d_interval number of simulation file after processing 1 file
      // For now, we have different method to check when compute set includes
      // crack tip calculation
      bool skip_nout = false;
      if (d_currentData->d_computeJInt_p) {
        if ((d_nOut - d_currentData->d_start) %
                (dt_interval_factor * d_currentData->d_interval) !=
            0)
          skip_nout = true;
      } else {
        if (d_nOut % (dt_interval_factor * d_currentData->d_interval) != 0)
          skip_nout = true;
      }

      if (skip_nout) continue;

      std::cout << "  PP_fe2D: Processing compute set = " << d_nC + 1 << "\n";

      // filename for writing postprocessed data
      d_outFilename = d_outPreTag + d_currentData->d_tagFilename + "_" +
                      std::to_string(d_nOut);

      // writer
      rw::writer::Writer writer;
      d_writerReady = false;

      // get damage at nodes if required
      if (d_needDamageZ) {
        if (!damage_computed) {
          if (!d_Z.empty()) d_Z.clear();
          // check if the input file has damage data
          if (rw::reader::vtuHasPointData(sim_out_filename, "Damage_Z"))
            rw::reader::readVtuFilePointData(sim_out_filename, "Damage_Z",
                                             &d_Z);
          else {
            // see if file for damage is provided
            if (!d_fileZ.empty()) {
              std::string fn = d_fileZ + "_" + std::to_string(d_nOut) + ".vtu";
              if (!rw::reader::readVtuFilePointData(fn, d_tagZ, &d_Z)) {
                std::cerr << "Error: Can not read file = " << fn
                          << " or data = " << d_tagZ
                          << " not available in vtu file.\n";
                exit(1);
              }
            } else
              computeDamage(&writer, &d_Z, false);
          }

          damage_computed = true;
        }
      }

      // apply postprocessing
      transformU(&writer);
      transformV(&writer);
      computeStrain(&writer);

      findCrackTip(&d_Z, &writer);

      computeJIntegral();

      if (d_currentData->d_damageAtNodes) {
        if (!d_writerReady) initWriter(&writer, &d_u);

        writer.appendPointData("Damage_Z", &d_Z);
      }

      // close file
      if (d_writerReady) writer.close();
    }  // loop compute set
  }    // loop simulation output

  finalize();
}

void tools::pp::Compute::init() {
  auto config = YAML::LoadFile(d_inpFilename);

  // get simulation filename
  std::string source_path;
  if (config["Source_Path"])
    source_path = config["Source_Path"].as<std::string>();
  else
    source_path = ".";

  if (config["Simulation_Input_File"])
    d_simInpFilename =
        source_path + "/" + config["Simulation_Input_File"].as<std::string>();
  else {
    std::cerr << "Error: Simulation input file is not provided.\n";
    exit(1);
  }

  // read input data
  std::cout << "PP_fe2D: Reading simulation input file.\n";
  d_input_p = new inp::Input(d_simInpFilename);
  d_dataManager_p->setModelDeckP(d_input_p->getModelDeck());
  d_outputDeck_p = d_input_p->getOutputDeck();
  d_fractureDeck_p = d_input_p->getFractureDeck();

  // to get path for simulation output files, append output path present in
  // the simulation input file
  d_simOutFilename = source_path + "/" + d_outputDeck_p->d_path;
  if (config["Filename_To_Read"])
    d_simOutFilename += config["Filename_To_Read"].as<std::string>() + "_";
  else
    d_simOutFilename += "output_";

  // get output path directory
  if (config["Output"]["Path"])
    d_outPath = config["Output"]["Path"].as<std::string>();
  else
    d_outPath = "./";  // default

  //  if (config["Output"]["Filename"])
  //    d_outPreTag =
  //        d_outPath + "/" + config["Output"]["Filename"].as<std::string>();
  //  else
  d_outPreTag = d_outPath + "/";

  // maximum number of output files to process
  d_dtN = d_dataManager_p->getModelDeckP()->d_Nt / d_outputDeck_p->d_dtOutCriteria;

  // read global start and end output step if provided
  d_dtStart = 1;
  if (config["Dt_Start"]) d_dtStart = config["Dt_Start"].as<int>();
  d_dtEnd = d_dtN;
  if (config["Dt_End"]) d_dtEnd = config["Dt_End"].as<int>();

  // handle change in output interval
  d_dtOutChange = d_dataManager_p->getModelDeckP()->d_Nt;
  if (config["Dt_Out_Change"])
    d_dtOutChange = config["Dt_Out_Change"].as<size_t>();

  // create mesh
  std::cout << "PP_fe2D: Creating mesh.\n";
  d_dataManager_p->setMeshP(new fe::Mesh(d_input_p->getMeshDeck()));

  // material deck and material
  {
    std::cout << "PP_fe2D: Initializing material object.\n";
    auto &material_deck = *d_input_p->getMaterialDeck();
    if (material_deck.d_materialType == "RNPBond")
      d_material_p = new material::pd::RNPBond(&material_deck, d_dataManager_p);
    // else if (material_deck.d_materialType == "PMBBond")
    // d_material_p = new material::pd::PmbMaterial(
    //     &material_deck, d_modelDeck_p->d_dim, d_modelDeck_p->d_horizon);
    // else if (material_deck.d_materialType == "PDElasticBond")
    //   d_material_p = new material::pd::PdElastic(
    //      &material_deck, d_modelDeck_p->d_dim, d_modelDeck_p->d_horizon);
    // else if (material_deck.d_materialType == "PDState")
    // d_material_p = new material::pd::PdState(
    //   &material_deck, d_modelDeck_p->d_dim, d_modelDeck_p->d_horizon);
  }

  d_matDeck_p = d_input_p->getMaterialDeck();

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

  // check if damage data file is provided
  if (config["Compute"]["File_Z"])
    d_fileZ = config["Compute"]["File_Z"].as<std::string>();

  if (config["Compute"]["Tag_Z"])
    d_tagZ = config["Compute"]["Tag_Z"].as<std::string>();
  else
    d_tagZ = "Damage_Z";

  if (config["Compute"]["Take_U_Plus"])
    d_uPlus = config["Compute"]["Take_U_Plus"].as<bool>();

  // read compute instruction
  auto num_compute = config["Compute"]["Sets"].as<size_t>();
  for (size_t c = 0; c < num_compute; c++) {
    std::string set = "Set_" + std::to_string(c + 1);
    auto data = tools::pp::InstructionData();
    readComputeInstruction(set, &data);
    if (data.d_start == -1) data.d_start = d_dtStart;
    if (data.d_end == -1) data.d_end = d_dtEnd;

    d_computeData.emplace_back(data);
  }

  // copy edge crack set for findCrackTip() function
  for (auto &d : d_computeData) {
    if (d.d_findCrackTip_p) {
      d.d_findCrackTip_p->d_cracks = d_fractureDeck_p->d_cracks;

      // check which point, top/bottom, we need to track
      auto bbox = d_dataManager_p->getMeshP()->getBoundingBox();
      for (auto &ck : d.d_findCrackTip_p->d_cracks) {
        // point should be in smaller box inside bounding box
        ck.d_trackt = util::geometry::isPointInsideRectangle(
            ck.d_pt, bbox.first[0] + d_dataManager_p->getModelDeckP()->d_horizon,
            bbox.second[0] - d_dataManager_p->getModelDeckP()->d_horizon,
            bbox.first[1] + d_dataManager_p->getModelDeckP()->d_horizon,
            bbox.second[1] - d_dataManager_p->getModelDeckP()->d_horizon);

        ck.d_trackb = util::geometry::isPointInsideRectangle(
            ck.d_pb, bbox.first[0] + d_dataManager_p->getModelDeckP()->d_horizon,
            bbox.second[0] - d_dataManager_p->getModelDeckP()->d_horizon,
            bbox.first[1] + d_dataManager_p->getModelDeckP()->d_horizon,
            bbox.second[1] - d_dataManager_p->getModelDeckP()->d_horizon);
      }  // modify crack track status
    }
  }

  // read crack tip data for J integral calculation
  for (auto &d : d_computeData) {
    auto data = d.d_computeJInt_p;
    if (data) {
      // read data from provided file
      readCrackTipData(data->d_crackTipFile, data->d_crackId,
                       &(data->d_crackTipData), data->d_isCrackInclined);

      // find the start, end, and interval of time in data
      if (!data->d_crackTipData.empty()) {
        d.d_start = data->d_crackTipData[0].d_n;
        d.d_end = data->d_crackTipData[data->d_crackTipData.size() - 1].d_n;
        if (data->d_crackTipData.size() > 1)
          d.d_interval =
              data->d_crackTipData[1].d_n - data->d_crackTipData[0].d_n;
      }
    }
  }
}

void tools::pp::Compute::finalize() {
  // close open file
  for (const auto &c : d_computeData) {
    if (c.d_computeJInt_p) {
      if (c.d_computeJInt_p->d_file) fclose(c.d_computeJInt_p->d_file);

      if (c.d_computeJInt_p->d_fileNew) fclose(c.d_computeJInt_p->d_fileNew);
    }

    if (c.d_findCrackTip_p) {
      if (c.d_findCrackTip_p->d_filet) fclose(c.d_findCrackTip_p->d_filet);
      if (c.d_findCrackTip_p->d_fileb) fclose(c.d_findCrackTip_p->d_fileb);
    }
  }
}

void tools::pp::Compute::safeExit(const std::string &err_message,
                                  rw::writer::Writer *writer) {
  // finalize
  finalize();

  // close writer if open
  if (writer != nullptr) writer->close();

  // print message
  std::cerr << err_message << "\n";

  exit(1);
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

  if (config["Compute"][set]["File_Format"])
    data->d_outFormat = config["Compute"][set]["File_Format"].as<std::string>();

  // read format of output file

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
                 " Dt_Interval of compute "
              << set << ".\n";
    exit(1);
  }

  // check if we need to perform calculations in reference or current
  // configuration

  // check in global space first
  if (config["Compute"]["Calculate_In_Ref_Config"])
    data->d_calculateInRefConfig =
        config["Compute"]["Calculate_In_Ref_Config"].as<bool>();

  // now check in local space (override the value provided in global space)
  if (config["Compute"][set]["Calculate_In_Ref_Config"])
    data->d_calculateInRefConfig =
        config["Compute"][set]["Calculate_In_Ref_Config"].as<bool>();

  // Scale displacement
  if (config["Compute"][set]["Scale_U_Ouptut"]) {
    if (!data->d_transformU_p)
      data->d_transformU_p = new tools::pp::TransformU();
    data->d_transformU_p->d_scale =
        config["Compute"][set]["Scale_U_Ouptut"].as<double>();
  }

  // Compute damage at nodes
  if (config["Compute"][set]["Damage_Z"]) {
    data->d_damageAtNodes = config["Compute"][set]["Damage_Z"].as<bool>();
    d_needDamageZ = true;
    d_needNeighborList = true;

    if (data->d_damageAtNodes && !d_fileZ.empty()) {
      std::cerr << "Error: Damage file is specified in Compute:File_Z and at "
                << "the same time damage is asked to be computed in compute "
                << " set = " << set << ".\n";
      exit(1);
    }
  }

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
        for (auto j : pt) locs.push_back(j.as<double>());

        if (locs.size() == 2) locs.push_back(0.0);

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

    if (config["Compute"][set]["Crack_Tip"]["Max_Z_Allowed"])
      data->d_findCrackTip_p->d_maxZAllowed =
          config["Compute"][set]["Crack_Tip"]["Max_Z_Allowed"].as<double>();

    if (config["Compute"][set]["Crack_Tip"]["Min_Z_Allowed"])
      data->d_findCrackTip_p->d_minZAllowed =
          config["Compute"][set]["Crack_Tip"]["Min_Z_Allowed"].as<double>();

    d_needDamageZ = true;
    d_needNeighborList = true;
  }

  // J integral
  if (config["Compute"][set]["J_Integral"]) {
    if (!data->d_computeJInt_p)
      data->d_computeJInt_p = new tools::pp::ComputeJIntegral();

    d_needDamageZ = true;
    d_needNeighborList = true;

    auto e = config["Compute"][set]["J_Integral"];
    if (!e["Crack_Orient"]) {
      data->d_computeJInt_p->d_crackOrient = d_fractureDeck_p->d_cracks[0].d_o;
    } else {
      data->d_computeJInt_p->d_crackOrient = e["Crack_Orient"].as<int>();

      if (data->d_computeJInt_p->d_crackOrient !=
          d_fractureDeck_p->d_cracks[0].d_o) {
        if (data->d_computeJInt_p->d_crackOrient != 0) {
          std::cerr << "Error: Setting crack orient as "
                    << data->d_computeJInt_p->d_crackOrient
                    << ", however the crack orient in fracture deck is "
                    << d_fractureDeck_p->d_cracks[0].d_o
                    << ".\n Crack could initially be horizontal or vertical "
                       "and later it is arbitrarily oriented.\n "
                       "However the case when crack is opposite orientation "
                       "of what it's initial orientation seems wrong.\n";
          exit(1);
        } else {
          std::cout << "Warning: Setting crack orient as "
                    << data->d_computeJInt_p->d_crackOrient
                    << ", however the crack orient in fracture deck is "
                    << d_fractureDeck_p->d_cracks[0].d_o << ".\n";
        }
      }
    }

    // see if crack inclined flag is provided. This flag is important if
    // crack propagation is not along horizontal/vertical direction but along
    // arbitrary direction.
    if (e["Is_Crack_Inclined"])
      data->d_computeJInt_p->d_isCrackInclined =
          e["Is_Crack_Inclined"].as<bool>();

    if (e["Crack_Id"])
      data->d_computeJInt_p->d_crackId = e["Crack_Id"].as<int>();

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

      data->d_computeJInt_p->d_contourGiven = false;

    } else if (e["Contour"]) {
      std::vector<double> locs;
      for (auto f : e["Contour"]) locs.push_back(f.as<double>());

      if (locs.size() != 6) {
        std::cerr << "Error: Need 6 parameters corresponding to x,y,z "
                     "coordinate of two corner points to define the contour.\n";
        exit(1);
      }

      data->d_computeJInt_p->d_contour.first =
          util::Point3(locs[0], locs[1], locs[2]);
      data->d_computeJInt_p->d_contour.second =
          util::Point3(locs[3], locs[4], locs[5]);

      data->d_computeJInt_p->d_contourFactor = {-1., -1.};

      data->d_computeJInt_p->d_contourGiven = true;
    } else {
      std::cerr
          << "Error: Either factors to create contour for J integral or"
             " contour rectangle is required for J-integral calculation.\n";
      exit(1);
    }

    if (e["Set_V_Lateral_Comp_Zero"])
      data->d_computeJInt_p->d_setLateralCompVZero =
          e["Set_V_Lateral_Comp_Zero"].as<bool>();

    if (e["Set_U_Lateral_Comp_Zero"])
      data->d_computeJInt_p->d_setLateralCompUZero =
          e["Set_U_Lateral_Comp_Zero"].as<bool>();

    if (e["Set_X_Lateral_Comp"])
      data->d_computeJInt_p->d_setLateralCompX =
          e["Set_X_Lateral_Comp"].as<double>();
    else {
      if (data->d_computeJInt_p->d_setLateralCompUZero) {
        std::cerr << "Error: Expecting value of lateral component of crack "
                     "tip location as the flag Set_U_Lateral_Comp_Zero set to"
                     " true.\n";
        exit(1);
      }
    }
  }
}

void tools::pp::Compute::readCrackTipData(
    const std::string &filename, int crack_id,
    std::vector<tools::pp::CrackTipData> *data, bool is_inclined) {
  data->clear();
  if (!is_inclined) {
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
  } else {
    // expected format of file:
    // <crack id>, <output step>, <tip x>, <tip y>, <tip vx>, <tip vy>, <tip
    // dx>, <tip dy>
    io::CSVReader<8> in(filename);
    in.read_header(io::ignore_extra_column, "'crack_id'", "'dt_out'", "'x'",
                   "'y'", "'vx'", "'vy'", "'dx'", "'dy'");

    double px, py, vx, vy, dx, dy;
    int n, ci;
    while (in.read_row(ci, n, px, py, vx, vy, dx, dy)) {
      if (ci == crack_id)
        data->emplace_back(tools::pp::CrackTipData(n, util::Point3(px, py, 0.),
                                                   util::Point3(vx, vy, 0.),
                                                   util::Point3(dx, dy, 0.)));
    }
  }
}

void tools::pp::Compute::initWriter(rw::writer::Writer *writer,
                                    std::vector<util::Point3> *u) {
  if (d_writerReady) return;

  writer->open(d_outFilename, d_currentData->d_outFormat,
               d_currentData->d_compressType);
  // append mesh (check if only nodes need to be written)
  if (d_currentData->d_outOnlyNodes)
    writer->appendNodes(d_dataManager_p->getMeshP()->getNodesP(), u);
  else
    writer->appendMesh(d_dataManager_p->getMeshP()->getNodesP(),
                       d_dataManager_p->getMeshP()->getElementType(),
                       d_dataManager_p->getMeshP()->getElementConnectivitiesP(),
                       u);

  // append current time
  writer->addTimeStep(d_time);

  // set flag
  d_writerReady = true;
}

//
// compute methods
//
void tools::pp::Compute::transformU(rw::writer::Writer *writer) {
  if (!d_currentData->d_transformU_p) return;

  std::vector<util::Point3> u_temp(d_dataManager_p->getMeshP()->getNumNodes(),
                                   util::Point3());
  auto scale = d_currentData->d_transformU_p->d_scale;
  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_dataManager_p->getMeshP()->getNumNodes(),
      [&u_temp, scale, this](boost::uint64_t i) {
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

void tools::pp::Compute::transformV(rw::writer::Writer *writer) {
  if (!d_currentData->d_transformV_p) return;

  std::vector<util::Point3> v_mark;

  auto data = d_currentData->d_transformV_p;

  // mark v operation (this preceeds symmetrizing of v)
  if (data->d_markVAsZero) {
    v_mark = d_v;

    if (data->d_markVInRectGiven) {
      auto f = hpx::parallel::for_loop(
          hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
          d_dataManager_p->getMeshP()->getNumNodes(),
          [&v_mark, data, this](boost::uint64_t i) {
            if (util::geometry::isPointInsideRectangle(
                    d_dataManager_p->getMeshP()->getNode(i),
                    data->d_markVRect.first.d_x, data->d_markVRect.second.d_x,
                    data->d_markVRect.first.d_y, data->d_markVRect.second.d_y))
              v_mark[i] = util::Point3();
          });
      f.get();
    }

    if (!data->d_markVPts.empty())
      for (auto x : data->d_markVPts) {
        size_t i_found;
        if (data->d_markVPtsAreInCurrentConfig)
          i_found = findNode(x, d_dataManager_p->getMeshP()->getNodesP(), &d_u);
        else
          i_found = findNode(x, d_dataManager_p->getMeshP()->getNodesP());

        // modify v_new
        v_mark[i_found] = util::Point3();
      }

    if (!data->d_markVNodes.empty())
      for (auto i : data->d_markVNodes) v_mark[i] = util::Point3();

    if (!d_writerReady) initWriter(writer, &d_u);

    // append velocity
    writer->appendPointData("Mark_Velocity", &v_mark);
  }  // mark v operation

  // symmetrize v operation
  if (data->d_symmetrizeV) {
    // use v_mark for modification
    // Reset it to current velocity if we are not combining this
    // with mark_v operation

    if (!data->d_combineMarkV) v_mark = d_v;

    auto f = hpx::parallel::for_loop(
        hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
        d_dataManager_p->getMeshP()->getNumNodes(),
        [&v_mark, data, this](boost::uint64_t i) {
          auto x = d_dataManager_p->getMeshP()->getNode(i);

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
            size_t i_found =
                findNode(search_x, d_dataManager_p->getMeshP()->getNodesP());

            // write velocity
            v_mark[i] = v_mark[i_found];
            if (data->d_symmAxis == "y") v_mark[i].d_x *= -1.;
            if (data->d_symmAxis == "x") v_mark[i].d_y *= -1.;
          }
        });  // parallel for loop
    f.get();

    if (!d_writerReady) initWriter(writer, &d_u);

    // append velocity
    writer->appendPointData("Symm_Velocity", &v_mark);
  }
}

void tools::pp::Compute::computeStrain(rw::writer::Writer *writer) {
  if (!d_currentData->d_compStrain_p) return;

  auto data = d_currentData->d_compStrain_p;
  if (!data->d_computeStrain) return;

  std::vector<util::SymMatrix3> strain(
      d_dataManager_p->getMeshP()->getNumElements(), util::SymMatrix3());
  std::vector<util::SymMatrix3> stress(
      d_dataManager_p->getMeshP()->getNumElements(), util::SymMatrix3());
  std::vector<float> magS;

  // get Quadrature
  fe::BaseElem *quad;
  if (d_dataManager_p->getMeshP()->getElementType() == util::vtk_type_triangle)
    quad = new fe::TriElem(1);
  else if (d_dataManager_p->getMeshP()->getElementType() == util::vtk_type_quad)
    quad = new fe::QuadElem(1);
  else {
    std::cerr << "Error: Can not compute strain/stress as the element "
                 "type is not implemented.\n";
    exit(1);
  }

  // compute strain and stress
  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_dataManager_p->getMeshP()->getNumElements(),
      [&strain, &stress, data, quad, this](boost::uint64_t e) {
        auto ssn = util::SymMatrix3();
        auto sss = util::SymMatrix3();

        // get ids of nodes of element, coordinate of nodes, 1st order
        // quad data, and first quad data
        auto id_nds = d_dataManager_p->getMeshP()->getElementConnectivity(e);
        auto nds = d_dataManager_p->getMeshP()->getElementConnectivityNodes(e);
        auto qds = quad->getQuadDatas(nds);
        auto qd0 = qds[0];

        // compute strain in xy plane
        for (size_t i = 0; i < id_nds.size(); i++) {
          auto id = id_nds[i];
          auto ui = d_u[id];

          ssn(0, 0) += ui.d_x * qd0.d_derShapes[i][0];
          ssn(1, 1) += ui.d_y * qd0.d_derShapes[i][1];
          ssn(0, 1) += 0.5 * ui.d_x * qd0.d_derShapes[i][1];
          ssn(0, 1) += 0.5 * ui.d_y * qd0.d_derShapes[i][0];
        }

        if (d_matDeck_p->d_isPlaneStrain)
          ssn(2, 2) = -d_matDeck_p->d_matData.d_nu * (ssn(0, 0) + ssn(1, 1)) /
                      (1. - d_matDeck_p->d_matData.d_nu);

        // compute stress
        auto trace = ssn(0, 0) + ssn(1, 1) + ssn(2, 2);
        sss(0, 0) = d_matDeck_p->d_matData.d_lambda * trace * 1. +
                    2. * d_matDeck_p->d_matData.d_mu * ssn(0, 0);
        sss(0, 1) = 2. * d_matDeck_p->d_matData.d_mu * ssn(0, 1);
        sss(0, 2) = 2. * d_matDeck_p->d_matData.d_mu * ssn(0, 2);

        sss(1, 1) = d_matDeck_p->d_matData.d_lambda * trace * 1. +
                    2. * d_matDeck_p->d_matData.d_mu * ssn(1, 1);
        sss(1, 2) = 2. * d_matDeck_p->d_matData.d_mu * ssn(1, 2);
        if (!d_matDeck_p->d_isPlaneStrain)
          sss(2, 2) = d_matDeck_p->d_matData.d_nu * (sss(0, 0) + sss(1, 1));

        strain[e] = ssn;
        stress[e] = sss;
      });  // parallel loop over elements
  f.get();

  // compute magnitude of strain
  if (data->d_magStrainTensor) {
    magS = std::vector<float>(strain.size(), 0.);
    auto f2 = hpx::parallel::for_loop(
        hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
        d_dataManager_p->getMeshP()->getNumElements(),
        [&magS, strain, data](boost::uint64_t e) {
          if (data->d_magStrainComp.empty()) {
            for (size_t i = 0; i < 6; i++) magS[e] = std::abs(strain[e].get(i));
          } else if (data->d_magStrainComp == "xx") {
            magS[e] = std::abs(strain[e](0, 0));
          } else if (data->d_magStrainComp == "yy") {
            magS[e] = std::abs(strain[e](1, 1));
          }
        });
    f2.get();
  }

  // output strain/stress data
  if (!d_currentData->d_outOnlyNodes) {
    // append mesh
    if (!d_writerReady) initWriter(writer, &d_u);

    writer->appendCellData("Strain_Tensor", &strain);
    writer->appendCellData("Stress_Tensor", &stress);
    if (data->d_magStrainTensor) writer->appendCellData("Mag_Strain", &magS);

    // mark magnitude of strain if asked
    if (!data->d_markMagStrainCells.empty() && data->d_magStrainTensor) {
      for (auto cell : data->d_markMagStrainCells)
        magS[cell.first] = cell.second;

      writer->appendCellData("Mark_Mag_Strain", &magS);
    }
  } else {
    // compute 1st order quad points and store them (quad points in
    // current configuration)
    std::vector<util::Point3> elem_quads = std::vector<util::Point3>(
        d_dataManager_p->getMeshP()->getNumElements(), util::Point3());

    auto f2 = hpx::parallel::for_loop(
        hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
        d_dataManager_p->getMeshP()->getNumElements(),
        [&elem_quads, quad, this](boost::uint64_t e) {
          auto nds = d_dataManager_p->getMeshP()->getElementConnectivity(e);
          std::vector<util::Point3> nds_current(nds.size(), util::Point3());
          for (size_t j = 0; j < nds.size(); j++)
            nds_current[j] =
                d_dataManager_p->getMeshP()->getNode(nds[j]) + d_u[nds[j]];

          auto qds = quad->getQuadPoints(nds_current);
          // store first quad point
          elem_quads[e] = qds[0].d_p;
        });
    f2.get();

    // create unstructured vtk output
    std::string fname = d_outPreTag + d_currentData->d_tagFilename + "_quads_" +
                        std::to_string(d_nOut);
    auto writer1 = rw::writer::Writer(fname, d_currentData->d_outFormat,
                                      d_currentData->d_compressType);
    writer1.appendNodes(&elem_quads);
    writer1.appendPointData("Strain_Tensor", &strain);
    writer1.appendPointData("Stress_Tensor", &stress);
    if (data->d_magStrainTensor) writer1.appendPointData("Mag_Strain", &magS);

    // mark magnitude of strain if asked
    if (!data->d_markMagStrainCells.empty() && data->d_magStrainTensor) {
      for (auto cell : data->d_markMagStrainCells)
        magS[cell.first] = cell.second;

      writer1.appendPointData("Mark_Mag_Strain", &magS);
    }

    writer1.addTimeStep(d_time);
    writer1.close();
  }
}

void tools::pp::Compute::computeDamage(rw::writer::Writer *writer,
                                       std::vector<double> *Z, bool perf_out) {
  //  if (!d_currentData->d_damageAtNodes)
  //    return;

  if (Z->size() != d_dataManager_p->getMeshP()->getNumNodes())
    Z->resize(d_dataManager_p->getMeshP()->getNumNodes());

  if (d_dataManager_p->getNeighborP() == nullptr)
    safeExit(
        "Error: Need neighbor list to compute damage. This should have "
        "been created at the beginning.\n");

  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_dataManager_p->getMeshP()->getNumNodes(),
      [&Z, this](boost::uint64_t i) {
        auto xi = d_dataManager_p->getMeshP()->getNode(i);

        double locz = 0.;
        auto i_neighs = d_dataManager_p->getNeighborP()->getNeighbors(i);
        for (const auto &j : i_neighs) {
          if (util::compare::definitelyGreaterThan(
                  xi.dist(d_dataManager_p->getMeshP()->getNode(j)),
                  d_dataManager_p->getModelDeckP()->d_horizon) ||
              j == i)
            continue;

          auto xj = d_dataManager_p->getMeshP()->getNode(j);
          if (util::compare::definitelyGreaterThan(xj.dist(xi), 1.0E-10)) {
            auto Sr = std::abs(d_material_p->getS(xj - xi, d_u[j] - d_u[i])) /
                      d_material_p->getSc(xj.dist(xi));

            if (util::compare::definitelyLessThan(locz, Sr)) locz = Sr;
          }
        }  // loop over neighbors

        (*Z)[i] = locz;
      });  // parallel loop over nodes
  f.get();

  if (!perf_out) return;

  if (!d_writerReady) initWriter(writer, &d_u);

  writer->appendPointData("Damage_Z", Z);
}

void tools::pp::Compute::findCrackTip(std::vector<double> *Z,
                                      rw::writer::Writer *writer) {
  auto data = d_currentData->d_findCrackTip_p;
  if (!data) return;

  // compute crack tip location and crack tip velocity
  updateCrack(d_time, Z);

  // perform output
  crackOutput();
}

void tools::pp::Compute::computeJIntegral() {
  auto data = d_currentData->d_computeJInt_p;
  const bool calc_in_ref = d_currentData->d_calculateInRefConfig;
  if (!data) return;

  // compute necessary quantities if the material is state-based
  if (d_material_p->isStateActive()) {
    if (d_dataManager_p->getFractureP() == nullptr) {
      d_dataManager_p->setFractureP(new geometry::Fracture(
          d_input_p->getFractureDeck(),
          d_dataManager_p->getMeshP()->getNodesP(),
          d_dataManager_p->getNeighborP()->getNeighborsListP()));
    }

    // initialization
    if (d_thetaX.size() != d_dataManager_p->getMeshP()->getNumNodes())
      d_thetaX =
          std::vector<double>(d_dataManager_p->getMeshP()->getNumNodes(), 0.);

    /*
    if (d_material_p->name() == "PDState") {
      if (d_mX.size() != d_mesh_p->getNumNodes()) {
        d_mX = std::vector<double>(d_mesh_p->getNumNodes(), 0.);

        material::computeStateMx(
            d_mesh_p->getNodes(), d_mesh_p->getNodalVolumes(),
            d_neighbor_p->getNeighborsList(), d_mesh_p->getMeshSize(),
            d_material_p, d_mX, true);
      }
    }
    */

    // update theta for given displacement
    /*
    if (d_material_p->name() == "RNPState")
      material::computeHydrostaticStrain(
          d_mesh_p->getNodes(), d_u, d_mesh_p->getNodalVolumes(),
          d_neighbor_p->getNeighborsList(), d_mesh_p->getMeshSize(),
          d_material_p, d_fracture_p, d_thetaX, d_modelDeck_p->d_dim, true);
    else if (d_material_p->name() == "PDState") {
      // need to update the fracture state of bonds
      material::updateBondFractureData(d_mesh_p->getNodes(), d_u,
                                       d_neighbor_p->getNeighborsList(),
                                       d_material_p, d_fracture_p, true);

      material::computeStateThetax(
          d_mesh_p->getNodes(), d_u, d_mesh_p->getNodalVolumes(),
          d_neighbor_p->getNeighborsList(), d_mesh_p->getMeshSize(),
          d_material_p, d_fracture_p, d_mX, d_thetaX, true);
    }
    */
  }

  // std::cout << "computeJIntegral processing step = " << d_nOut << "\n";

  // to hold energy into crack
  auto j_energy = JEnergy();

  // get crack tip data
  auto ctip = data->d_crackTipData[(d_nOut - d_currentData->d_start) /
                                   d_currentData->d_interval];

  // If this is inclined crack, then we are already provided the crack direction
  // in the input file. If this is not a inclined crack than compute the
  // direction based on orientation of crack
  if (!data->d_isCrackInclined) {
    if (data->d_crackOrient == -1)
      ctip.d_d = util::Point3(0., 1., 0.);
    else if (data->d_crackOrient == 1)
      ctip.d_d = util::Point3(1., 0., 0.);
  }

  if (ctip.d_n != d_nOut) {
    oss.str("");
    oss << "Error: Output step (" << ctip.d_n << ") of crack tip data is"
        << " not matching current output step (" << d_nOut << ").\n";
    safeExit(oss.str());
  }

  // set lateral component of crack velocity zero
  if (data->d_setLateralCompVZero) {
    if (data->d_crackOrient == -1)
      ctip.d_v.d_x = 0.;
    else if (data->d_crackOrient == 1)
      ctip.d_v.d_y = 0.;
  }

  // set lateral component of displacement of crack tip zero
  if (data->d_setLateralCompUZero) {
    if (data->d_crackOrient == -1)
      ctip.d_p.d_x = data->d_setLateralCompX;
    else if (data->d_crackOrient == 1)
      ctip.d_p.d_y = data->d_setLateralCompX;
  }

  // Schematic for horizontal crack (similar for vertical crack)
  //
  // Case when vertical or horizontal crack
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
  // Case when inclined crack
  //
  //                         D                    C
  //                         + + + + + + + + + + +
  //                         +                   +
  //                         +    o              +
  //                         + /                 +
  //                      /  +                   +
  //                  /      +                   +
  //              /          + + + + + + + + + + +
  //                        A                    B
  //
  // Contour is formed by lines A-B, B-C, C-D, D-A
  std::pair<util::Point3, util::Point3> cd;
  if (data->d_contourGiven)
    cd = data->d_contour;
  else {
    double tolx = 0.5 * data->d_contourFactor[0] * d_dataManager_p->getModelDeckP()->d_horizon;
    double toly = 0.5 * data->d_contourFactor[1] * d_dataManager_p->getModelDeckP()->d_horizon;
    cd = std::make_pair(
        util::Point3(ctip.d_p.d_x - tolx, ctip.d_p.d_y - toly, 0.),
        util::Point3(ctip.d_p.d_x + tolx, ctip.d_p.d_y + toly, 0.));
  }

  // compute nodes and elements list for search
  std::vector<size_t> search_nodes;
  std::vector<size_t> search_elems;
  listElemsAndNodesInDomain(cd,
                            d_dataManager_p->getModelDeckP()->d_horizon +
                                2. * d_dataManager_p->getMeshP()->getMeshSize(),
                            d_dataManager_p->getMeshP()->getMeshSize(),
                            &search_nodes, &search_elems, calc_in_ref);

  //
  // Compute contour integral
  //
  // create second order quadrature class for 1-d line element
  auto line_quad = fe::LineElem(2);
  auto h = d_dataManager_p->getMeshP()->getMeshSize();

  // in the expression of contour integrals, we have dot product of
  // normal to the edges in contour with the crack velocity direction.
  // Thus, if crack velocity is along horizontal direction, then we only need
  // to evaluate contour integrals on vertical edges D-A and B-C.
  // Similarly, if velocity is along the vertical direction then we only need
  // to evaluate contour integrals on horizontal edge A-B and C-D.
  //
  // However if this is inclined crack, then we have to integrate over all
  // four edges
  bool integrate_horizontal = true;

  // skip horizontal edge for horizontal crack
  if (data->d_crackOrient == 1) integrate_horizontal = false;

  // Debug
  //  Must delete after testing
  // d_material_p->getMaterialDeckP()->d_irreversibleBondBreak = false;

  // Assumption:
  //
  // For vertical crack, we assume (convention) it has +ve velocity in +y
  // direction
  // therefore n dot n_c = -1 on bottom edge A-B and n dot n_c = +1 on top
  // edge C-D
  //
  // For horizontal crack, we assume (convention) it has +ve velocity in +x
  // direction therefore n dot n_c = -1 on left edge D-A and n dot n_c = +1 on
  // right edge B-C

  // loop over horizontal and vertical edge of contour
  // E = 0 : horizontal edge A-B and C-D
  // E = 1 : vertical edge D-A and B-C

  double h_contour = 0.5 * h;
  for (size_t E = 0; E < 2; E++) {
    // skip horizontal edge if required
    if (E == 0 && !integrate_horizontal && !data->d_isCrackInclined) continue;

    // skip vertical edge if required
    if (E == 1 && integrate_horizontal && !data->d_isCrackInclined) continue;

    long int N = 0;
    if (E == 0) {
      // number of elements for horizontal edge
      N = (cd.second.d_x - cd.first.d_x) / h_contour;
      if (util::compare::definitelyLessThan(
              cd.first.d_x + double(N) * h_contour, cd.second.d_x))
        N++;
    } else {
      // number of elements for vertical edge
      N = (cd.second.d_y - cd.first.d_y) / h_contour;
      if (util::compare::definitelyLessThan(
              cd.first.d_y + double(N) * h_contour, cd.second.d_y))
        N++;
    }

    auto ced = ContourDataLambdaFn();
    ced.d_contourPdStrainEnergies = std::vector<double>(N, 0.);
    ced.d_contourPdStrainEnergiesRate = std::vector<double>(N, 0.);
    ced.d_contourKineticEnergiesRate = std::vector<double>(N, 0.);
    ced.d_contourElasticInternalWorksRate = std::vector<double>(N, 0.);

    // check
    if (E == 0 && data->d_crackOrient == 1 && !data->d_isCrackInclined)
      safeExit("For horizontal crack, we should not reach this point\n");
    if (E == 1 && data->d_crackOrient == -1 && !data->d_isCrackInclined)
      safeExit("For vertical crack, we should not reach this point\n");

    // Debug
    //  collect all quad points
    // std::vector<util::Point3> quad_points(2 * N, util::Point3());

    //    auto f = hpx::parallel::for_loop(
    //        hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
    //        N,
    //        [&ced, N, h, cd, ctip, &line_quad, search_nodes, search_elems,
    //         E, this, calc_in_ref](boost::uint64_t I) {
    for (size_t I = 0; I < N; I++) {
      double kinetic_energy_q = 0.;
      double pd_energy_q = 0.;
      double elastic_energy_q = 0.;
      util::Point3 dot_u_q = util::Point3();

      double pd_strain_energy = 0.;
      double pd_strain_energy_rate = 0.;
      double kinetic_energy_rate = 0.;
      double elastic_internal_work_rate = 0.;

      // normal to edge
      util::Point3 edge_normal = util::Point3();

      // dot product of normal to edge and crack velocity direction
      double n_dot_n_c = 0.;

      // dot product of normal to edge and crack velocity
      double n_dot_v_c = 0.;

      // line element
      auto x1 = 0.;
      auto x2 = 0.;
      if (E == 0) {
        // discretization of horizontal line
        x1 = cd.first.d_x + double(I) * h_contour;
        x2 = cd.first.d_x + double(I + 1) * h_contour;
        if (I == N - 1) x2 = cd.second.d_x;
      } else {
        // discretization of vertical line
        x1 = cd.first.d_y + double(I) * h_contour;
        x2 = cd.first.d_y + double(I + 1) * h_contour;
        if (I == N - 1) x2 = cd.second.d_y;
      }

      // get quadrature points
      auto qds = line_quad.getQuadPoints(std::vector<util::Point3>{
          util::Point3(x1, 0., 0.), util::Point3(x2, 0., 0.)});

      // loop over quad points
      for (auto qd : qds) {
        for (int top_side = 0; top_side < 2; top_side++) {
          // process data
          util::Point3 qp = qd.d_p;
          processQuadPointForContour(E, top_side, cd, ctip.d_v, ctip.d_d, qp,
                                     edge_normal, n_dot_n_c, n_dot_v_c);

          // quad_points[2 * I + top_side] = qp;

          // get energy density
          getContourContribJInt(qp, &search_nodes, &search_elems, edge_normal,
                                pd_energy_q, kinetic_energy_q, elastic_energy_q,
                                dot_u_q, calc_in_ref, ctip);

          // debug information
          if (false) {
            std::cout
                << "--------------------------------------------------------\n";
            std::cout << "E: " << E << ", I: " << I << ", qp: " << qp.printStr()
                      << ", qw: " << qd.d_w << ", top_side: " << top_side
                      << "\n";

            std::cout << "tip: " << ctip.d_p.printStr()
                      << ", v: " << ctip.d_v.printStr()
                      << ", d: " << ctip.d_d.printStr()
                      << ", edge_normal: " << edge_normal.printStr()
                      << ", n_dot_n_c: " << n_dot_n_c
                      << ", n_dot_v_c: " << n_dot_v_c << "\n";

            std::cout << "pd_energy_q: " << pd_energy_q
                      << ", kinetic_energy_q: " << kinetic_energy_q
                      << ", elastic_energy_q: " << elastic_energy_q
                      << ", pd_strain_energy: "
                      << pd_energy_q * n_dot_n_c * qd.d_w
                      << ", pd_strain_energy_rate: "
                      << pd_energy_q * n_dot_v_c * qd.d_w
                      << ", elastic_internal_work_rate: "
                      << elastic_energy_q * qd.d_w << "\n\n";
          }

          //          // Debug
          //          //  Must be deleted after testing
          //          if ((data->d_crackOrient == -1 and
          //               util::compare::definitelyGreaterThan(qd.d_p.d_x,
          //                                                    ctip.d_p.d_x))
          //                                                    or
          //              (data->d_crackOrient == 1 and
          //               util::compare::definitelyGreaterThan(qd.d_p.d_y,
          //                                                    ctip.d_p.d_y)))
          //                                                    {
          // pd energy
          pd_strain_energy += pd_energy_q * n_dot_n_c * qd.d_w;

          // pd energy rate
          pd_strain_energy_rate += pd_energy_q * n_dot_v_c * qd.d_w;

          // kinetic energy rate
          kinetic_energy_rate += kinetic_energy_q * n_dot_v_c * qd.d_w;

          // elastic internal work rate
          elastic_internal_work_rate += elastic_energy_q * qd.d_w;
          //          }
        }

      }  // loop over quad points

      // add energy
      ced.d_contourPdStrainEnergies[I] = pd_strain_energy;
      ced.d_contourPdStrainEnergiesRate[I] = pd_strain_energy_rate;
      ced.d_contourKineticEnergiesRate[I] = kinetic_energy_rate;
      ced.d_contourElasticInternalWorksRate[I] = elastic_internal_work_rate;
    }
    //    );
    //    f.get();

    // {
    //   // write to csv file and exit
    //   std::ofstream qd_pt_file;
    //   qd_pt_file.open(d_outPreTag + "quad_points.csv");
    //   for (auto p : quad_points)
    //     qd_pt_file << p.d_x << ", " << p.d_y << "\n";
    //   qd_pt_file.close();
    //   safeExit("Output debug quad points and exit");
    // }

    // sum energies
    j_energy.d_contourPdStrainEnergy +=
        util::methods::add(ced.d_contourPdStrainEnergies);
    j_energy.d_contourPdStrainEnergyRate +=
        util::methods::add(ced.d_contourPdStrainEnergiesRate);
    j_energy.d_contourKineticEnergyRate +=
        util::methods::add(ced.d_contourKineticEnergiesRate);
    j_energy.d_contourElasticInternalWorkRate +=
        util::methods::add(ced.d_contourElasticInternalWorksRate);
  }

  //
  // Compute work done by internal forces
  //
  // decompose search_nodes list in two parts: one list for nodes outside
  // domain A formed by contour and other for nodes on contour and inside
  // domain A.
  std::vector<size_t> search_node_comp;
  decomposeSearchNodes(cd, &search_nodes, &search_node_comp);

  {
    auto pd_internal_works = std::vector<double>(search_node_comp.size(), 0.);
    auto pd_internal_works_rate =
        std::vector<double>(search_node_comp.size(), 0.);

    // loop over nodes in compliment of domain A
    auto f = hpx::parallel::for_loop(
        hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
        search_node_comp.size(),
        [&pd_internal_works, &pd_internal_works_rate, h, cd, search_nodes,
         search_node_comp, this](boost::uint64_t i) {
          auto id = search_node_comp[i];
          auto xi = d_dataManager_p->getMeshP()->getNode(id);
          auto ui = d_u[id];
          auto vi = d_v[id];
          auto voli = d_dataManager_p->getMeshP()->getNodalVolume(id);

          double pd_internal_work = 0.;
          double pd_internal_work_rate = 0.;

          // loop over nodes in domain A
          for (auto j : search_nodes) {
            auto xj = d_dataManager_p->getMeshP()->getNode(j);
            auto rji = xj.dist(xi);
            if (util::compare::definitelyGreaterThan(
                    rji, d_dataManager_p->getModelDeckP()->d_horizon) ||
                j == id)
              continue;

            auto v_sum = d_v[j] + vi;
            auto xji = xj - xi;
            auto Sji = d_material_p->getS(xji, d_u[j] - ui);

            // get corrected volume of node j
            auto volj = d_dataManager_p->getMeshP()->getNodalVolume(j);
            if (util::compare::definitelyGreaterThan(
                    rji, d_dataManager_p->getModelDeckP()->d_horizon - 0.5 * h))
              volj *= (d_dataManager_p->getModelDeckP()->d_horizon + 0.5 * h - rji) / h;

            // get bond force
            bool fracture_state = false;
            auto ef = d_material_p->getBondEF(i, j);

            // pd internal work rate
            // need factor half (see the formula for internal work rate and
            // formula for value returned by getBondEF())
            pd_internal_work_rate +=
                0.5 * ef.second * volj *
                (xji.d_x * v_sum.d_x + xji.d_y * v_sum.d_y +
                 xji.d_z * v_sum.d_z) /
                rji;
          }  // loop over neighboring nodes

          pd_internal_works[i] += pd_internal_work * voli;
          pd_internal_works_rate[i] += pd_internal_work_rate * voli;
        });
    f.get();

    // add energies
    j_energy.d_pdInternalWork = util::methods::add(pd_internal_works);
    j_energy.d_pdInternalWorkRate = util::methods::add(pd_internal_works_rate);
  }

  //
  // Compute total pd strain and kinetic energy inside the domain formed by
  // contour
  //

  // in outer loop, we loop over nodes inside contour
  // in inner loop, we loop over nodes inside contour plus nodes in the
  // horizon radius of boundary of contour outside contour
  {
    auto pd_strain_energies = std::vector<double>(search_nodes.size(), 0.);
    auto kinetic_energies = std::vector<double>(search_nodes.size(), 0.);

    // loop over nodes in compliment of domain A
    auto f = hpx::parallel::for_loop(
        hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
        search_nodes.size(),
        [&pd_strain_energies, &kinetic_energies, h, cd, search_nodes,
         search_node_comp, this](boost::uint64_t i) {
          auto id = search_nodes[i];
          auto xi = d_dataManager_p->getMeshP()->getNode(id);
          auto ui = d_u[id];
          auto vi = d_v[id];
          auto voli = d_dataManager_p->getMeshP()->getNodalVolume(id);

          // kinetic energy
          kinetic_energies[i] =
              0.5 * d_material_p->getDensity() * vi.dot(vi) * voli;

          // compute pd strain energy
          double strain_energy = 0.;

          // add peridynamic energy
          auto i_neighs = d_dataManager_p->getNeighborP()->getNeighbors(id);
          for (auto j : i_neighs) {
            auto xj = d_dataManager_p->getMeshP()->getNode(j);
            if (util::compare::definitelyGreaterThan(
                    xj.dist(xi), d_dataManager_p->getModelDeckP()->d_horizon) ||
                util::compare::definitelyLessThan(xj.dist(xi), 1.0E-10))
              continue;

            // get displacement and strain
            auto uj = d_u[j];
            auto rji = xj.dist(xi);
            auto Sji = d_material_p->getS(xj - xi, uj - ui);

            // get volume correction
            auto volj = d_dataManager_p->getMeshP()->getNodalVolume(j);
            if (util::compare::definitelyGreaterThan(
                    rji, d_dataManager_p->getModelDeckP()->d_horizon -
                             0.5 * d_dataManager_p->getMeshP()->getMeshSize()))
              volj *= (d_dataManager_p->getModelDeckP()->d_horizon +
                       0.5 * d_dataManager_p->getMeshP()->getMeshSize() - rji) /
                      d_dataManager_p->getMeshP()->getMeshSize();

            bool fracture_state = false;
            auto ef = d_material_p->getBondEF(i, j);

            // add contribution to energy
            strain_energy += ef.second * volj;
          }  // loop over nodes for pd energy density

          pd_strain_energies[i] += strain_energy * voli;
        });
    f.get();

    // add energies
    j_energy.d_pdStrainEnergyInsideContour =
        util::methods::add(pd_strain_energies);
    j_energy.d_kineticEnergyInsideContour =
        util::methods::add(kinetic_energies);
  }

  //
  // Compute peridynamic fracture energy
  //

  // loop over nodes and compute fracture energy
  {
    auto pd_fracture_energies =
        std::vector<double>(d_dataManager_p->getMeshP()->getNumNodes(), 0.);

    // loop over nodes in compliment of domain A
    auto f = hpx::parallel::for_loop(
        hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
        d_dataManager_p->getMeshP()->getNumNodes(),
        [&pd_fracture_energies, this](boost::uint64_t i) {
          auto xi = d_dataManager_p->getMeshP()->getNode(i);
          auto ui = d_u[i];
          auto voli = d_dataManager_p->getMeshP()->getNodalVolume(i);
          auto Zi = this->d_Z[i];

          // compute pd strain energy
          double fracture_energy = 0.;

          // add peridynamic energy
          auto i_neighs = d_dataManager_p->getNeighborP()->getNeighbors(i);
          for (auto j : i_neighs) {
            auto xj = d_dataManager_p->getMeshP()->getNode(j);
            if (util::compare::definitelyGreaterThan(
                    xj.dist(xi), d_dataManager_p->getModelDeckP()->d_horizon) ||
                util::compare::definitelyLessThan(xj.dist(xi), 1.0E-10))
              continue;

            // get displacement and strain
            auto uj = d_u[j];
            auto rji = xj.dist(xi);
            auto Sji = d_material_p->getS(xj - xi, uj - ui);

            // get volume correction
            auto volj = d_dataManager_p->getMeshP()->getNodalVolume(j);
            if (util::compare::definitelyGreaterThan(
                    rji, d_dataManager_p->getModelDeckP()->d_horizon -
                             0.5 * d_dataManager_p->getMeshP()->getMeshSize()))
              volj *= (d_dataManager_p->getModelDeckP()->d_horizon +
                       0.5 * d_dataManager_p->getMeshP()->getMeshSize() - rji) /
                      d_dataManager_p->getMeshP()->getMeshSize();

            bool fracture_state = false;
            auto ef = d_material_p->getBondEF(i, j);

            // add contribution to energy
            fracture_energy += ef.second * volj;
          }  // loop over nodes for pd energy density

          if (util::compare::definitelyGreaterThan(Zi, 1.0 - 1.0E-10))
            pd_fracture_energies[i] += fracture_energy * voli;
        });
    f.get();

    // add energies
    j_energy.d_pdFractureEnergy = util::methods::add(pd_fracture_energies);
  }

  // compute remaining energies
  {
    auto vmag = ctip.d_v.length();

    // compute lefm energy and j-integrals
    j_energy.d_lefmEnergyRate = vmag * d_matDeck_p->d_matData.d_Gc;
  }

  // old version of output file
  {
    // create file in first call
    if (!data->d_file) {
      std::string filename =
          d_outPreTag + d_currentData->d_tagFilename + ".csv";
      data->d_file = fopen(filename.c_str(), "w");

      // write header
      fprintf(data->d_file,
              "dt_out, vmag, E, v*Gc, Gc_experiment, Gc_theory\n");
    }

    // write data
    double Gcompute = 0.;
    auto vmag = ctip.d_v.length();
    if (util::compare::definitelyGreaterThan(vmag, 1.E-10))
      Gcompute = -j_energy.d_pdInternalWorkRate / vmag;
    fprintf(data->d_file, "%u, %4.6e, %4.6e, %4.6e, %4.6e, %4.6e\n", d_nOut,
            vmag, -j_energy.d_pdInternalWorkRate, j_energy.d_lefmEnergyRate,
            Gcompute, d_matDeck_p->d_matData.d_Gc);
  }  // old print format

  // new version of output file
  {
    // create file in first call
    if (!data->d_fileNew) {
      std::string filename =
          d_outPreTag + d_currentData->d_tagFilename + "_new" + ".csv";
      data->d_fileNew = fopen(filename.c_str(), "w");

      // write header
      fprintf(data->d_fileNew,
              "dt_out, V_mag, "
              "contour_strain_energy, "
              "contour_strain_energy_rate, "
              "contour_kinetic_energy_rate, "
              "contour_elastic_internal_work_rate, "
              "pd_internal_work, "
              "pd_internal_work_rate, "
              "lefm_energy_rate, "
              "pd_strain_energy_inside, "
              "kinetic_energy_inside, "
              "pd_fracture_energy\n");
    }

    // write data
    fprintf(data->d_fileNew,
            "%u, %4.6e, "
            "%4.6e, "
            "%4.6e, "
            "%4.6e, "
            "%4.6e, "
            "%4.6e, "
            "%4.6e, "
            "%4.6e, "
            "%4.6e, "
            "%4.6e, "
            "%4.6e\n",
            d_nOut, ctip.d_v.length(), j_energy.d_contourPdStrainEnergy,
            j_energy.d_contourPdStrainEnergyRate,
            j_energy.d_contourKineticEnergyRate,
            j_energy.d_contourElasticInternalWorkRate,
            j_energy.d_pdInternalWork, j_energy.d_pdInternalWorkRate,
            j_energy.d_lefmEnergyRate, j_energy.d_pdStrainEnergyInsideContour,
            j_energy.d_kineticEnergyInsideContour, j_energy.d_pdFractureEnergy);
  }
}

void tools::pp::Compute::computeJIntegralAngledCrack() {
  // TODO
  return;
}

//
// utility methods
//
void tools::pp::Compute::addUniqueToList(size_t i, std::vector<size_t> *list) {
  for (auto a : *list) {
    if (a == i) return;
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

  if (i_found < 0) safeExit("Error: Could not find the node.\n");

  return size_t(i_found);
}

void tools::pp::Compute::listElemsAndNodesInDomain(
    const std::pair<util::Point3, util::Point3> &cd, const double &tol,
    const double &tol_elem, std::vector<size_t> *nodes,
    std::vector<size_t> *elements, bool calc_in_ref) {
  // nodes list
  nodes->clear();
  for (size_t i = 0; i < d_dataManager_p->getMeshP()->getNumNodes(); i++) {
    auto x = d_dataManager_p->getMeshP()->getNode(i);
    if (!calc_in_ref) x += d_u[i];

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
  for (size_t e = 0; e < d_dataManager_p->getMeshP()->getNumElements(); e++) {
    auto ids = d_dataManager_p->getMeshP()->getElementConnectivity(e);
    bool add_e = false;

    // idea: if any of the nodes of this element is inside a zone then we add
    // the element to list.
    // Zone is given by difference of two rectangles enveloping contour

    for (auto i : ids) {
      if (add_e) break;

      auto x = d_dataManager_p->getMeshP()->getNode(i);
      if (!calc_in_ref) x += d_u[i];

      // check if node is in the bigger domain and not in smaller domain
      if (util::geometry::isPointInsideRectangle(
              x, cd.first.d_x - tol_elem, cd.second.d_x + tol_elem,
              cd.first.d_y - tol_elem, cd.second.d_y + tol_elem) &&
          !util::geometry::isPointInsideRectangle(
              x, cd.first.d_x + tol_elem, cd.second.d_x - tol_elem,
              cd.first.d_y + tol_elem, cd.second.d_y - tol_elem))
        add_e = true;
    }

    // add e to list
    if (add_e) addUniqueToList(e, elements);
  }
}

void tools::pp::Compute::decomposeSearchNodes(
    const std::pair<util::Point3, util::Point3> &cd, std::vector<size_t> *nodes,
    std::vector<size_t> *nodes_new) {
  std::vector<size_t> nodes_temp = *nodes;
  nodes_new->clear();
  nodes->clear();
  for (auto i : nodes_temp) {
    auto xi = d_dataManager_p->getMeshP()->getNode(i);
    if (util::geometry::isPointInsideRectangle(
            d_dataManager_p->getMeshP()->getNode(i), cd.first.d_x,
            cd.second.d_x, cd.first.d_y, cd.second.d_y))
      nodes->emplace_back(i);
    else
      nodes_new->emplace_back(i);
  }
}

bool tools::pp::Compute::triCheckAndInterpolateUV(
    const util::Point3 &p, util::Point3 &up, util::Point3 &vp,
    const std::vector<size_t> &ids, bool calc_in_ref, bool check_only) {
  // get triangle element object
  auto tri_quad = fe::TriElem(0);
  auto nodes = std::vector<util::Point3>{
      d_dataManager_p->getMeshP()->getNode(ids[0]) +
          (calc_in_ref ? util::Point3() : d_u[ids[0]]),
      d_dataManager_p->getMeshP()->getNode(ids[1]) +
          (calc_in_ref ? util::Point3() : d_u[ids[1]]),
      d_dataManager_p->getMeshP()->getNode(ids[2]) +
          (calc_in_ref ? util::Point3() : d_u[ids[2]])};
  double area = std::abs(tri_quad.elemSize(nodes));
  double sum_area = 0.;
  for (size_t a = 0; a < nodes.size(); a++)
    sum_area += std::abs(tri_quad.elemSize(std::vector<util::Point3>{
        p, nodes[a], nodes[(a == nodes.size() - 1 ? 0 : a + 1)]}));

  if (util::compare::definitelyGreaterThan(sum_area, area)) return false;

  if (check_only) return true;

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
                                       const std::vector<size_t> *elements,
                                       bool calc_in_ref) {
  // check if element data is available
  bool elem_interp = false;
  if (elements != nullptr)
    if (!elements->empty()) elem_interp = true;
  if (!elem_interp) {
    // use piecewise constant interpolation
    double dist = 1000.;
    long int loc_i = -1;
    // search for closest node
    for (auto i : *nodes) {
      auto xi = d_dataManager_p->getMeshP()->getNode(i);
      if (!calc_in_ref) xi += d_u[i];

      if (util::compare::definitelyLessThan(p.dist(xi), dist)) {
        dist = p.dist(xi);
        loc_i = i;
      }
    }
    if (loc_i == -1) {
      oss.str("");
      oss << "Error: Can not find node closer to point p = (" << p.d_x << ", "
          << p.d_y << ").\n";
      safeExit(oss.str());
    }
    up = d_u[loc_i];
    vp = d_v[loc_i];

    return;
  } else {
    for (auto e : *elements) {
      auto ids = d_dataManager_p->getMeshP()->getElementConnectivity(e);

      // cases
      if (d_dataManager_p->getMeshP()->getElementType() ==
          util::vtk_type_triangle) {
        if (triCheckAndInterpolateUV(p, up, vp, ids, calc_in_ref)) return;
      } else if (d_dataManager_p->getMeshP()->getElementType() ==
                 util::vtk_type_quad) {
        // check in triangle {v1, v2, v3}
        if (triCheckAndInterpolateUV(
                p, up, vp, std::vector<size_t>{ids[0], ids[1], ids[2]},
                calc_in_ref))
          return;

        // check in triangle {v1, v3, v4}
        if (triCheckAndInterpolateUV(
                p, up, vp, std::vector<size_t>{ids[0], ids[2], ids[3]},
                calc_in_ref))
          return;
      }
    }

    // issue error since element containing the point is not found
    oss.str("");
    oss << "Error: Can not find element for point p = (" << p.d_x << ", "
        << p.d_y << ") for interpolation.\n";
    oss << "Num elems = " << elements->size() << "\n";
    // for debug
    std::vector<util::Point3> enodes;
    for (auto e : *elements) {
      for (auto n : d_dataManager_p->getMeshP()->getElementConnectivity(e))
        enodes.emplace_back(d_dataManager_p->getMeshP()->getNode(n) +
                            (calc_in_ref ? util::Point3() : d_u[n]));
    }
    enodes.emplace_back(p);
    std::vector<size_t> etags(enodes.size(), 1);
    etags[etags.size() - 1] = 100;
    auto writer1 =
        rw::writer::Writer(d_outPreTag + d_currentData->d_tagFilename +
                           "_debug_nodes_" + std::to_string(d_nOut));
    writer1.appendNodes(&enodes);
    writer1.appendPointData("Tag", &etags);
    writer1.close();

    safeExit(oss.str());
  }
}

size_t tools::pp::Compute::interpolateUVNodes(const util::Point3 &p,
                                              util::Point3 &up,
                                              util::Point3 &vp,
                                              const std::vector<size_t> *nodes,
                                              bool calc_in_ref) {
  // use piecewise constant interpolation
  double dist = 1000.;
  long int loc_i = -1;
  // search for closest node
  for (auto i : *nodes) {
    auto xi = d_dataManager_p->getMeshP()->getNode(i);
    if (!calc_in_ref) xi += d_u[i];

    if (util::compare::definitelyLessThan(p.dist(xi), dist)) {
      dist = p.dist(xi);
      loc_i = i;
    }
  }
  if (loc_i == -1) {
    oss.str("");
    oss << "Error: Can not find node closer to point p = (" << p.d_x << ", "
        << p.d_y << ").\n";
    safeExit(oss.str());
  }
  up = d_u[loc_i];
  vp = d_v[loc_i];

  return loc_i;
}

void tools::pp::Compute::getContourContribJInt(
    const util::Point3 &p, const std::vector<size_t> *nodes,
    const std::vector<size_t> *elements, const util::Point3 &normal,
    double &pd_energy, double &kinetic_energy, double &elastic_energy,
    util::Point3 &dot_u, bool calc_in_ref,
    const tools::pp::CrackTipData &ctip) {
  const auto crack_orient = d_currentData->d_computeJInt_p->d_crackOrient;

  // reset data
  pd_energy = 0.;
  kinetic_energy = 0.;

  // get displacement and velocity at point
  auto uq = util::Point3();
  auto vq = util::Point3();
  // interpolateUV(p, uq, vq, nodes, elements, calc_in_ref);
  auto node_p = interpolateUVNodes(p, uq, vq, nodes, calc_in_ref);

  // add kinetic energy
  kinetic_energy = 0.5 * d_material_p->getDensity() * vq.dot(vq);

  // add peridynamic energy
  double loc_pd_energy = 0.;
  if (!d_material_p->isStateActive()) {
    // upper and lower bound for volume correction
    auto h = d_dataManager_p->getMeshP()->getMeshSize();
    auto check_up = d_dataManager_p->getModelDeckP()->d_horizon + 0.5 * h;
    auto check_low = d_dataManager_p->getModelDeckP()->d_horizon - 0.5 * h;

    auto i_neighs = d_dataManager_p->getNeighborP()->getNeighbors(node_p);
    for (size_t j = 0; j < i_neighs.size(); j++) {
      auto j_id = i_neighs[j];

      auto xj = d_dataManager_p->getMeshP()->getNode(j_id);
      if (util::compare::definitelyGreaterThan(xj.dist(p),
                                               d_dataManager_p->getModelDeckP()->d_horizon) ||
          util::compare::definitelyLessThan(xj.dist(p), 1.0E-10))
        continue;

      // get displacement and strain
      auto uj = d_u[j_id];
      auto rjq = xj.dist(p);
      auto Sjq = d_material_p->getS(xj - p, uj - uq);

      // get corrected volume of node j
      auto volj = d_dataManager_p->getMeshP()->getNodalVolume(j_id);
      if (util::compare::definitelyGreaterThan(rjq, check_low))
        volj *= (check_up - rjq) / h;

      bool fracture_state = false;
      auto ef = d_material_p->getBondEF(node_p, j);

      // add contribution to energy
      // Debug
      //  Add contribution to strain energy only if the bond is strain above
      //  critical strain
      //  This does not give good result
      if (false) {
        double sr = std::abs(Sjq) / this->d_material_p->getSc(rjq);
        if (util::compare::definitelyGreaterThan(sr, 1.))
          loc_pd_energy += ef.second * volj;
      }

      // Debug
      //  Add contribution to strain energy only if the bond intersects the
      //  crack line
      if (false) {
        if (doesBondIntersectCrack(xj, p, ctip, crack_orient))
          loc_pd_energy += ef.second * volj;
      }

      loc_pd_energy += ef.second * volj;
    }  // loop over nodes for pd energy density

  } else {
    auto thetap = d_thetaX[node_p];
    auto mp = d_mX[node_p];

    // upper and lower bound for volume correction
    auto h = d_dataManager_p->getMeshP()->getMeshSize();
    auto check_up = d_dataManager_p->getModelDeckP()->d_horizon + 0.5 * h;
    auto check_low = d_dataManager_p->getModelDeckP()->d_horizon - 0.5 * h;

    auto i_neighs = d_dataManager_p->getNeighborP()->getNeighbors(node_p);
    for (size_t j = 0; j < i_neighs.size(); j++) {
      auto j_id = i_neighs[j];

      auto xj = d_dataManager_p->getMeshP()->getNode(j_id);

      if (util::compare::definitelyGreaterThan(xj.dist(p),
                                               d_dataManager_p->getModelDeckP()->d_horizon) ||
          util::compare::definitelyLessThan(xj.dist(p), 1.0E-10))
        continue;

      // get displacement and strain
      auto uj = d_u[j_id];
      auto thetaj = d_thetaX[j_id];
      auto mj = d_mX[j_id];
      auto rjq = xj.dist(p);
      auto Sjq = d_material_p->getS(xj - p, uj - uq);

      // get corrected volume of node j
      auto volj = d_dataManager_p->getMeshP()->getNodalVolume(j_id);
      if (util::compare::definitelyGreaterThan(rjq, check_low))
        volj *= (check_up - rjq) / h;

      // auto fs = this->d_fracture_p->getBondState(node_p, j);
      auto ef_i = this->d_material_p->getBondEF(node_p, j);
      auto ef_j = this->d_material_p->getBondEF(node_p, j);

      // add contribution to energy
      loc_pd_energy += (ef_i.second + ef_j.second) * volj;
    }  // loop over nodes for pd energy density
  }

  pd_energy = loc_pd_energy;

  // get elastic work done
  dot_u = vq;

  // if there is no element-node connectivity data, we can not compute
  // elastic work done, so return
  if (d_dataManager_p->getMeshP()->getNumElements() == 0) return;

  // get Quadrature
  fe::BaseElem *quad;
  if (d_dataManager_p->getMeshP()->getElementType() == util::vtk_type_triangle)
    quad = new fe::TriElem(1);
  else if (d_dataManager_p->getMeshP()->getElementType() == util::vtk_type_quad)
    quad = new fe::QuadElem(1);
  else {
    return;
  }

  // loop over elements and find the element which has quadrature point
  // closest to the point
  size_t e_found = 0;
  auto dist = DBL_MAX;
  for (const auto e : *elements) {
    // get ids of nodes of element, coordinate of nodes, 1st order
    // quad data, and first quad data
    auto id_nds = d_dataManager_p->getMeshP()->getElementConnectivity(e);
    auto nds = d_dataManager_p->getMeshP()->getElementConnectivityNodes(e);

    if (!calc_in_ref) {
      for (size_t i = 0; i < id_nds.size(); i++) nds[i] += d_u[id_nds[i]];
    }

    auto qds = quad->getQuadDatas(nds);
    auto dx = qds[0].d_p - p;

    if (dist < dx.length()) {
      dist = dx.length();
      e_found = e;
    }
  }

  // now compute strain and stress at the quadrature point of found element
  {
    auto ssn = util::SymMatrix3();
    auto sss = util::SymMatrix3();

    // get ids of nodes of element, coordinate of nodes, 1st order
    // quad data, and first quad data
    auto id_nds = d_dataManager_p->getMeshP()->getElementConnectivity(e_found);
    auto nds =
        d_dataManager_p->getMeshP()->getElementConnectivityNodes(e_found);
    if (!calc_in_ref) {
      for (size_t i = 0; i < id_nds.size(); i++) nds[i] += d_u[id_nds[i]];
    }

    auto qds = quad->getQuadDatas(nds);
    auto qd0 = qds[0];

    // compute strain in xy plane
    for (size_t i = 0; i < id_nds.size(); i++) {
      auto id = id_nds[i];
      auto ui = d_u[id];

      ssn(0, 0) += ui.d_x * qd0.d_derShapes[i][0];
      ssn(1, 1) += ui.d_y * qd0.d_derShapes[i][1];
      ssn(0, 1) += 0.5 * ui.d_x * qd0.d_derShapes[i][1] +
                   0.5 * ui.d_y * qd0.d_derShapes[i][0];
    }

    if (!d_matDeck_p->d_isPlaneStrain)
      ssn(2, 2) = -d_matDeck_p->d_matData.d_nu * (ssn(0, 0) + ssn(1, 1)) /
                  (1. - d_matDeck_p->d_matData.d_nu);

    // compute stress
    auto trace = ssn(0, 0) + ssn(1, 1) + ssn(2, 2);
    sss(0, 0) = d_matDeck_p->d_matData.d_lambda * trace * 1. +
                2. * d_matDeck_p->d_matData.d_mu * ssn(0, 0);

    sss(0, 1) = 2. * d_matDeck_p->d_matData.d_mu * ssn(0, 1);

    sss(1, 1) = d_matDeck_p->d_matData.d_lambda * trace * 1. +
                2. * d_matDeck_p->d_matData.d_mu * ssn(1, 1);

    if (d_matDeck_p->d_isPlaneStrain)
      sss(2, 2) = d_matDeck_p->d_matData.d_nu * (sss(0, 0) + sss(1, 1));

    // static int debug = -1;
    // if (debug < 0) {
    //   std::cout << "Lambda = " << d_matDeck_p->d_matData.d_lambda
    //             << ", G = " << d_matDeck_p->d_matData.d_mu << "\n";
    //   debug = 0;
    // }

    // compute Stress * dot_u * normal
    util::Point3 stress_dot_v = sss.dot(dot_u);
    elastic_energy = stress_dot_v * normal;
  }
}

void tools::pp::Compute::updateCrack(const double &time,
                                     const std::vector<double> *Z) {
  auto compute_data = d_currentData->d_findCrackTip_p;

  // loop over crack lines
  for (auto &crack : compute_data->d_cracks) {
    if (crack.d_o == 0) continue;

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
      auto pnew = findTipInRects(crack, max_Z, rects_b, nodes_b, Z_b, Z, false);
      addNewCrackTip(crack, pnew, time, false);
    }
  }  // loop over cracks
}

void tools::pp::Compute::crackOutput() {
  auto compute_data = d_currentData->d_findCrackTip_p;

  // stop when update exceeds upper bound
  int up_bound = 10000;
  if (compute_data->d_updateCount >= up_bound) {
    std::cout << "Warning: Number of times crack data output requested "
                 "exceeds the upper limit 10000.\n";
    if (compute_data->d_filet) fclose(compute_data->d_filet);
    if (compute_data->d_fileb) fclose(compute_data->d_fileb);
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
    fprintf(compute_data->d_filet,
            "'crack_id', 'dt_out', 'x', 'y', 'vx', "
            "'vy'\n");
    fprintf(compute_data->d_fileb,
            "'crack_id', 'dt_out', 'x', 'y', 'vx', "
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

  auto h = d_dataManager_p->getMeshP()->getMeshSize();
  auto horizon = d_dataManager_p->getModelDeckP()->d_horizon;
  auto bbox = d_dataManager_p->getMeshP()->getBoundingBox();
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
      if (util::compare::definitelyLessThan(pt.d_x + Nt * seq_size,
                                            bbox.second[0]))
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
      if (util::compare::definitelyLessThan(pt.d_y + Nt * seq_size,
                                            bbox.second[1]))
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
        rects_t.emplace_back(std::make_pair(
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
        rects_b.emplace_back(std::make_pair(
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

  // create a bigger rectangle which has all the rectangles inside
  // we use this bigger rectangle to filter out the nodes not on it (to speed
  // up)
  // Not used currently
  std::pair<util::Point3, util::Point3> glob_rect;
  if (crack.d_o == 1)
    glob_rect = std::make_pair(
        util::Point3(bbox.first[0], crack.d_initPb.d_y - 2. * horizon, 0.),
        util::Point3(bbox.second[0], crack.d_initPb.d_y + 2. * horizon, 0.));
  else if (crack.d_o == -1)
    glob_rect = std::make_pair(
        util::Point3(crack.d_initPb.d_x - 2. * horizon, bbox.first[1], 0.),
        util::Point3(crack.d_initPb.d_x + 2. * horizon, bbox.second[1], 0.));

  // resize nodes list and damage list
  nodes_t.resize(rects_t.size());
  nodes_b.resize(rects_b.size());
  Z_t.resize(rects_t.size());
  Z_b.resize(rects_b.size());

  for (size_t i = 0; i < d_dataManager_p->getMeshP()->getNumNodes(); i++) {
    auto xi = d_dataManager_p->getMeshP()->getNode(i);
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
              break;
            }
          }
          if (!found) {
            nodes_t[r].emplace_back(i);
            Z_t[r].emplace_back(damage);
          }
        }
      }
    }  // top point
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
              break;
            }
          }
          if (!found) {
            nodes_b[r].emplace_back(i);
            Z_b[r].emplace_back(damage);
          }
        }
      }
    }  // bottom point
  }    // loop over nodes

  std::cout << "Max Z allowed = " << compute_data->d_maxZAllowed
            << ", min Z allowed = " << compute_data->d_minZAllowed << "\n";
}

void tools::pp::Compute::addNewCrackTip(inp::EdgeCrack &crack,
                                        util::Point3 pnew, double time,
                                        bool is_top) {
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
    crack.d_vb = util::Point3(diff.d_x / delta_t, diff.d_y / delta_t,
                              diff.d_z / delta_t);
    compute_data->d_timeb = time;
  }
}

util::Point3 tools::pp::Compute::findTipInRects(
    inp::EdgeCrack &crack, const double &max_Z,
    const std::vector<std::pair<util::Point3, util::Point3>> &rects,
    const std::vector<std::vector<size_t>> &nodes,
    const std::vector<std::vector<double>> &Zs, const std::vector<double> *Z,
    bool is_top) {
  const bool calc_in_ref = d_currentData->d_calculateInRefConfig;

  if (is_top)
    std::cout << "    Tip information: Top     \n";
  else
    std::cout << "    Tip information: Bottom     \n";
  std::cout << "    max_Z = " << max_Z << "\n";
  std::cout << "    size rects = " << rects.size() << "\n";
  // find size of rectangle with atleast one node
  size_t rect_size_with_node = 0;
  for (size_t r = 0; r < rects.size(); r++) {
    if (nodes[r].size() > 0) rect_size_with_node++;
  }
  std::cout << "    size valid rects = " << rect_size_with_node << "\n";

  // get old tip
  auto pold = crack.d_pt;
  if (!is_top) pold = crack.d_pb;

  // if there are no rectangles with desired nodes, return old tip
  if (rect_size_with_node == 0) return pold;

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
    if (!Zs[r].empty()) {
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
    if (r > 2) continue;

    auto mz = sortZ[r];
    auto y0 = d_dataManager_p->getMeshP()->getNode(mz.d_i);
    if (!calc_in_ref) y0 += d_u[mz.d_i];
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
      if (x == mz.d_i) continue;

      auto y = d_dataManager_p->getMeshP()->getNode(x);
      if (!calc_in_ref) y += d_u[x];
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
  auto p_choices =
      std::vector<util::Point3>{crack.d_initPt, crack.d_oldPt, crack.d_pt};
  if (!is_top)
    p_choices =
        std::vector<util::Point3>{crack.d_initPb, crack.d_oldPb, crack.d_pb};
  for (size_t r = 0; r < sortZ.size(); r++) {
    // only consider first and second rectangle in sorted rectangle list
    if (r > 2) break;

    // get data for this rectangle
    auto mz = sortZ[r];

    // do not consider node with large damage difference to the minimum
    // damage in mz.d_Z
    // diff_z will be further reduced to find node with closest damage to mz.d_Z
    double diff_z = 1.0E-02;

    // loop over crack line choices
    for (auto p_search : p_choices) {
      // check if we already found symmetrically opposite node for crack line
      // before this crack line
      if (sym_rect[r] >= 0) break;

      for (auto x : nodes[mz.d_r]) {
        auto y = d_dataManager_p->getMeshP()->getNode(x);
        auto y0 = d_dataManager_p->getMeshP()->getNode(mz.d_i);
        if (!calc_in_ref) {
          y += d_u[x];
          y0 += d_u[mz.d_i];
        }
        if (util::compare::definitelyLessThan(std::abs((*Z)[x] - mz.d_Z),
                                              diff_z) &&
            x != mz.d_i) {
          // find if it is symmetrically opposite to the node i
          if (crack.d_o == -1) {
            // compare x coordinate
            auto d0 = p_search.d_x - y0.d_x;
            auto d = p_search.d_x - y.d_x;

            // using not(definitely greater than) to allow value to be 0 as well
            if (!util::compare::definitelyGreaterThan(d * d0, 0.)) {
              sym_rect[r] = x;
              diff_z = std::abs((*Z)[x] - mz.d_Z);
            }
          } else if (crack.d_o == 1) {
            // compare y coordinate
            auto d0 = p_search.d_y - y0.d_y;
            auto d = p_search.d_y - y.d_y;

            // using not(definitely greater than) to allow value to be 0 as well
            if (!util::compare::definitelyGreaterThan(d * d0, 0.)) {
              sym_rect[r] = x;
              diff_z = std::abs((*Z)[x] - mz.d_Z);
            }
          }
        }
      }  // loop over nodes within rectangle
    }
  }  // loop over rectangle

  // output data for debugging
  std::cout << "    sortZ (r, i, Z) = {";
  for (size_t i = 0; i < sortZ.size(); i++) {
    std::cout << "(" << sortZ[i].d_r << "," << sortZ[i].d_i << ","
              << sortZ[i].d_Z << ")";
    if (i < sortZ.size() - 1) std::cout << ", ";
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
    auto y1 = d_dataManager_p->getMeshP()->getNode(i1);
    size_t i2 = sym_rect[0];
    auto y2 = d_dataManager_p->getMeshP()->getNode(i2);

    if (!calc_in_ref) {
      y1 += d_u[i1];
      y2 += d_u[i2];
    }

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
  }  // rect 0 is symmetric
  else {
    if (sym_rect[1] >= 0) {
      size_t i1 = sortZ[1].d_i;
      auto y1 = d_dataManager_p->getMeshP()->getNode(i1);
      size_t i2 = sym_rect[1];
      auto y2 = d_dataManager_p->getMeshP()->getNode(i2);

      if (!calc_in_ref) {
        y1 += d_u[i1];
        y2 += d_u[i2];
      }

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
    }  // rect 1 is symmetric
    else {
      // use average between best node in rectangle and old tip
      size_t i1 = sortZ[0].d_i;
      auto y1 = d_dataManager_p->getMeshP()->getNode(i1);
      auto y2 = pold;

      if (!calc_in_ref) y1 += d_u[i1];

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
    }  // rect 1 is not symmetric
  }    // rect 0 is not symmetric

  // output data for debugging
  std::cout << "    Old tip = (" << pold.d_x << ", " << pold.d_y << "), "
            << "New tip = (" << pnew.d_x << ", " << pnew.d_y << ").\n";

  return pnew;
}
