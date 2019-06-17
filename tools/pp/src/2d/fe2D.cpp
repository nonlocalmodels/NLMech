// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include "fe2D.h"
#include "external/csv.h"           // csv reader
#include "fe/mesh.h"                // definition of Mesh
#include "fe/quadrature.h"          //
#include "fe/quadElem.h"            // definition of QuadElem
#include "fe/triElem.h"             // definition of TriElem
#include "fe/lineElem.h"            // definition of LineElem
#include "geometry/fracture.h"      // definition of Fracture
#include "geometry/neighbor.h"      // definition of Neighbor
#include "inp/decks/materialDeck.h" // definition of MaterialDeck
#include "inp/decks/modelDeck.h"    // definition of ModelDeck
#include "inp/decks/outputDeck.h"   // definition of OutputDeck
#include "inp/input.h"              // definition of Input
#include "inp/policy.h"             // definition of Policy
#include "material/pdMaterial.h"    // definition of Material
#include "rw/reader.h"              // definition of readVtuFileRestart
#include "rw/writer.h"              // definition of VtkWriterInterface
#include "util/compare.h"           // compare real numbers
#include "util/feElementDefs.h"     // definition of fe element type
#include "util/matrix.h"            // definition of SymMatrix3
#include "util/point.h"             // definition of Point3
#include "util/utilGeom.h"          // definition of isPointInsideRectangle
#include <cmath>
#include <geometry/fracture.h>
#include <hpx/include/parallel_algorithm.hpp>
#include <iostream>
#include <yaml-cpp/yaml.h> // YAML reader

namespace {

struct CrackTipData {
  size_t d_n;
  util::Point3 d_p;
  util::Point3 d_v;

  CrackTipData() : d_n(0), d_p(util::Point3()), d_v(util::Point3()){};
  CrackTipData(size_t n, util::Point3 p, util::Point3 v)
      : d_n(n), d_p(p), d_v(v){};
};

struct InstructionData {
  std::string d_tagFilename;

  int d_start;
  int d_end;

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

  bool d_magStrainTensor;
  std::string d_magStrainComp;
  std::vector<std::pair<size_t, double>> d_markMagStrainCells;

  bool d_symmetrizeV;
  bool d_combineMarkV;
  std::string d_symmAxis;
  double d_symmLine;

  bool d_crackTip;
  bool d_crackSameDtOut;

  bool d_compJIntegral;
  int d_crackOrient;
  std::string d_crackTipFile;
  std::vector<double> d_contourFactor;
  std::vector<CrackTipData> d_crackTipData;
  size_t d_startJIntegral;
  size_t d_endJIntegral;

  InstructionData()
      : d_start(-1), d_end(-1), d_scaleUOut(false), d_scaleU(1.),
        d_damageAtNodes(false), d_outOnlyNodes(true), d_removeElements(false),
        d_markVAsZero(false), d_markVInRectGiven(false),
        d_markVRect(std::make_pair(util::Point3(), util::Point3())),
        d_markVPtsAreInCurrentConfig(false), d_computeStrain(false),
        d_magStrainTensor(false), d_symmetrizeV(false), d_combineMarkV(false),
        d_symmLine(0.), d_crackTip(false), d_crackSameDtOut(true),
        d_compJIntegral(false), d_crackOrient(0), d_startJIntegral(0),
        d_endJIntegral(0){};
};

void readInputFile(YAML::Node config, const std::string &set,
                   InstructionData *data) {

  if (config["Compute"][set]["Tag_Filename"]) {
    data->d_tagFilename =
        config["Compute"][set]["Tag_Filename"].as<std::string>();
  } else {
    // we work with default filename of form
    // pp_set_1_time_step_0.vtu
    data->d_tagFilename = set;
  }

  //
  // Output only nodes
  //
  if (config["Compute"][set]["Output_Only_Nodes"])
    data->d_outOnlyNodes =
        config["Compute"][set]["Output_Only_Nodes"].as<bool>();

  //
  // check if start and end time step are specified
  //
  if (config["Compute"][set]["Dt_Start"])
    data->d_start = config["Compute"][set]["Dt_Start"].as<int>();

  if (config["Compute"][set]["Dt_End"])
    data->d_end = config["Compute"][set]["Dt_End"].as<int>();

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
  if (config["Compute"][set]["Strain_Stress"])
    data->d_computeStrain = config["Compute"][set]["Strain_Stress"].as<bool>();

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
    data->d_crackTip = true;
    if (config["Compute"][set]["Crack_Tip"]["Same_Dt_Out"])
      data->d_crackSameDtOut =
          config["Compute"][set]["Crack_Tip"]["Same_Dt_Out"].as<bool>();
  }

  //
  // J integral
  //
  if (config["Compute"][set]["J_Integral"]) {
    auto e = config["Compute"][set]["J_Integral"];
    data->d_compJIntegral = true;
    if (e["Crack_Orient"])
      data->d_crackOrient = e["Crack_Orient"].as<int>();
    else {
      std::cerr << "Error: Crack orientation is not provided.\n";
      exit(1);
    }
    if (e["Crack_Tip_File"])
      data->d_crackTipFile = e["Crack_Tip_File"].as<std::string>();
    else {
      std::cerr << "Error: Crack tip information filename is not provided.\n";
      exit(1);
    }
    if (e["Contour_Size"]) {
      for (auto f : e["Contour_Size"])
        data->d_contourFactor.push_back(f.as<double>());

      if (data->d_contourFactor.size() == 1)
        data->d_contourFactor.push_back(data->d_contourFactor[0]);
    } else {
      std::cerr << "Error: Factors to create contour for J integral not "
                   "provided.\n";
      exit(1);
    }
  }
}

size_t findNode(const util::Point3 &x, const std::vector<util::Point3> *nodes,
                const std::vector<util::Point3> *u = nullptr) {

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

void computeDamage(inp::Input *deck, inp::ModelDeck *model_deck,
                   material::pd::Material *material, fe::Mesh *mesh,
                   const std::vector<util::Point3> *u, std::vector<float> *Z) {
  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      mesh->getNumNodes(),
      [Z, mesh, material, model_deck, u](boost::uint64_t i) {
        auto xi = mesh->getNode(i);
        for (size_t j = 0; j < mesh->getNumNodes(); j++) {
          if (util::compare::definitelyGreaterThan(xi.dist(mesh->getNode(j)),
                                                   model_deck->d_horizon) ||
              j == i)
            continue;

          auto xj = mesh->getNode(j);
          if (util::compare::definitelyGreaterThan(xj.dist(xi), 1.0E-10)) {
            auto Sr = std::abs(material->getS(xj - xi, (*u)[j] - (*u)[i])) /
                      material->getSc(xj.dist(xi));

            if (util::compare::definitelyLessThan((*Z)[i], Sr))
              (*Z)[i] = Sr;
          }
        } // loop over neighbors
      }); // parallel loop over nodes
  f.get();
}

void readCrackTipData(const std::string &filename,
                      std::vector<CrackTipData> *data) {
  // expected format of file:
  // <output step>, <tip x>, <tip y>, <tip vx>, <tip vy>
  io::CSVReader<5> in(filename);
  in.read_header(io::ignore_extra_column, "id", "x", "y", "z", "volume");

  double px, py, vx, vy;
  int n;
  while (in.read_row(n, px, py, vx, vy)) {
    data->emplace_back(CrackTipData(size_t(n), util::Point3(px, py, 0.),
                                    util::Point3(vx, vy, 0.)));
  }
}

void interpolateUV(util::Point3 p, util::Point3 &uq, util::Point3 &vq,
                   std::vector<util::Point3> *u, std::vector<util::Point3> *v,
                   fe::Mesh *mesh) {
  // check if element data is available
  if (mesh->getNumElements() == 0) {
    // use piecewise constant interpolation
    long int loc_i = -1;
    double dist = 1000.0;
    for (size_t i=0; i<mesh->getNumNodes(); i++) {
      if (util::compare::definitelyLessThan(p.dist(mesh->getNode(i)), dist)) {
        dist = p.dist(mesh->getNode(i));
        loc_i = i;
      }
    }

    if (loc_i == -1) {
      std::cerr << "Error: Can not locate node near to point p = (" << p.d_x <<
      ","<<p.d_y<<").\n";
      exit(1);
    }

    uq = (*u)[loc_i];
    vq = (*v)[loc_i];
  }
  else {
    // find element containing
  }
}

void computeJIntegral(const size_t &out, const InstructionData &data,
                      inp::ModelDeck *model_deck, inp::MaterialDeck *mat_deck,
                      material::pd::Material *material, fe::Mesh *mesh,
                      geometry::Neighbor * neighbor_list,
                      std::vector<util::Point3> *u,
                      std::vector<util::Point3> *v, double &energy) {

  if (out < data.d_startJIntegral || out > data.d_endJIntegral)
    return;

  auto ctip = data.d_crackTipData[out - data.d_startJIntegral];
  //
  // create rectangle domain of specified size at the tip of crack
  //
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
  // 1. Contour is formed by lines A-B, B-C, C-D, D-A
  //
  // 2. Discretization scheme for integral over contour
  // a. Discretize edges with uniform number of points and use linear line
  // element for the interpolation
  // b. Quadrature approximation of second order for integration over line
  // elements
  // c. Mesh size of uniform discretization is fixed to half of the mesh size
  // of the simulation
  //
  std::pair<util::Point3, util::Point3> cd(
      std::make_pair(util::Point3(), util::Point3()));
  if (data.d_crackOrient == -1) {
    // vertical crack
    cd.first = util::Point3(ctip.d_p.d_x - 0.5 * data.d_contourFactor[0] *
                                               model_deck->d_horizon,
                            ctip.d_p.d_y, 0.);
    cd.second = util::Point3(
        ctip.d_p.d_x + 0.5 * data.d_contourFactor[0] * model_deck->d_horizon,
        ctip.d_p.d_y + data.d_contourFactor[1] * model_deck->d_horizon, 0.);
  } else if (data.d_crackOrient == 1) {
    // horizontal crack
    cd.first = util::Point3(ctip.d_p.d_x,
                            ctip.d_p.d_y - 0.5 * data.d_contourFactor[1] *
                                               model_deck->d_horizon,
                            0.);
    cd.second = util::Point3(
        ctip.d_p.d_x + data.d_contourFactor[0] * model_deck->d_horizon,
        ctip.d_p.d_y + 0.5 * data.d_contourFactor[1] * model_deck->d_horizon,
        0.);
  }

  // create second order quadrature class for 1-d line element
  auto quad = fe::LineElem(2);

  // Discretize edge A - B and C - D
  auto h = mesh->getMeshSize();
  size_t N = (cd.second.d_x - cd.first.d_x) / h;
  if (util::compare::definitelyLessThan(cd.fist.d_x + double(N) * h, cd
  .second.d_x))
    N++;
  for (size_t I=0; I<N; I++) {

    // line element
    auto x1 = cd.first.d_x + double(I) * h;
    auto x2 = cd.first.d_x + double(I+1) * h;
    if (I == N - 1)
      x2 = cd.second.d_x;

    // get quadrature points
    auto qds = quad.getQuadPoints(std::vector<util::Point3>{
        util::Point3(x1, 0., 0.), util::Point3(x2, 0., 0.)});

    // Edge A - B
    auto y = cd.first.d_y;

    // get normal
    auto n = util::Point3(0.,-1., 0.);

    // loop over quad points
    for (auto qd : qds) {
      // modify the y coordinate of quad point
      qd.d_p.d_y = y;

      // either interpolate or use piecewise constant approximation to get the
      // displacement and velocity at the quadrature point
      auto uq = util::Point3();
      auto vq = util::Point3();
      interpolateUV(qd.d_p, uq, vq, u, v, mesh);
    }
  }


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
  auto *deck = new inp::Input(sim_filename);
  auto model_deck = deck->getModelDeck();
  auto output_deck = deck->getOutputDeck();
  auto fracture_deck = deck->getFractureDeck();

  // get policy deck
  auto policy = inp::Policy::getInstance(deck->getPolicyDeck());

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
  std::string out_filename;
  if (config["Output"]["Filename"])
    out_filename =
        out_path + "/" + config["Output"]["Filename"].as<std::string>();
  else
    out_filename = out_path + "/pp";

  int start_file = 1;
  int end_file = model_deck->d_Nt / output_deck->d_dtOut;

  // create mesh
  std::cout << "PP_fe2D: Creating mesh.\n";
  auto *mesh = new fe::Mesh(deck->getMeshDeck());

  // material deck and material
  auto *material = new material::pd::Material(
      deck->getMaterialDeck(), model_deck->d_dim, model_deck->d_horizon);
  auto mat_deck = material->getMaterialDeck();
  // if material deck does not have valid material properties,
  // search the properties in the pp input file
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

  // to hold compute input data
  auto num_compute = config["Compute"]["Sets"].as<size_t>();
  std::vector<InstructionData> compute_data;
  for (size_t c = 0; c < num_compute; c++) {
    std::string set = "Set_" + std::to_string(c + 1);
    auto data = InstructionData();
    readInputFile(config, set, &data);
    if (data.d_start == -1)
      data.d_start = int(start_file);
    if (data.d_end == -1)
      data.d_end = int(end_file);

    compute_data.emplace_back(data);
  }

  // vector of null pointer Fracture class
  std::vector<geometry::Fracture *> fractures(num_compute, nullptr);
  std::vector<inp::FractureDeck> fracture_decks(num_compute, *fracture_deck);

  // read crack tip data for J integral calculation
  for (auto &data : compute_data) {
    if (data.d_compJIntegral) {
      readCrackTipData(data.d_crackTipFile, &(data.d_crackTipData));
      data.d_startJIntegral = data.d_crackTipData[0].d_n;
      data.d_endJIntegral =
          data.d_crackTipData[data.d_crackTipData.size() - 1].d_n;
    }
  }

  //
  // compute neighbor list (if required)
  //
  bool compute_neighbor_list = false;
  for (const auto& d : compute_data)
    compute_neighbor_list = d.d_compJIntegral;

  geometry::Neighbor * neighbor_list = nullptr;
  if (compute_neighbor_list)
    neighbor_list = new geometry::Neighbor(model_deck->d_horizon,
                                           deck->getNeighborDeck(),
                                           mesh->getNodesP());


    //
  // loop over output files
  //
  for (int out = start_file; out <= end_file; out++) {
    std::cout << "PP_fe2D: Processing output file = " << out << "\n";

    // append path to filename
    auto sim_results_file = source_path + "/" + filename_to_read + "_" +
                            std::to_string(out) + ".vtu";

    // see if we need to read the data
    bool read_file = false;
    for (const auto &data : compute_data)
      read_file = out >= data.d_start && out <= data.d_end;

    // just read one file and do operations and move to the next compute set
    std::vector<util::Point3> u;
    std::vector<util::Point3> v;

    // get displacement and velocity from the simulation
    if (read_file)
      rw::reader::readVtuFileRestart(sim_results_file, &u, &v);

    //    // output file for debug
    //    // open a output vtu file
    //    auto writer_debug = rw::writer::VtkWriterInterface(out_filename +
    //        "_debug" + "_" + std::to_string(out));
    //    writer_debug.appendNodes(mesh->getNodesP(), &u);
    //    writer_debug.appendPointData("Displacement", &u);
    //    writer_debug.appendPointData("Velocity", &v);
    //    writer_debug.close();

    // loop over compute sets and do as instructed in input file
    for (size_t c = 0; c < num_compute; c++) {
      auto data = compute_data[c];

      // continue if in the specified bound
      if (out < data.d_start || out > data.d_end)
        continue;

      std::cout << "  PP_fe2D: Processing compute set = " << c + 1 << "\n";

      // open a output vtu file
      std::string compute_filename =
          out_filename + "_" + data.d_tagFilename + "_" + std::to_string(out);
      rw::writer::VtkWriterInterface *writer = nullptr;

      //
      // operation : Scale displacement and write to output mesh
      //
      // Steps: Create displacement vector and scale it and write to mesh
      //
      if (data.d_scaleUOut) {

        std::vector<util::Point3> u_temp(mesh->getNumNodes(), util::Point3());
        //        auto f = hpx::parallel::for_loop(
        //            hpx::parallel::execution::par(hpx::parallel::execution::task),
        //            0, mesh->getNumNodes(), [&u_temp, data, u](boost::uint64_t
        //            i) {
        for (size_t i = 0; i < mesh->getNumNodes(); i++)
          u_temp[i] =
              util::Point3(data.d_scaleU * u[i].d_x, data.d_scaleU * u[i].d_y,
                           data.d_scaleU * u[i].d_z);
        //            });
        //        f.get();

        if (!writer) {
          writer = new rw::writer::VtkWriterInterface(compute_filename);
          // append mesh (check if only nodes need to be written)
          if (data.d_outOnlyNodes)
            writer->appendNodes(mesh->getNodesP(), &u_temp);
          else
            writer->appendMesh(mesh->getNodesP(), mesh->getElementType(),
                               mesh->getElementConnectivitiesP(), &u_temp);
        }

        // append original displacement
        writer->appendPointData("Displacement", &u);

        // append velocity
        writer->appendPointData("Velocity", &v);
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

        if (data.d_markVInRectGiven) {
          auto f = hpx::parallel::for_loop(
              hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
              mesh->getNumNodes(), [&v_mark, data, mesh](boost::uint64_t i) {
                if (util::geometry::isPointInsideRectangle(
                        mesh->getNode(i), data.d_markVRect.first.d_x,
                        data.d_markVRect.second.d_x, data.d_markVRect.first.d_y,
                        data.d_markVRect.second.d_y))
                  v_mark[i] = util::Point3();
              });
          f.get();
        }

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

        if (!writer) {
          writer = new rw::writer::VtkWriterInterface(compute_filename);
          if (data.d_outOnlyNodes)
            writer->appendNodes(mesh->getNodesP(), &u);
          else
            writer->appendMesh(mesh->getNodesP(), mesh->getElementType(),
                               mesh->getElementConnectivitiesP(), &u);
        }

        // append velocity
        writer->appendPointData("Mark_Velocity", &v_mark);
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

        auto f = hpx::parallel::for_loop(
            hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
            mesh->getNumNodes(), [&v_mark, data, mesh](boost::uint64_t i) {
              auto x = mesh->getNode(i);

              bool proceed = true;
              if (data.d_symmAxis == "y" &&
                  util::compare::definitelyLessThan(x.d_x,
                                                    data.d_symmLine + 1.0E-8))
                proceed = false;

              if (data.d_symmAxis == "x" &&
                  util::compare::definitelyLessThan(x.d_y,
                                                    data.d_symmLine + 1.0E-8))
                proceed = false;

              if (proceed) {
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
            }); // parallel for loop
        f.get();

        if (!writer) {
          writer = new rw::writer::VtkWriterInterface(compute_filename);
          if (data.d_outOnlyNodes)
            writer->appendNodes(mesh->getNodesP(), &u);
          else
            writer->appendMesh(mesh->getNodesP(), mesh->getElementType(),
                               mesh->getElementConnectivitiesP(), &u);
        }

        // append velocity
        writer->appendPointData("Symm_Velocity", &v_mark);
      }
      // clear the v_mark data
      if (!v_mark.empty())
        v_mark.shrink_to_fit();

      //
      // operation : Compute strain and magnitude of strain
      //
      // Steps: Call CCMFE class and compute strain and stress
      //
      if (data.d_computeStrain) {
        std::vector<util::SymMatrix3> strain(mesh->getNumElements(),
                                             util::SymMatrix3());
        std::vector<util::SymMatrix3> stress(mesh->getNumElements(),
                                             util::SymMatrix3());
        std::vector<float> magS;

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

        // compute strain and stress
        auto f = hpx::parallel::for_loop(
            hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
            mesh->getNumElements(),
            [&strain, &stress, data, mesh, quad, mat_deck,
             u](boost::uint64_t e) {
              auto ssn = util::SymMatrix3();
              auto sss = util::SymMatrix3();

              // get ids of nodes of element, coordinate of nodes, 1st order
              // quad data, and first quad data
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

              strain[e] = ssn;
              stress[e] = sss;
            }); // parallel loop over elements
        f.get();

        // compute magnitude of strain
        if (data.d_magStrainTensor) {
          magS = std::vector<float>(strain.size(), 0.);
          auto f2 = hpx::parallel::for_loop(
              hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
              mesh->getNumElements(), [&magS, strain, data](boost::uint64_t e) {
                if (data.d_magStrainComp.empty()) {
                  magS[e] = std::abs(strain[e].d_xx);
                  magS[e] += std::abs(strain[e].d_yy);
                  magS[e] += std::abs(strain[e].d_zz);
                  magS[e] += std::abs(strain[e].d_xy);
                  magS[e] += std::abs(strain[e].d_xz);
                  magS[e] += std::abs(strain[e].d_yz);
                } else if (data.d_magStrainComp == "xx") {
                  magS[e] = std::abs(strain[e].d_xx);
                } else if (data.d_magStrainComp == "yy") {
                  magS[e] = std::abs(strain[e].d_yy);
                }
              });
          f2.get();
        }

        //
        // output strain/stress data
        //
        if (!data.d_outOnlyNodes) {

          // append mesh
          if (!writer) {
            writer = new rw::writer::VtkWriterInterface(compute_filename);
            writer->appendMesh(mesh->getNodesP(), mesh->getElementType(),
                               mesh->getElementConnectivitiesP(), &u);
          }

          writer->appendCellData("Strain_Tensor", &strain);
          writer->appendCellData("Stress_Tensor", &stress);
          if (data.d_magStrainTensor)
            writer->appendCellData("Mag_Strain", &magS);

          //          // mark magnitude of strain if asked
          //          if (!data.d_markMagStrainCells.empty()) {
          //            for (auto cell : data.d_markMagStrainCells)
          //              magS[cell.first] = cell.second;
          //
          //            if (data.d_magStrainTensor)
          //              writer.appendCellData("Mark_Mag_Strain", &magS);
          //          }

        } else {

          // compute 1st order quad points and store them (quad points in
          // current configuration)
          std::vector<util::Point3> elem_quads =
              std::vector<util::Point3>(mesh->getNumElements(), util::Point3());

          auto f2 = hpx::parallel::for_loop(
              hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
              mesh->getNumElements(),
              [&elem_quads, mesh, quad, u](boost::uint64_t e) {
                auto nds = mesh->getElementConnectivity(e);
                std::vector<util::Point3> nds_current(nds.size(),
                                                      util::Point3());
                for (size_t j = 0; j < nds.size(); j++)
                  nds_current[j] = mesh->getNode(nds[j]) + u[nds[j]];

                std::vector<fe::QuadData> qds =
                    quad->getQuadPoints(nds_current);
                // store first quad point
                elem_quads[e] = qds[0].d_p;
              });
          f2.get();

          // create unstructured vtk output
          std::string fname = out_filename + "_" + data.d_tagFilename +
                              "_quads_" + std::to_string(out);
          auto writer1 = rw::writer::VtkWriterInterface(fname);
          writer1.appendNodes(&elem_quads);
          writer1.appendPointData("Strain_Tensor", &strain);
          writer1.appendPointData("Stress_Tensor", &stress);
          if (data.d_magStrainTensor)
            writer1.appendPointData("Mag_Strain", &magS);

          // mark magnitude of strain if asked
          if (!data.d_markMagStrainCells.empty()) {
            for (auto cell : data.d_markMagStrainCells)
              magS[cell.first] = cell.second;

            if (data.d_magStrainTensor)
              writer1.appendPointData("Mark_Mag_Strain", &magS);
          }
          writer1.close();
        }
      }

      //
      // operation : Compute crack tip location and velocity
      //
      std::vector<float> damage_Z;
      if (data.d_crackTip) {

        if (fractures[c] == nullptr) {
          // modify fracture deck for crack tip calculation
          if (!data.d_crackSameDtOut) {
            fracture_decks[c].d_dtCrackOut = output_deck->d_dtOut;
            fracture_decks[c].d_dtCrackVelocity = 1;
          } else {
            fracture_decks[c].d_dtCrackOut = output_deck->d_dtOut;
            fracture_decks[c].d_dtCrackVelocity = output_deck->d_dtOut;
          }
          fracture_decks[c].d_crackOutFilename = data.d_tagFilename;

          // get fracture class
          fractures[c] = new geometry::Fracture(&(fracture_decks[c]));
        }

        auto n = out * output_deck->d_dtOut;
        auto time = n * model_deck->d_dt;

        // compute displacement, damage, and crack tip at k-1 step if
        // crack update is not same as simulation output interval
        if (!data.d_crackSameDtOut) {
          n -= 1;
          time -= model_deck->d_dt;
          auto f = hpx::parallel::for_loop(
              hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
              mesh->getNumNodes(), [&u, v, model_deck](boost::uint64_t i) {
                u[i] -= util::Point3(model_deck->d_dt * v[i].d_x,
                                     model_deck->d_dt * v[i].d_y,
                                     model_deck->d_dt * v[i].d_z);
              });
          f.get();

          // compute damage at n-1
          if (damage_Z.empty())
            damage_Z = std::vector<float>(mesh->getNumNodes(), 0.);
          computeDamage(deck, model_deck, material, mesh, &u, &damage_Z);

          // compute crack tip location and crack tip velocity
          fractures[c]->updateCrackAndOutput(n, time, out_path,
                                             model_deck->d_horizon,
                                             mesh->getNodesP(), &u, &damage_Z);
        }

        // get current displacement
        n = out * output_deck->d_dtOut;
        time = n * model_deck->d_dt;
        auto f2 = hpx::parallel::for_loop(
            hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
            mesh->getNumNodes(), [&u, v, model_deck](boost::uint64_t i) {
              u[i] += util::Point3(model_deck->d_dt * v[i].d_x,
                                   model_deck->d_dt * v[i].d_y,
                                   model_deck->d_dt * v[i].d_z);
            });
        f2.get();

        // compute damage at current displacement
        if (damage_Z.empty())
          damage_Z = std::vector<float>(mesh->getNumNodes(), 0.);
        computeDamage(deck, model_deck, material, mesh, &u, &damage_Z);
        fractures[c]->updateCrackAndOutput(n, time, out_path,
                                           model_deck->d_horizon,
                                           mesh->getNodesP(), &u, &damage_Z);
      }

      //
      // operation : Compute damage at node
      //
      if (data.d_damageAtNodes) {

        // compute if it is not computed already
        if (damage_Z.empty())
          computeDamage(deck, model_deck, material, mesh, &u, &damage_Z);

        if (!writer) {
          writer = new rw::writer::VtkWriterInterface(compute_filename);
          if (data.d_outOnlyNodes)
            writer->appendNodes(mesh->getNodesP(), &u);
          else
            writer->appendMesh(mesh->getNodesP(), mesh->getElementType(),
                               mesh->getElementConnectivitiesP(), &u);
        }
        // append data to file
        writer->appendPointData("Damage_Z", &damage_Z);
      }

      // clear data
      if (!damage_Z.empty())
        damage_Z.shrink_to_fit();

      //
      // operation : Compute J integral
      //
      if (data.d_compJIntegral) {

        double energy_into_crack = 0.;
        computeJIntegral(out, data, model_deck, mat_deck, material, mesh,
            neighbor_list, &u, &v, energy_into_crack);
      }

      //
      // close file
      //
      if (writer)
        writer->close();
    } // loop over compute sets
  }   // processing output files
}