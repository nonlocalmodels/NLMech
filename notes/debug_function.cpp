////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//  Copyright (c) 2019 Patrick Diehl
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

static int wait(void) {

  std::cin.ignore();

  return 0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
// Debug Input function setLoadingDeck()
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void debug_Input_setLoadingDeck() {

  // read and create d_loadingDeck_p and then run following code to debug
  // the output
  for (size_t s = 0; s < d_loadingDeck_p->d_uBCData.size(); s++) {
    io::BCData bc = d_loadingDeck_p->d_uBCData[s];
    std::cout << "Displacement bc set = " << s + 1 << "\n";
    std::cout << "Region type = " << bc.d_regionType << " Region = " << bc.d_x1
              << "," << bc.d_y1 << "," << bc.d_x2 << "," << bc.d_y2 << "\n";
    std::cout << "Direction = ";
    for (auto j : bc.d_direction)
      std::cout << j << ",";
    std::cout << "\n";
    std::cout << "Time function type = " << bc.d_timeFnType << " Parameters = ";
    for (auto j : bc.d_timeFnParams)
      std::cout << j << ",";
    std::cout << "\n";
    std::cout << "Spatial function type = " << bc.d_spatialFnType
              << " Parameters = ";
    for (auto j : bc.d_spatialFnParams)
      std::cout << j << ",";
    std::cout << "\n";
  }

  for (size_t s = 0; s < d_loadingDeck_p->d_fBCData.size(); s++) {
    io::BCData bc = d_loadingDeck_p->d_fBCData[s];
    std::cout << "Force bc set = " << s + 1 << "\n";
    std::cout << "Region type = " << bc.d_regionType << " Region = " << bc.d_x1
              << "," << bc.d_y1 << "," << bc.d_x2 << "," << bc.d_y2 << "\n";
    std::cout << "Direction = ";
    for (auto j : bc.d_direction)
      std::cout << j << ",";
    std::cout << "\n";
    std::cout << "Time function type = " << bc.d_timeFnType << " Parameters = ";
    for (auto j : bc.d_timeFnParams)
      std::cout << j << ",";
    std::cout << "\n";
    std::cout << "Spatial function type = " << bc.d_spatialFnType
              << " Parameters = ";
    for (auto j : bc.d_spatialFnParams)
      std::cout << j << ",";
    std::cout << "\n";
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
// MESH 1
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
  // check getElementConnectivity() method
  for (size_t i=0; i<5; i++) {

    size_t j = (i+1) * d_numElems/6;

    // get connectivity using standard method
    std::vector<size_t> j_con_1;
    for (size_t k=0; k<d_eNumVertex; k++)
      j_con_1.emplace_back(d_enc[j*d_eNumVertex + k]);

    // get connectivity using function
    std::vector<size_t> j_con_2 = getElementConnectivity(j);

    if (j_con_1 != j_con_2)
      std::cerr << "Error in connectivity of element = "<<j<<"\n";
  }

  std::cout<<"Num nodes = "<<d_nodes.size()<<" "<<d_enc.size()<<" "
                                                                ""<<d_nec
                                                                .size()<<"\n";

 std::cout<<d_nodes[0].d_x<<" "<<d_nodes[0].d_y<<" "<<d_nodes[0].d_z<<"\n";
  std::cout<<d_nodes[9].d_x<<" "<<d_nodes[9].d_y<<" "<<d_nodes[9].d_z<<"\n";

    std::ofstream myfile("nodes.csv");
  myfile.precision(6);
  for (size_t i=0; i<10; i++)
    myfile << i <<","<< d_nodes[i].d_x << "," << d_nodes[i].d_y << "," <<
    d_nodes[i].d_z <<"\n";
  myfile.close();

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
// MESH 2
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

// check fixity
for (size_t i=0; i<d_fix.size(); i++) {

if (!(d_fix[i] & FIX_X_MASK) and i%1000 == 0)
std::cout << "x is free "<<d_fix[i]<<"\n";

d_fix[i] |= FIX_X_MASK;

if (d_fix[i] & FIX_X_MASK and i%1000 == 0)
std::cout << "x is fixed "<<d_fix[i]<<"\n";

if (!(d_fix[i] & FIX_Y_MASK) and i%1000 == 0)
std::cout << "y is free "<<d_fix[i]<<"\n";

d_fix[i] |= FIX_Y_MASK;

if (d_fix[i] & FIX_Y_MASK and i%1000 == 0)
std::cout << "y is fixed "<<d_fix[i]<<"\n";

if (!(d_fix[i] & FIX_Z_MASK) and i%1000 == 0)
std::cout << "z is free "<<d_fix[i]<<"\n";

d_fix[i] |= FIX_Z_MASK;

if (d_fix[i] & FIX_Z_MASK and i%1000 == 0)
std::cout << "z is fixed "<<d_fix[i]<<"\n";

}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
// MESH 3 Output mesh in csv file
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
std::cout<<"Num elems = "<<d_numElems<<"\n";

std::ofstream writeCsv;
writeCsv.open("meshFeTest.csv");
writeCsv<<d_nodes.size()<<"\n";
for (auto nd : d_nodes)
writeCsv<<nd.d_x<<","<<nd.d_y<<","<<nd.d_z<<"\n";
writeCsv<<d_numElems<<"\n";
for (size_t e=0; e< d_numElems; e++) {
for (size_t i=0; i<d_eNumVertex; i++) {
writeCsv << d_enc[d_eNumVertex*e + i];
if (i == d_eNumVertex - 1) {
if (e < d_numElems-1)
writeCsv << "\n";
}
else
writeCsv<<",";
}
}
writeCsv.close();

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
// MESH 4 Debug quadrature
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
  // debug print
  std::ofstream writeCsv;
  std::string writeFile = "triMesh_quads_ex2_" + std::to_string(n) + ".csv";
  writeCsv.open(writeFile.c_str());
  for (size_t e = 0; e < num_elems; e++) {
    std::vector<util::Point3> enodes = {nodes[elements[num_vertex * e + 0]],
                                        nodes[elements[num_vertex * e + 1]],
                                        nodes[elements[num_vertex * e + 2]]};

    std::cout << "*** e = " << e << " nodes = "
              <<enodes[0].d_x << "," << enodes[0].d_y << ","
              <<enodes[1].d_x << "," << enodes[1].d_y << ","
              <<enodes[2].d_x << "," << enodes[2].d_y << "\n";

    std::vector<fe::QuadData> qds = quad.getQuadPoints(enodes);
    double sum = 0.;
    for (auto qd : qds) {
      writeCsv << qd.d_w << "," << qd.d_p.d_x << "," << qd.d_p.d_y << "\n";
      sum += qd.d_w;

      std::cout << "  quads = " << qd.d_w << "," << qd.d_p.d_x << "," << qd.d_p.d_y << "\n";
    }
//    if (std::abs(sum - 0.125) > tol)
//      std::cout << "Error: Sum of quad points do not match the area of "
//                   "elements.\n";
  }
  writeCsv.close();

//++++

    //
    // Check quad points for triangle {(1,0), (2,0), (0,2)}
    //
    nodes = {util::Point3(1., 0., 0.), util::Point3(2., 0., 0.),
             util::Point3(0., 2., 0.)};
    qds = quad.getQuadPoints(nodes);
    std::cout<<"Check quads on new triangle\n";
    for (auto qd: qds)
      std::cout<< qd.d_w << "," << qd.d_p.d_x << "," << qd.d_p.d_y << "\n";

    if (n==1) {
      fe::TriElem *triE = new fe::TriElem(n);
      util::Point3 p = util::Point3();
      std::vector<double> shapes = {1. - 1./3. - 1./3., 1./3., 1./3.};
      std::vector<std::vector<double>> der_shapes;
      double j = triE->mapRefElemToElem(p, shapes, der_shapes, nodes);

      std::cout<< j*qds[0].d_w << "," << p.d_x << "," << p.d_y << "\n";
    }

// ++++++++

std::string filename = "tri_quads_" + std::to_string(n) + ".txt";
std::ofstream myfile(filename.c_str());
myfile.precision(10);
size_t counter = 0;
double sum = 0.;
for (auto qd : qds) {
myfile << ++counter << " " << qd.d_p.d_x << " " << qd.d_p.d_y << " "
<< qd.d_w << " " << qd.d_shapes[0] << " " << qd.d_shapes[1] << " "
<< qd.d_shapes[2] << " " << qd.d_derShapes[0][0] << " "
<< qd.d_derShapes[0][1] << " " << qd.d_derShapes[1][0] << " "
<< qd.d_derShapes[0][1] << " " << qd.d_derShapes[2][0] << " "
<< qd.d_derShapes[2][1] << "\n";
sum += qd.d_w;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
// Debug bitwise operations
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
// Example program
#include <bitset>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

#define FREE_MASK 0x000
#define FIX_X_MASK 0x001
#define FIX_Y_MASK 0x002
#define FIX_Z_MASK 0x004

struct fixity {
  bool x;
  bool y;
  bool z;
  fixity() : x(false), y(false), z(false){};
};

int main() {
  vector<fixity> f1(100, fixity());
  vector<bool> f2(300, false);
  vector<size_t> f3(100, 0);
  vector<short int> f4(100, 0);

  vector<char> f5(100, 0);

  cout << "Size of fixity = " << sizeof(fixity) << "\n";
  cout << "Size of bool = " << sizeof(bool) << "\n";
  cout << "Size of int = " << sizeof(int) << "\n";
  cout << "Size of short int = " << sizeof(short int) << "\n";
  cout << "Size of char = " << sizeof(char) << "\n";

  cout << "Size of f1 = " << sizeof(fixity) * f1.capacity() << "\n";
  cout << "Size of f2 = " << sizeof(bool) * f2.capacity() << "\n";
  cout << "Size of f3 = " << sizeof(size_t) * f3.capacity() << "\n";
  cout << "Size of f4 = " << sizeof(short int) * f4.capacity() << "\n";
  cout << "Size of f5 = " << sizeof(char) * f5.capacity() << "\n";

  // set fixity of all f5
  f5.resize(10);
  cout << "*****************\n";
  for (size_t i = 0; i < f5.size(); i++) {
    bitset<8> b(f5[i]);
    cout << b << "\n";
    if (i % 3 == 0)
      f5[i] |= FIX_X_MASK;
    else if (i % 3 == 1)
      f5[i] |= FIX_Y_MASK;
    else if (i % 3 == 2)
      f5[i] |= FIX_Z_MASK;
  }

  cout << "*****************\n";
  for (size_t i = 0; i < f5.size(); i++) {
    bitset<8> b(f5[i]);
    cout << b << "\n";
  }
}


++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
From FDModel.cpp
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void model::FDModel::integrate() {

  create_neighbors(d_mesh_p, neighbors);

  // at the beginning compute forces and apply initial and boundary condition

  // initial condition
  if (d_n == 0)
    d_initialCondition_p->apply(&d_u, &d_v, d_mesh_p);

  // boundary condition
  d_uLoading_p->apply(d_time, &d_u, &d_v, d_mesh_p);
  //  apply_u_bc(d_time, &d_u, &d_v, d_mesh_p);
  //  d_fLoading_p->apply(d_time, &d_f, d_mesh_p);

  // internal forces
  computeForces();

  // perform output at the beginning
  if (d_n == 0) {
    if (d_policy_p->enablePostProcessing())
      computePostProcFields();

    output();
  }

  // start time integration
  size_t i = d_n;
  for (i; i < d_modelDeck_p->d_Nt; i++) {
    //    if (d_modelDeck_p->d_timeDiscretization == "central_difference")
    integrateCD();
    //    else if (d_modelDeck_p->d_timeDiscretization == "velocity_verlet")
    //      integrateVerlet();

    if ((d_n % d_outputDeck_p->d_dtOut == 0) &&
        (d_n >= d_outputDeck_p->d_dtOut)) {
      if (d_policy_p->enablePostProcessing())
        computePostProcFields();

      output();
    }
  } // loop over time steps
}

void model::FDModel::integrateVerlet() {

  // step 1 and 2 : Compute v_mid and u_new
  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_mesh_p->getNumNodes(), [this](boost::uint64_t i) {
        size_t dim = this->d_mesh_p->getDimension();
        double delta_t = this->d_modelDeck_p->d_dt;
        double density = this->d_material_p->getDensity();

        // modify dofs which are not marked fixed
        if (this->d_mesh_p->isNodeFree(i, 0)) {

          this->d_v[i].d_x += 0.5 * delta_t * (this->d_f[i].d_x +
              this->d_fext[i].d_x)/ density;

          this->d_u[i].d_x += delta_t * this->d_v[i].d_x;
        }

        if (dim > 1)
          if (this->d_mesh_p->isNodeFree(i, 1)) {

            this->d_v[i].d_y += 0.5 * delta_t * (this->d_f[i].d_y +
                this->d_fext[i].d_y)/ density;

            this->d_u[i].d_y += delta_t * this->d_v[i].d_y;
          }

        if (dim > 2)
          if (this->d_mesh_p->isNodeFree(i, 2)) {

            this->d_v[i].d_z += 0.5 * delta_t * this->d_f[i].d_z / density;

            this->d_u[i].d_z += delta_t * this->d_v[i].d_z;
          }

//        // reset force
//        this->d_f[i].d_x = 0.;
//        this->d_f[i].d_y = 0.;
//        this->d_f[i].d_z = 0.;
      }); // end of parallel for loop

  f.get();

  // compute forces and energy due to new displacement field (this will be
  // used in next time step)
  d_n++;
  d_time += d_modelDeck_p->d_dt;

  // boundary condition
  d_uLoading_p->apply(d_time, &d_u, &d_v, d_mesh_p);
  d_fLoading_p->apply(d_time, &d_fext, d_mesh_p);

  // internal forces
  computeForces();

  // Step 3: Compute v_new
  f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_mesh_p->getNumNodes(), [this](boost::uint64_t i) {
        size_t dim = this->d_mesh_p->getDimension();
        double delta_t = this->d_modelDeck_p->d_dt;
        double density = this->d_material_p->getDensity();

        // modify dofs which are not marked fixed
        if (this->d_mesh_p->isNodeFree(i, 0))
          this->d_v[i].d_x += 0.5 * delta_t * (this->d_f[i].d_x +
              this->d_fext[i].d_x)/ density;

        if (dim > 1)
          if (this->d_mesh_p->isNodeFree(i, 1))
            this->d_v[i].d_y += 0.5 * delta_t * (this->d_f[i].d_y +
                this->d_fext[i].d_y) / density;

        if (dim > 2)
          if (this->d_mesh_p->isNodeFree(i, 2))
            this->d_v[i].d_z += 0.5 * delta_t * this->d_f[i].d_z / density;
      }); // end of parallel for loop

  f.get();
}

void model::FDModel::computeForces() {

  if (d_material_p->isStateActive())
    computeHydrostaticStrains();

  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_mesh_p->getNumNodes(),
      [this](boost::uint64_t i) {
        // local variable to hold force
        util::Point3 force_i = util::Point3();

        // reference coordinate and displacement at the node
        util::Point3 xi = this->d_mesh_p->getNode(i);
        util::Point3 ui = this->d_u[i];

        // get hydrostatic energy and force
        std::pair<double, double> gi;
        if (this->d_material_p->isStateActive())
          gi = d_material_p->getStateEF(this->d_hS[i]);

        // get interior flag
        auto node_i_interior = this->d_interiorFlags_p->getInteriorFlag(i, xi);

        // upper and lower bound for volume correction
        double h = d_mesh_p->getMeshSize();
        double check_up = d_modelDeck_p->d_horizon + 0.5 * h;
        double check_low = d_modelDeck_p->d_horizon - 0.5 * h;

        // inner loop over neighbors
        auto i_neighs = this->d_neighbor_p->getNeighbors(i);
        for (size_t j = 0; j < i_neighs.size(); j++) {

          size_t j_id = i_neighs[j];

          // there are two contributions to force at node i
          // 1. From bond j-i due to bond-based forces
          // 2. From hydrostatic strains at node j and node i due to
          // hydrostatic forces

          // compute bond-based contribution
          util::Point3 xj = this->d_mesh_p->getNode(j_id);
          util::Point3 uj = this->d_u[j_id];
          util::Point3 xji = xj - xi;
          util::Point3 uji = uj - ui;
          double rji = xji.length();
          double Sji = this->d_material_p->getS(xji, uji);

          // get corrected volume of node j
          double volj = this->d_mesh_p->getNodalVolume(j_id);
//          if (util::compare::definitelyGreaterThan(rji, check_low))
//            volj *= (check_up - rji) / h;

          // check if bond is in no-fail region
          bool break_bonds = true;
          if (!node_i_interior ||
              !this->d_interiorFlags_p->getInteriorFlag(j_id, xj))
            break_bonds = false;

          /*

          // get peridynamics force and energy density between bond i and j
          bool fs = this->d_fracture_p->getBondState(i, j);
          std::pair<double, double> ef =
              this->d_material_p->getBondEF(rji, Sji, fs, break_bonds);

          // update the fractured state of bond
          //          this->d_fracture_p->setBondState(i, j, fs);

          // compute the contribution of bond force to force at i
          double scalar_f = ef.second * volj / rji;

          force_i.d_x += scalar_f * xji.d_x;
          force_i.d_y += scalar_f * xji.d_y;
          force_i.d_z += scalar_f * xji.d_z;

          */

          auto bond_fractured = this->d_fracture_p->getBondState(i, j);
          auto ef =
              getFandFprimeNew(rji, Sji, bond_fractured, break_bonds);

          // update the fractured state of bond
//          this->d_fracture_p->setBondState(i, j, bond_fractured);

          // compute influence function
          auto wji = getInfFn(rji);

          // compute the contribution of bond force f_ij to force at i

          // compute the contribution of bond force to force at i
          double scalar_f = wji * ef.second * volj;

          force_i.d_x += scalar_f * xji.d_x;
          force_i.d_y += scalar_f * xji.d_y;
          force_i.d_z += scalar_f * xji.d_z;

          // compute state-based contribution
          if (d_material_p->isStateActive() &&
              d_material_p->addBondContribToState(Sji, rji)) {

            // Compute gj while noting that gi is already computed
            std::pair<double, double> gj =
                d_material_p->getStateEF(this->d_hS[j_id]);

            double scalar_g =
                d_material_p->getStateForce(gi.second + gj.second, rji) / rji;

            force_i.d_x += scalar_g * xji.d_x;
            force_i.d_y += scalar_g * xji.d_y;
            force_i.d_z += scalar_g * xji.d_z;
          }
        } // loop over neighboring nodes

        // update force and energy
        this->d_f[i] = force_i;
      } // loop over nodes

  ); // end of parallel for loop

  f.get();
}

void model::FDModel::computeHydrostaticStrains() {

  std::cout << "Here hs\n";

  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_mesh_p->getNumNodes(),
      [this](boost::uint64_t i) {
        // local variable to hold strain
        double hydro_strain_i = 0.;

        // reference coordinate and displacement at the node
        util::Point3 xi = this->d_mesh_p->getNode(i);
        util::Point3 ui = this->d_u[i];

        // upper and lower bound for volume correction
        double h = d_mesh_p->getMeshSize();
        double check_up = d_modelDeck_p->d_horizon + 0.5 * h;
        double check_low = d_modelDeck_p->d_horizon - 0.5 * h;

        // inner loop over neighbors
        auto i_neighs = this->d_neighbor_p->getNeighbors(i);
        for (size_t j = 0; j < i_neighs.size(); j++) {

          size_t j_id = i_neighs[j];

          // compute bond-based contribution
          util::Point3 uji = this->d_u[j_id] - ui;
          util::Point3 xji = this->d_mesh_p->getNode(j_id) - xi;
          double rji = xji.length();
          double Sji = this->d_material_p->getS(xji, uji);

          // get corrected volume of node j
          double volj = this->d_mesh_p->getNodalVolume(j_id);
          if (util::compare::definitelyGreaterThan(rji, check_low))
            volj *= (check_up - rji) / h;

          // compute the contribution of bond to hydrostatic strain at i
          if (d_material_p->addBondContribToState(Sji, rji)) {

            hydro_strain_i +=
                volj * d_material_p->getBondContribToHydroStrain(Sji, rji);
          }
        } // loop over neighboring nodes

        // update
        this->d_hS[i] = hydro_strain_i;
      } // loop over nodes

  ); // end of parallel for loop

  f.get();
}

void model::FDModel::computePostProcFields() {

  // if work done is to be computed, get the external forces
  std::vector<util::Point3> fext(d_mesh_p->getNumNodes(), util::Point3());
  if (d_policy_p->populateData("Model_d_w"))
    d_fLoading_p->apply(d_time, &fext, d_mesh_p);

  // local data for kinetic energy
  std::vector<float> vec_ke(d_mesh_p->getNumNodes(), 0.);

  auto f = hpx::parallel::for_loop(
      hpx::parallel::execution::par(hpx::parallel::execution::task), 0,
      d_mesh_p->getNumNodes(),
      [this, fext, &vec_ke](boost::uint64_t i) {
        // local variable
        double energy_i = 0.0;
        double hydro_energy_i = 0.0;
        double a = 0.; // for damage
        double b = 0.; // for damage
        double z = 0.; // for damage

        // reference coordinate and displacement at the node
        util::Point3 xi = this->d_mesh_p->getNode(i);
        util::Point3 ui = this->d_u[i];

        // get volume of node i
        double voli = this->d_mesh_p->getNodalVolume(i);

        // get hydrostatic energy and force
        std::pair<double, double> gi;
        if (this->d_material_p->isStateActive())
          gi = d_material_p->getStateEF(this->d_hS[i]);

        // get interior flag
        auto node_i_interior = this->d_interiorFlags_p->getInteriorFlag(i, xi);

        // upper and lower bound for volume correction
        double h = d_mesh_p->getMeshSize();
        double check_up = d_modelDeck_p->d_horizon + 0.5 * h;
        double check_low = d_modelDeck_p->d_horizon - 0.5 * h;

        // inner loop over neighbors
        auto i_neighs = this->d_neighbor_p->getNeighbors(i);
        for (size_t j = 0; j < i_neighs.size(); j++) {

          size_t j_id = i_neighs[j];

          // there are two contributions to force at node i
          // 1. From bond j-i due to bond-based forces
          // 2. From hydrostatic strains at node j and node i due to
          // hydrostatic forces

          // compute bond-based contribution
          util::Point3 xj = this->d_mesh_p->getNode(j_id);
          util::Point3 uj = this->d_u[j_id];
          util::Point3 xji = xj - xi;
          double rji = xji.length();
          double Sji = this->d_material_p->getS(xji, uj - ui);

          // get corrected volume of node j
          double volj = this->d_mesh_p->getNodalVolume(j_id);
          if (util::compare::definitelyGreaterThan(rji, check_low))
            volj *= (check_up - rji) / h;

          // check if bond is in no-fail region
          bool break_bonds = true;
          if (!node_i_interior ||
              !this->d_interiorFlags_p->getInteriorFlag(j_id, xj))
            break_bonds = false;

          // get peridynamics force and energy density between bond i and j
          bool fs = this->d_fracture_p->getBondState(i, j);
          std::pair<double, double> ef =
              this->d_material_p->getBondEF(rji, Sji, fs, break_bonds);

          // energy
          energy_i += ef.first * volj;

          // parameters for damage function \phi
          if (!fs)
            a += volj;
          b += volj;

          // parameters for damage function Z
          double sr = 0.;
          if (util::compare::definitelyGreaterThan(rji, 1.0E-12))
            sr = std::abs(Sji) / this->d_material_p->getSc(rji);
          if (util::compare::definitelyLessThan(z, sr))
            z = sr;
        } // loop over neighboring nodes

        // compute hydrostatic energy
        if (this->d_material_p->isStateActive())
          hydro_energy_i =
              gi.first / (d_modelDeck_p->d_horizon * d_modelDeck_p->d_horizon);

        if (this->d_policy_p->populateData("Model_d_e"))
          this->d_e[i] = (energy_i + hydro_energy_i) * voli;

        if (this->d_policy_p->populateData("Model_d_w"))
          this->d_w[i] = ui.dot(fext[i]);

        if (this->d_policy_p->populateData("Model_d_eFB") &&
            util::compare::definitelyGreaterThan(z, 1.0 - 1.0E-10))
          this->d_eFB[i] = energy_i * voli;

        if (this->d_policy_p->populateData("Model_d_eF") &&
            util::compare::definitelyGreaterThan(z, 1.0 - 1.0E-10))
          this->d_eF[i] = (energy_i + hydro_energy_i) * voli;

        if (this->d_policy_p->populateData("Model_d_phi"))
          this->d_phi[i] = 1. - a / b;

        if (this->d_policy_p->populateData("Model_d_Z"))
          this->d_Z[i] = z;

        // compute kinetic energy
        vec_ke[i] = 0.5 * this->d_material_p->getDensity() *
            this->d_v[i].dot(this->d_v[i]);
      } // loop over nodes

  ); // end of parallel for loop

  f.get();

  // add energies to get total energy
  if (this->d_policy_p->populateData("Model_d_e"))
    d_te = util::methods::add(d_e);
  if (this->d_policy_p->populateData("Model_d_w"))
    d_tw = util::methods::add(d_w);
  if (this->d_policy_p->populateData("Model_d_eF"))
    d_teF = util::methods::add(d_eF);
  if (this->d_policy_p->populateData("Model_d_eFB"))
    d_teFB = util::methods::add(d_eFB);

  d_tk = util::methods::add(vec_ke);
}
