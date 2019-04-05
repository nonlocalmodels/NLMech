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
