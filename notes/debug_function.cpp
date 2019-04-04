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
std::ofstream writeCsv;
writeCsv.open("meshFeTest.csv");
writeCsv<<d_nodes.size()<<"\n";
for (auto nd : d_nodes)
writeCsv<<nd.d_x<<","<<","<<nd.d_y<<","<<nd.d_z<<"\n";
d_numElems = d_enc.size()/d_eNumVertex;
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
