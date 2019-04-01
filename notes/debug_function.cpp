//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
// Debug Input function setLoadingDeck()
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void debug_Input_setLoadingDeck() {

	// read and create d_loadingDeck_p and then run following code to debug
	// the output
	for (size_t s=0; s<d_loadingDeck_p->d_uBCData.size(); s++) {
    io::BCData bc = d_loadingDeck_p->d_uBCData[s];
    std::cout<<"Displacement bc set = "<<s+1<<"\n";
    std::cout << "Region type = " << bc.d_regionType << " Region = " << bc.d_x1
              << "," << bc.d_y1 << "," << bc.d_x2 << "," << bc.d_y2 << "\n";
    std::cout<<"Direction = ";
    for (auto j : bc.d_direction)
      std::cout<<j<<",";
    std::cout<<"\n";
    std::cout<<"Time function type = "<<bc.d_timeFnType<<" Parameters = ";
    for (auto j : bc.d_timeFnParams)
      std::cout<<j<<",";
    std::cout<<"\n";
    std::cout<<"Spatial function type = "<<bc.d_spatialFnType<<" Parameters = ";
    for (auto j : bc.d_spatialFnParams)
      std::cout<<j<<",";
    std::cout<<"\n";
  }

  for (size_t s=0; s<d_loadingDeck_p->d_fBCData.size(); s++) {
    io::BCData bc = d_loadingDeck_p->d_fBCData[s];
    std::cout<<"Force bc set = "<<s+1<<"\n";
    std::cout << "Region type = " << bc.d_regionType << " Region = " << bc.d_x1
              << "," << bc.d_y1 << "," << bc.d_x2 << "," << bc.d_y2 << "\n";
    std::cout<<"Direction = ";
    for (auto j : bc.d_direction)
      std::cout<<j<<",";
    std::cout<<"\n";
    std::cout<<"Time function type = "<<bc.d_timeFnType<<" Parameters = ";
    for (auto j : bc.d_timeFnParams)
      std::cout<<j<<",";
    std::cout<<"\n";
    std::cout<<"Spatial function type = "<<bc.d_spatialFnType<<" Parameters = ";
    for (auto j : bc.d_spatialFnParams)
      std::cout<<j<<",";
    std::cout<<"\n";
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
// Debug bitwise operations
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
// Example program
#include <iostream>
#include <string>
#include <vector>
#include <bitset>

using namespace std;

#define FREE_MASK 0x000
#define FIX_X_MASK 0x001
#define FIX_Y_MASK 0x002
#define FIX_Z_MASK 0x004

struct fixity {
    bool x;
    bool y;
    bool z;
    fixity() :x(false), y(false), z(false){};
};

int main()
{
    vector<fixity> f1(100,fixity());
    vector<bool> f2(300, false);
    vector<size_t> f3(100,0);
    vector<short int> f4(100,0);
    
    vector<char> f5(100,0);
    
    cout<<"Size of fixity = "<<sizeof(fixity)<<"\n";
    cout<<"Size of bool = "<<sizeof(bool)<<"\n";
    cout<<"Size of int = "<<sizeof(int)<<"\n";
    cout<<"Size of short int = "<<sizeof(short int)<<"\n";
    cout<<"Size of char = "<<sizeof(char)<<"\n";
    
    cout<<"Size of f1 = "<<sizeof(fixity) * f1.capacity()<<"\n";
    cout<<"Size of f2 = "<<sizeof(bool) * f2.capacity()<<"\n";
    cout<<"Size of f3 = "<<sizeof(size_t) * f3.capacity()<<"\n";
    cout<<"Size of f4 = "<<sizeof(short int) * f4.capacity()<<"\n";
    cout<<"Size of f5 = "<<sizeof(char) * f5.capacity()<<"\n";
    
    // set fixity of all f5
    f5.resize(10);
    cout<<"*****************\n";
    for (size_t i=0; i<f5.size(); i++) {
        bitset<8> b(f5[i]);
        cout<<b<<"\n";
        if (i%3 == 0)
            f5[i] |= FIX_X_MASK;
        else if (i%3 == 1)
            f5[i] |= FIX_Y_MASK;
        else if (i%3 == 2)
            f5[i] |= FIX_Z_MASK;
    }
    
    cout<<"*****************\n";
    for (size_t i=0; i<f5.size(); i++) {
        bitset<8> b(f5[i]);
        cout<<b<<"\n";
    }
    
}

