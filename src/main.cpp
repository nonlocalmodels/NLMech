// Copyright (c) 2019    Prashant K. Jha
//
// Distributed under the GNU GENERAL PUBLIC LICENSE, Version 3.0.
// (See accompanying file LICENSE.txt)

#include <hpx/hpx_main.hpp> // Need main source file

#include "inp/input.h"        // Input class
#include "model/model.h" // Model class
#include "model/fd/fDModel.h" // Model class
#include "fe/mesh.h"          // Mesh class
#include <iostream>

#include <algorithm>
#include <boost/program_options.hpp> // program options

//! Main driver
int main(int argc, char *argv[]) {

  boost::program_options::options_description desc("Allowed options");
  desc.add_options()("help", "produce help message")(
      "input-file,i", boost::program_options::value<std::string>(),
      "Configuration file");

  boost::program_options::variables_map vm;
  boost::program_options::store(
      boost::program_options::parse_command_line(argc, argv, desc), vm);
  boost::program_options::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 1;
  }

  // read input file
  std::string filename = "";
  if (vm.count("input-file"))
    filename = vm["input-file"].as<std::string>();

  if (filename == "") {
    std::cerr << argv[0] << " -i input.yaml --hpx:threads=n \n";
    exit(1);
  }

  // record current time
  //  std::uint64_t t = hpx::util::high_resolution_clock::now();

  // read input data
  inp::Input *deck = new inp::Input(filename);

  // check which model to run
  //  if (deck->getSpatialDiscretization() == "finite_difference")
  model::FDModel fdModel(deck);

//  //
//  // compute nodal volume
//  //
//  std::vector<double> a(1000,0);
//  auto f = hpx::parallel::for_loop(
//      hpx::parallel::execution::par(hpx::parallel::execution::task),
//      0, a.size(), [&a](boost::uint64_t i)
//      {
//
//        a[i] = double(i);
//      }
//  ); //end of parallel for loop
//
//  f.get();
//
//  for (size_t i=0; i<10; i++)
//    std::cout<<a[i]<<"\n";

  // get time elapsed
  //  std::uint64_t elapsed = hpx::util::high_resolution_clock::now() - t;
  // get number of threads
  //  std::uint64_t const os_thread_count = hpx::get_os_thread_count();

  //  std::cout << "Number of threads= " << os_thread_count
  //            << " Time elapsed= " << elapsed / 1e9 << "\n";

  return EXIT_SUCCESS;
}