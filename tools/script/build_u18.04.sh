#!/bin/bash

# local directory with external libraries
local="/home/prashant/Softwares/local"

# clean previous build files (optional)
./clean.sh

# hpx
hpx_ver="1.1.0"
hpx_dir="$local""/hpx/hpx-""$hpx_ver"

# vtk
# for ubuntu 18.04, can install vtk dev using
# sudo apt-get install libvtk7-dev

# yaml-cpp
yamlcpp_ver="0.5.3"
yamlcpp_dir="$local""/yaml-cpp/yaml-cpp-""$yamlcpp_ver"

# blaze
blaze_ver="3.5"
blaze_dir="$local""/blaze/blaze-""$blaze_ver"

# target directory where code will be built
target_build=$pwd

# source
source="../NLMech"

# run cmake with flags
cmake -DHPX_DIR="$hpx_dir""/lib/cmake/HPX" \
 	  -DYAML_CPP_DIR="$yamlcpp_dir" \
 	  -Dblaze_DIR="$blaze_dir""/share/blaze/cmake" \
 	  -Dblaze_INCLUDE_DIR="$blaze_dir""/include" \
	  -DCMAKE_INSTALL_PREFIX="$target_build" \
	  -DEnable_Documentation=ON \
	  -DCMAKE_BUILD_TYPE=Release \
	  "$source"

# make
make -j 8

# create documentation
make doc 
