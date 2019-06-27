#
# export paths
#
local="/home/prashant/Softwares/local"

./clean.sh

# hpx
hpx_ver="1.1.0"
hpx_dir="$local""/hpx/hpx-""$hpx_ver"

# yaml-cpp
yamlcpp_ver="0.5.3"
yamlcpp_dir="$local""/yaml-cpp/yaml-cpp-""$yamlcpp_ver"

# blaze
blaze_dir="$local""/blaze"

#
# target where we want to build PeridynamicsHPX
#
PWD=$pwd
target_build="$PWD"

# source
source="../NLMech-v0.1"

#
# run cmake with flags
#
cmake -DHPX_DIR="$hpx_dir""/lib/cmake/HPX" \
 	  -DYAML_CPP_DIR="$yamlcpp_dir" \
 	  -DBLAZE_DIR="$blaze_dir" \
	  -DCMAKE_INSTALL_PREFIX="$target_build" \
	  -DCMAKE_BUILD_TYPE=Release \
	  "$source"

make -j 8
