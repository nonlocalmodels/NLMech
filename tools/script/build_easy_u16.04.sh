#
# export paths
#
local="/home/prashant/Softwares/local_pd"

# hpx
hpx_ver="1.1.0"
hpx_dir="$local""/hpx/hpx-""$hpx_ver"

# yaml-cpp
yamlcpp_ver="0.5.3"
yamlcpp_dir="$local""/yaml-cpp/yaml-cpp-""$yamlcpp_ver"
export PATH=$PATH:"$yamlcpp_dir"
export PATH=$PATH:"$yamlcpp_dir""/include"
export PATH=$PATH:"$yamlcpp_dir""/lib"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"$yamlcpp_dir""/lib"

# vtk
vtk_ver="8.0.1"
vtk_dir="$local""/vtk/vtk-""$vtk_ver"
export PATH=$PATH:"$vtk_dir"
export PATH=$PATH:"$vtk_dir""/include"
export PATH=$PATH:"$vtk_dir""/lib"
export PATH=$PATH:"$vtk_dir""/bin"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"$vtk_dir""/lib"

#
#	target where we want to build PeridynamicsHPX
#
PWD=$pwd
target_build="$PWD"

# source
source="../."

#
# run cmake with flags
#
cmake -DHPX_DIR="$hpx_dir""/lib/cmake/HPX" \
		-DVTK_DIR="$vtk_dir" \
 	   -DYAML_CPP_DIR="yamlcpp_dir" \
	   -DCMAKE_INSTALL_PREFIX="$target_build" \
	   -DCMAKE_BUILD_TYPE=Release \
	   "$source"

make -j 8