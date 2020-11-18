# CMake options

## NLMech 

* Enable_Documentation : Generates target for generating the documentation (Default = False)
* Enable_Tools : Enables the tools to the build target (Default = False)
* Enable_Expensive_Tests : Enables the computation intense tests (Default = False)
* Enable_RPM: Enables to generate RPM packages (Default = False)

## General options

* CMAKE_BUILD_TYPE : Specifies the build type (Default = Release)
* CMAKE_INSTALL_PREFIX : Specifies the install directory (Default = /usr/local)

## Dependencies

* VTK_DIR : Path the the VTK installation
* BLAZEITERATIVE_INCLUDE : Path to the Blaze_Iterative installation 
* GMSH_DIR : Path to the Gmsh installation 
* YAML_CPP_DIR : Path to the YAML CPP installation 
* HPX_DIR : Path to the HPX installation 

Note that HPX's CMake will load the same Boost version as HPX was compiled with. Note that it is important 
to use the same Boost version for this code. 

