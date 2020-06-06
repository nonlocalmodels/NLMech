# NLMech
Welcome to NLMech repository. In this project we implement 
nonlocal fracture theory, referred to as Peridynamics, 
using both finite element and finite difference discretization.  

Logo below has been obtained by running Peridynamic simulation on the logo mesh, [see](https://nonlocalmodels.github.io/examples/fd-logo-soft-material.html).


<p style="text-align:center;"><img src="assets/logo/logo_sim.png?raw=true" alt="logo" width="400"/></p>


## Building 

The code uses [CMake](https://cmake.org/) as the build system and has following dependencies: [HPX](https://github.com/STEllAR-GROUP/hpx), [Blaze](https://bitbucket.org/blaze-lib/blaze/src/master/), [Blaze_Iterative](https://github.com/STEllAR-GROUP/BlazeIterative), [Boost](https://www.boost.org/), [VTK](https://www.vtk.org), and [YAML-CPP](https://github.com/jbeder/yaml-cpp).

We provide following support to build the code

* Bash scripts to build the code on [Ubuntu](https://github.com/nonlocalmodels/buildscripts/tree/master/Ubuntu) and [Fedora](https://github.com/nonlocalmodels/buildscripts/tree/master/Fedora),
* Docker files to build the code using the [Fedora packages](https://github.com/nonlocalmodels/buildscripts/blob/master/Docker/Fedora) or [build all dependencies](https://github.com/nonlocalmodels/buildscripts/blob/master/Docker/FedoraAll),
* [Scripts](https://github.com/nonlocalmodels/HPCBuildInfrastructure) to build the code on HPC systems.

If you just interested the code, we provide a Doker image with the latest stable version and the latest sucessful build of the master branch. For writing your own 
[YAML](https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html) files, we provide a bunch of [examples](https://nonlocalmodels.github.io/examples/).

## Documentation


## Releases

The current stable version is [*0.1.0*]() and the current developement happens in the [master branch](). For more details, we refer to the [Changelog](ChANGELOG.md) file.

## Code of conduct

We have adopted a [code of conduct](.github/CODE_OF_CONDUCT.md) for this project. Please refer to this document if you would like to know more about the expectations for members of our community, with regard to how they will behave toward each other.

## Contributing
