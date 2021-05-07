# NLMech 

[![CircleCI](https://circleci.com/gh/nonlocalmodels/NLMech.svg?style=shield)](https://circleci.com/gh/nonlocalmodels/NLMech) [![Codacy Badge](https://app.codacy.com/project/badge/Grade/118379d7d745464584b73e9e06f60462)](https://www.codacy.com/gh/nonlocalmodels/NLMech?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=nonlocalmodels/NLMech&amp;utm_campaign=Badge_Grade) [![Coverage Status](https://coveralls.io/repos/github/nonlocalmodels/NLMech/badge.svg?branch=main)](https://coveralls.io/github/nonlocalmodels/NLMech?branch=main) [![status](https://joss.theoj.org/papers/271dd66ea91b7fbfdccb4b10a7ba462c/status.svg)](https://joss.theoj.org/papers/271dd66ea91b7fbfdccb4b10a7ba462c)

Welcome to NLMech repository. In this project we implement 
nonlocal fracture theory, referred to as Peridynamics, 
using both finite element and finite difference discretization.  

Logo below has been obtained by running Peridynamic simulation on the logo mesh, [see](https://nonlocalmodels.github.io/examples/fd-logo-soft-material.html).

<p style="text-align:center;"><img src="https://github.com/nonlocalmodels/NLMech/blob/main/assets/logo/logo_sim.png?raw=true" alt="logo" width="400"/></p>

## Building 

The code uses 

  * [CMake](https://cmake.org/) 3.16

as the build system and has following dependencies: 

  * [HPX](https://github.com/STEllAR-GROUP/hpx) 1.6.0
  * [Blaze](https://bitbucket.org/blaze-lib/blaze/src/master/) 3.8, 
  * [Blaze_Iterative](https://github.com/STEllAR-GROUP/BlazeIterative) master, 
  * [Boost](https://www.boost.org/) 1.75,  
  * [gmsh](https://gmsh.info/) 4.7 ,
  * [VTK](https://www.vtk.org) 9.0, and 
  * [YAML-CPP](https://github.com/jbeder/yaml-cpp) 0.6.

Following dependencies are optional, but are recommended for large simulations:

  * [PCL](https://github.com/PointCloudLibrary/pcl) 1.11.

Since we use the latest stable [Fedora](https://getfedora.org/) workstation spin for our continuous testing, we recommend
to use a Docker image to build the code. Thus, these are the minimal required versions are meet. For the master branch the same versions
should work. We strongly recommend to use the same CMake version to build HPX and NLMech.

We provide following support to build the code

* Bash scripts to build the code on [Ubuntu](https://github.com/nonlocalmodels/buildscripts/tree/main/Ubuntu) and [Fedora](https://github.com/nonlocalmodels/buildscripts/tree/main/Fedora),
* Docker files to build the code using the [Fedora packages](https://github.com/nonlocalmodels/buildscripts/blob/main/Docker/Fedora) or [build all dependencies](https://github.com/nonlocalmodels/buildscripts/blob/main/Docker/FedoraAll),
* [Scripts](https://github.com/nonlocalmodels/HPCBuildInfrastructure) to build the code on HPC systems.

If you just interested in running the application, we provide a Doker image with the [latest Successful build](https://github.com/nonlocalmodels/NLMech/packages/384758) of the main branch. For writing your own 
[YAML](https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html) files, we provide a bunch of [examples](https://nonlocalmodels.github.io/examples/).

If you encounter any bug or want to submit a feature request, please generate an [issue](https://github.com/nonlocalmodels/NLMech/issues) on GitHub. 
Right now, we have [![GitHub issues](https://img.shields.io/github/issues/nonlocalmodels/nlmech.svg)](https://github.com/nonlocalmodels/NLMech/issues). For contributing, we refer to the Contributing section below.

## Documentation

The documentation of the the main branch is available [here](https://nonlocalmodels.github.io/documentation/).

## Releases

The current stable version is [![GitHub release](https://img.shields.io/github/release/nonlocalmodels/NLMech.svg)](https://GitHub.com/nonlocalmodels/NLMech/releases/) and the current development happens in the [main branch](https://github.com/nonlocalmodels/NLMech). For more details, we refer to the [Changelog](https://github.com/nonlocalmodels/NLMech/blob/main/CHANGELOG.md) file.

## Code of conduct

We have adopted a [code of conduct](https://github.com/nonlocalmodels/NLMech/blob/main/CODE_OF_CONDUCT.md) for this project. Please refer to this document if you would like to know more about the expectations for members of our community, with regard to how they will behave toward each other.

## Contributing

The source code is released under the [![GitHub license](https://img.shields.io/github/license/nonlocalmodels/nonlocalmodels.github.io.svg)](https://github.com/nonlocalmodels/nonlocalmodels.github.io/blob/main/LICENSE) license. If you like to contribute, we only accept your pull request using the same license. Please feel free to add your name to license header of the files you added or contributed to. If possible please add a test for your new feature using [CTest](https://gitlab.kitware.com/cmake/community/-/wikis/doc/ctest/Testing-With-CTest). We adapted the Google C++ [Style Guide](https://google.github.io/styleguide/cppguide.html) for this project. We use the [clang-format](https://clang.llvm.org/docs/ClangFormat.html) tool to format the source code with respect to this style guide. Please run the `format.sh' script before your do your pull request.

## Citing

In publications, please use our paper as the main citation for NLMech: 

* P. Diehl, P. K. Jha, H. Kaiser, R. Lipton, and M. Lévesque. An asynchronous and task-based implementation of peridynamics utilizing hpx—the C++ standard library for parallelism and concurrency. SN Applied Sciences, 2(12):2144, 2020, [10.1007/s42452-020-03784-x](https://doi.org/10.1007/s42452-020-03784-x), [Preprint](https://arxiv.org/abs/1806.06917).

For more references, we refer to NLMech's [publication list](https://nonlocalmodels.github.io/publications/).

## Acknowledgments

NLMech has been funded by:

* Army Research Office Grant # W911NF-16-1-0456 to PI Dr. Robert Lipton (Professor at Louisiana State University). This grant supported Prashant K. Jha on a postdoctoral position from October 2016 - July 2019.
*  Canada Research Chairs Program under the Canada Research Chair in Multiscale Modelling of Advanced Aerospace Materials held by M. Lévesque; Natural Sciences and Engineering Research Council of Canada (NSERC) Discovery Grants Program under Discovery Grant RGPIN-2016-06412.
* We are grateful for the support of the [Google Summer of Code program](https://summerofcode.withgoogle.com/) funding internships.



