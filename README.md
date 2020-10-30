# NLMech

[![CircleCI](https://circleci.com/gh/nonlocalmodels/NLMech.svg?style=shield)](https://circleci.com/gh/nonlocalmodels/nonlocalheatequation) [![Codacy Badge](https://app.codacy.com/project/badge/Grade/118379d7d745464584b73e9e06f60462)](https://www.codacy.com/gh/nonlocalmodels/NLMech?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=nonlocalmodels/NLMech&amp;utm_campaign=Badge_Grade) [![Coverage Status](https://coveralls.io/repos/github/nonlocalmodels/NLMech/badge.svg?branch=master)](https://coveralls.io/github/nonlocalmodels/NLMech?branch=master)

Welcome to NLMech repository. In this project we implement 
nonlocal fracture theory, referred to as Peridynamics, 
using both finite element and finite difference discretization.  

Logo below has been obtained by running Peridynamic simulation on the logo mesh, [see](https://nonlocalmodels.github.io/examples/fd-logo-soft-material.html).

<p style="text-align:center;"><img src="https://github.com/nonlocalmodels/NLMech/blob/master/assets/logo/logo_sim.png?raw=true" alt="logo" width="400"/></p>

## Building 

The code uses [CMake](https://cmake.org/) as the build system and has following dependencies: [HPX](https://github.com/STEllAR-GROUP/hpx), [Blaze](https://bitbucket.org/blaze-lib/blaze/src/master/), [Blaze_Iterative](https://github.com/STEllAR-GROUP/BlazeIterative), [Boost](https://www.boost.org/), [VTK](https://www.vtk.org), and [YAML-CPP](https://github.com/jbeder/yaml-cpp).

We provide following support to build the code

* Bash scripts to build the code on [Ubuntu](https://github.com/nonlocalmodels/buildscripts/tree/master/Ubuntu) and [Fedora](https://github.com/nonlocalmodels/buildscripts/tree/master/Fedora),
* Docker files to build the code using the [Fedora packages](https://github.com/nonlocalmodels/buildscripts/blob/master/Docker/Fedora) or [build all dependencies](https://github.com/nonlocalmodels/buildscripts/blob/master/Docker/FedoraAll),
* [Scripts](https://github.com/nonlocalmodels/HPCBuildInfrastructure) to build the code on HPC systems.

If you just interested the code, we provide a Doker image with the [latest sucessful build](https://github.com/nonlocalmodels/NLMech/packages/384758) of the master branch. For writing your own 
[YAML](https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html) files, we provide a bunch of [examples](https://nonlocalmodels.github.io/examples/).

If you encounter any bug or want to submit a feature request, please generate an [issue](https://github.com/nonlocalmodels/NLMech/issues) on GitHub. 
Right now, we have [![GitHub issues](https://img.shields.io/github/issues/nonlocalmodels/nlmech.svg)](https://github.com/nonlocalmodels/NLMech/issues). For contributing, we refer to the Contributing section below.

## Documentation

The documentation of the the master branch is available [here](https://nonlocalmodels.github.io/documentation/).

## Releases

The current stable version is [![GitHub release](https://img.shields.io/github/release/nonlocalmodels/NLMech.svg)](https://GitHub.com/nonlocalmodels/NLMech/releases/) and the current developement happens in the [master branch](https://github.com/nonlocalmodels/NLMech). For more details, we refer to the [Changelog]() file.

## Code of conduct

We have adopted a [code of conduct](https://github.com/nonlocalmodels/NLMech/blob/master/CODE_OF_CONDUCT.md) for this project. Please refer to this document if you would like to know more about the expectations for members of our community, with regard to how they will behave toward each other.

## Contributing

The source code is released under the [![GitHub license](https://img.shields.io/github/license/nonlocalmodels/nonlocalmodels.github.io.svg)](https://github.com/nonlocalmodels/nonlocalmodels.github.io/blob/master/LICENSE) license. If you like to contribute, we only accept your pull request using the same license. Please feel free to add your name to license header of the files you added or contributed to. If possible please add a test for your new feature using [CTest](https://gitlab.kitware.com/cmake/community/-/wikis/doc/ctest/Testing-With-CTest). We adapted the Google C++ [Style Guide](https://google.github.io/styleguide/cppguide.html) for this project. We use the [clang-format](https://clang.llvm.org/docs/ClangFormat.html) tool to format the soruce code with respect to this style guide. Please run the `format.sh' script before your do your pull request.
