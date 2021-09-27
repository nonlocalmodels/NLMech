# Installation 

This tutorial guides you how to install [NLMech](https://github.com/nonlocalmodels/NLMech) from scratch. In addition, we provide some [Docker files](https://github.com/nonlocalmodels/buildscripts/blob/main/Docker/) to build the code on [Fedora](https://getfedora.org/), since we use the same OS on [Circle-CI](https://app.circleci.com/pipelines/github/nonlocalmodels/NLMech). In addition, we provide some [bash script](https://github.com/nonlocalmodels/buildscripts/tree/main/bash) to install the code in one step. 

In this installation guide, we will use a set of [scripts](https://github.com/nonlocalmodels/HPCBuildInfrastructure) to build the code on HPC clusters.

## Prerequisites 

### Tools

To compile NLMech and its dependencies the following tools are needed:
* GCC compiler collection (gcc) > 4.9, however, gcc >= 8 is recommended
* [autoconf](https://www.gnu.org/software/autoconf/)
* [wget](https://www.gnu.org/software/wget/)
* [cmake](https://cmake.org/)
* [git](https://git-scm.com/)

### Dependencies

* [Lapack](http://www.netlib.org/lapack/)
* [BLAS](http://www.netlib.org/blas/)
* [Openssl](https://www.openssl.org/)
* [OpenGL](https://www.opengl.org/)

We recommend to install these libraries using the package manager. Bewlow you can find some
examples to install them using apt and dnf

```bash
apt-get install build-essential git wget cmake libssl-dev libblas-dev liblapack-dev autoconf freeglut3-dev
```

```bash
dnf install @development-tools cmake git wget blas-devel lapack-devel freeglut-devel
```

## Using the HPCBuildInfrastructure

First, we clone the [HPCBuildInfrastructure](https://github.com/nonlocalmodels/HPCBuildInfrastructure) 

```bash
git clone https://github.com/nonlocalmodels/HPCBuildInfrastructure.git
cd HPCBuildInfrastructure/
```
The uses version of each library is defined in the `config.sh`, if you need to change them. 

### Build cmake

To build all dependencies, the script `build-all.sh` is used and the last argument is always the dependency to be build.

```bash
./build-all.sh Release without-gcc cmake
```
The first command `Release` specifies the `CMAKE_BUILD_TYPE` of the build. Note that you should use the same build type for all of the dependencies. The second command specifies if you use the gcc from the system in `/usr/bin` by using `without-gcc`. If you do not have access to gcc > 4.9 on your system, we provide the option `with-gcc` to use your own build of gcc. Note you have to run `./build-all.sh Release with-gcc gcc` first to build your gcc.

### Build HPX

```bash
./build-all.sh Release without-gcc hwloc
./build-all.sh Release without-gcc jemalloc
./build-all.sh Release without-gcc boost
./build-all.sh Release without-gcc hpx
```

### Build Blaze and blaze_iterative

```bash
./build-all.sh Release without-gcc blaze
./build-all.sh Release without-gcc blazeIterative
```

### Build VTK

```bash
./build-all.sh Release without-gcc vtk
```

### Build yampl-cpp

```bash
./build-all.sh Release without-gcc yamlcpp
```

### Build NLMech

```bash
./build-all.sh Release without-gcc NLMech
```

All these instructions are available in a single [bash script](https://github.com/nonlocalmodels/buildscripts/blob/main/bash/buildAll.sh) and a [Dockerfile](https://github.com/nonlocalmodels/buildscripts/blob/main/Docker/FedoraAll).
