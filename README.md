## Status

[![Alquimia CI](https://github.com/LBL-EESA/alquimia-dev/actions/workflows/dev.yml/badge.svg)](https://github.com/LBL-EESA/alquimia-dev/actions/workflows/dev.yml)

## Description

Alquimia provides a generic interface between flow and transport simulators (drivers) and biogeochemical reaction capabilities (engines). Alquimia consists of two parts: an application programming interface, API, and a wrapper library. The API describes in detail how information is exchanged between the driver and the engine. This includes the function call signatures and data structures required for implementing alquimia in a mixed language (C/C++/Fortran) programming environment. The library is an implementation of the documented API, providing an interface into the biogeochemistry routines supplied by other codes: Alquimia does not do any geochemical calculations. Version 1.x provides interfaces to the open source codes (BSD) CrunchFlow/CrunchTope and (LGPL) PFLOTRAN. Version 1.x meets the policies set forth by (and is part of) the Extreme-scale Scientific Software Development Kit, xSDK version 0.6.0.

Originally developed as part of the DOE ASCEM project, it is now mantained and developed under the U.S. Department of Energy's [IDEAS Watersheds](https://ideas-watersheds.github.io/) project.

Currently, the geochemical engines suported are [CrunchFlow](https://bitbucket.org/crunchflow/crunchtope-dev) and [PFLOTRAN](https://bitbucket.org/pflotran/pflotran-dev).

## Legal

Alquimia Copyright (c) 2013-2021, The Regents of the University of
California, through Lawrence Berkeley National Laboratory (subject
to receipt of any required approvals from the U.S. Dept. of Energy). 
All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative 
works, and perform publicly and display publicly, and to permit others to do so.

## Citing Alquimia

The article:

[![DOI](https://img.shields.io/badge/doi-10.5194%2Fgmd--18--3241--2025-blue)](https://doi.org/10.5194/gmd-18-3241-2025)

Molins, S., Andre, B. J., Johnson, J. N., Hammond, G. E., Sulman, B. N., Lipnikov, K., et al. (2025). Alquimia v1.0: a generic interface to biogeochemical codes – a tool for interoperable development, prototyping and benchmarking for multiphysics simulators. Geoscientific Model Development, 18(11), 3241–3263. [doi:10.5194/gmd-18-3241-2025](https://doi.org/10.5194/gmd-18-3241-2025)

The software:

[![DOI](https://img.shields.io/badge/doi-10.11578%2Fdc%2E20210416.49-1abc9c.svg)](https://doi.org/10.11578/dc.20210416.49)

Andre, B., Molins, S., Johnson, J., and Steefel, C.I. Alquimia. Computer Software.
https://github.com/LBL-EESA/alquimia-dev. USDOE. 01 Aug. 2013. Web.
[doi:10.11578/dc.20210416.49](https://doi.org/10.11578/dc.20210416.49).

## Building

You'll need working C and Fortran compilers and CMake installed on your system.
For UNIX and UNIX-like operating systems, you'll need GNU Make or another 
capable version of Make installed as well. To build on Windows, you'll need 
some recent version of Visual Studio and a decent Fortran compiler such as 
Intel's.

### Required packages and versions

Because Alquimia is simply an interface, to be built, it requires PETSc and at least one of the two geochemical engines, either PFLOTRAN or CrunchFlow.

Alquimia is part of the [Extreme-scale Scientific Software Development Kit (xSDK)](https://xsdk.info), along with PETSc and PFLOTRAN. xSDK releases ensure that certain version of these software packages will build together. In addition to the instructions that follow, note that Alquimia -like the other xSDK packages- may also be built using the package manager [Spack](https://spack.io). 

To build Alquimia, use the version of the packages in
[the latest release of the xSDK](https://xsdk.info/releases/)
to ensure compatibility. Currently:

|            | Version    |
|------------|------------|
|xSDK        | 1.0.0      |
|Alquimia    | 1.1.0      |
|PETSc       | 3.20.0     |
|PFLOTRAN    | 5.0.0      |
|CrunchFlow  | dev        |

CrunchFlow is currently not part of the xSDK but generally the development
branch in [CrunchFlow](https://bitbucket.org/crunchflow/crunchtope-dev)
will work.

### Automated Build (Superbuild)

If you do not have PETSc, PFLOTRAN, or CrunchFlow installed, you can use the
automated build system to download and build them for you.

```bash
mkdir build && cd build
cmake -DALQUIMIA_SUPERBUILD=ON ..
make
```

This will install everything locally in the build/install directory.
You may use one of the -DXSDK_WITH_PFLOTRAN=OFF or -DXSDK_WITH_CRUNCHFLOW=OFF options to skip building one of the engines but not both at the same time. Default is ON for both.

##### Example Configuration:
```bash
mkdir build && cd build
cmake .. \
  -DALQUIMIA_SUPERBUILD=ON \
  -DXSDK_WITH_PFLOTRAN=ON \
  -DXSDK_WITH_CRUNCHFLOW=ON \
make -j$(nproc)
```

See what [Alquimia CMake configure options](#Alquimia-CMake-configure-options) below for additional CMake options.

### Manual Build

You may install the required packages manually and provide them to the CMake system to build Alquimia itself.

### PETSc

[PETSc](https://petsc.org) is a suite of data structures and routines for
the scalable (parallel) solution of scientific applications modeled by partial
differential equations. The PETSc requirement in Alquimia stems from the fact
that both engines, PFLOTRAN or CrunchFlow, require PETSc.

To download and install PETSc, please follow the instructions in
[petsc.org](https://petsc.org) or in
[pflotran.org](http://doc-dev.pflotran.org/user_guide/how_to/installation/installation.html). 
At the end of the installation, the PETSC_DIR and PETSC_ARCH environment variables
must set properly.

### PFLOTRAN engine

[PFLOTRAN](https://www.pflotran.org) is an open source, state-of-the-art
massively parallel subsurface flow and reactive transport code. Alquimia provides
access to the geochemical capabilities of PFLOTRAN; more speficically, the
geochemical capabilities available under the operator splitting mode.

Follow the instruction to download and build PFLOTRAN found
[here](http://doc-dev.pflotran.org/user_guide/how_to/installation/installation.html),
Do not build the pflotran target rather pflotran_rxn:

```bash
cd $PFLOTRAN_DIR/src/pflotran
make pflotran_rxn
```

To build PFLOTRAN on Windows, see the instructions 
[here](https://bitbucket.org/pflotran/pflotran-dev/wiki/Installation/Windows_with_Visual_Studio).

### CrunchFlow engine

[CrunchFlow](https://bitbucket.org/crunchflow/crunchtope-dev)
is a powerful software package for simulating reactive transport
developed by Carl Steefel and co-workers and applied since 1988 to a variety
of problems in the earth and environmental sciences. Alquimia provides access
to the geochemical capabilities of CrunchFlow; more speficically, the
geochemical capabilities available under the operator splitting mode.

Download the master branch of CrunchFlow, apply the makefile patch and build
the libcrunchchem.a target: 

```bash
cd $CRUNCHFLOW_DIR/source
git apply MakefileForAlquimia.patch
make libcrunchchem.a
```
### Alquimia CMake configure options

When you have built all the desired chemistry engines manually, you can build the 
Alquimia interface. On UNIX and UNIX-like systems, you can use the following 
command, which assumes you have set ALQUIMIA_DIR to the top of your Alquimia 
source tree. Note that you will need to create a build tree from which to 
invoke CMake.

```bash
cd $ALQUIMIA_DIR
mkdir build ; cd build
cmake .. \
  -DCMAKE_INSTALL_PREFIX=<installation prefix> \
  -DCMAKE_C_COMPILER=<C compiler> \
  -DCMAKE_CXX_COMPILER=<C++ compiler> \
  -DCMAKE_Fortran_COMPILER=<Fortran compiler> \
  -DXSDK_WITH_PFLOTRAN=ON \
  -DTPL_PFLOTRAN_LIBRARIES=$PFLOTRAN_DIR/src/pflotran/libpflotranchem.a \
  -DTPL_PFLOTRAN_INCLUDE_DIRS=$PFLOTRAN_DIR/src/pflotran \
  -DXSDK_WITH_CRUNCHFLOW=ON \
  -DTPL_CRUNCHFLOW_LIBRARIES=$CRUNCHFLOW_DIR/libcrunchchem.a \
  -DTPL_CRUNCHFLOW_INCLUDE_DIRS=$CRUNCHFLOW_DIR
make 
```

**NOTE**: you can omit either of the engines if you aren't building them both. 
If you don't specify any chemistry engine, Alquimia will halt and remind you 
that building it without an engine is pointless. So, for example, to build 
Alquimia with an install of PFlotran at $PFLOTRAN_DIR using MPI compilers, 
in Debug mode, to be installed in /usr/local:

```bash
cd $ALQUIMIA_DIR
mkdir build ; cd build
cmake .. \
  -DCMAKE_INSTALL_PREFIX=/usr/local \
  -DCMAKE_C_COMPILER=`which mpicc` \
  -DCMAKE_CXX_COMPILER=`which mpicxx` \
  -DCMAKE_Fortran_COMPILER=`which mpif90` \
  -DCMAKE_BUILD_TYPE=Debug \
  -DXSDK_WITH_PFLOTRAN=ON \
  -DTPL_PFLOTRAN_LIBRARIES=$PFLOTRAN_DIR/src/pflotran/libpflotranchem.a \
  -DTPL_PFLOTRAN_INCLUDE_DIRS=$PFLOTRAN_DIR/src/pflotran
make 
```

If you are using a geochemical engine that requires PETSc, and you want to 
specify the exact locations of its headers, and the method for linking against 
PETSc's libraries, you can specify these with the -DTPL_PETSC_INCLUDE_DIRS=<list of dirs> and 
-DTPL_PETSC_LDFLAGS=<link flags> arguments. Normally, these options are 
omitted and Alquimia automatically detects PETSc's location using the PETSC_DIR
and PETSC_ARCH environment variables.

## Testing

To run Alquimia's suite of tests from your build directory, just type

```bash
make test
```

See the CMakeLists.txt file for other available build options, including
optimization level, shared/static libraries, build prefix, etc. Alquimia 
supports all xSDK-compliant build options, which can be passed to CMake 
when configuring your build.

## Installation

You can install the Alquimia library and the demo drivers into your desired 
location, type

```bash
make install
```

This will install libraries into `${CMAKE_INSTALL_PREFIX}/lib`, headers into 
`${CMAKE_INSTALL_PREFIX}/include/alquimia`, and the demo drivers into 
`${CMAKE_INSTALL_PREFIX}/bin`. To run some basic sanity checks on these installed
drivers, you can type

```bash
make test_install
```

This will run a few benchmark tests to make sure that the executables have been 
properly linked and installed.

## Building Standalone Chemistry Engines with Alquimia

The Alquimia build system also supports building the standalone executables for its backend chemistry engines (CrunchFlow and PFLOTRAN) alongside the core Alquimia interface libraries as needed. 

### How to Enable Standalone Builds

When configuring Alquimia with CMake via the Superbuild, you can enable the standalone builds by passing the `-DALQUIMIA_BUILD_STANDALONE_ENGINES=ON` flag. 

### Example Configuration:
```bash
mkdir build && cd build
cmake .. \
  -DALQUIMIA_SUPERBUILD=ON \
  -DXSDK_WITH_PFLOTRAN=ON \
  -DXSDK_WITH_CRUNCHFLOW=ON \
  -DALQUIMIA_BUILD_STANDALONE_ENGINES=ON
make -j$(nproc)
```

### Where to Find the Executables

For this example configuration, once the compilation finishes, the standalone executables will be copied to the central install directory located in your CMake binary folder:

```text
build/
└── install/
    ├── bin/
    │   ├── crunchflow
    │   └── pflotran
    ├── lib/
    │   ├── libalquimia.a
    │   ├── libcrunchchem.a
    │   └── libpflotranchem.a
    └── include/
```

You can execute them directly from `build/install/bin/` or add this directory to your system's `PATH`.


