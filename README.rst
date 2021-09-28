Travis CI Status 
-----

Master:

.. image:: https://travis-ci.org/LBL-EESA/alquimia-dev.svg?branch=master
    :target: https://travis-ci.org/LBL-EESA/alquimia-dev

Stable:

.. image:: https://travis-ci.org/LBL-EESA/alquimia-dev.svg?branch=stable
    :target: https://travis-ci.org/LBL-EESA/alquimia-dev
	     
Description
-----------

Alquimia provides a generic interface between flow and transport simulators (drivers) and biogeochemical reaction capabilities (engines). Alquimia consists of two parts: an application programming interface, API, and a wrapper library. The API describes in detail how information is exchanged between the driver and the engine. This includes the function call signatures and data structures required for implementing alquimia in a mixed language (C/C++/Fortran) programming environment. The library is an implementation of the documented API, providing an interface into the biogeochemistry routines supplied by other codes: Alquimia does not do any geochemical calculations. Version 1.0 provides interfaces to the open source codes (BSD) CrunchFlow/CrunchTope and (LGPL) PFLOTRAN. Version 1.0 meets the policies set forth by (and is part of) the Extreme-scale Scientific Software Development Kit, xSDK version 0.6.0.

Originally developed as part of the DOE ASCEM project, it is now mantained and developed under the 
U.S. Department of Energy's `IDEAS Watersheds <https://ideas-productivity.org/>`_ project

Currently, the geochemical engines suported are `CrunchFlow <https://bitbucket.org/crunchflow/crunchtope-dev>`_ and
`PFLOTRAN <https://bitbucket.org/pflotran/pflotran-dev>`_.


Legal
-----

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

Citing Alquimia
---------------

Andre, B., Molins, S., Johnson, J., and Steefel, C.I. Alquimia. Computer Software. https://github.com/LBL-EESA/alquimia-dev. USDOE. 01 Aug. 2013. Web. `doi:10.11578/dc.20210416.49 <https://doi.org/10.11578/dc.20210416.49>`_.


Building
--------

You'll need working C and Fortran compilers and CMake installed on your system.
For UNIX and UNIX-like operating systems, you'll need GNU Make or another 
capable version of Make installed as well. To build on Windows, you'll need 
some recent version of Visual Studio and a decent Fortran compiler such as 
Intel's.

Both engines require PETSc to be installed, with the PETSC_DIR and 
PETSC_ARCH environment variables set properly. Currently, PETSc must be 
configured to use 32-bit indices.

PFlotran engine
===============

Currently, Alquimia only works with a particular version of PFlotran: 
hash 611092f80ddb from the pflotran-dev repository. You can download this 
revision directly as a ZIP file from 
https://bitbucket.org/pflotran/pflotran-dev/get/611092f80ddb.zip

*NOTE ABOUT BUILDING WITH PETSC 3.6 or later: This version of PFlotran was 
written to use PETSC 3.5.x, which is slightly different from the later minor 
releases of PETSc. If you use a later version of PETSc, please note the following:*

*1. You must create the following symbolic links within $PETSC_DIR:*

::

  ln -s $PETSC_DIR/lib/petsc/conf $PETSC_DIR/conf
  ln -s $PETSC_DIR/include/petsc/finclude $PETSC_DIR/include/finclude

*2. You will see a linking error (for a missing symbol _petsclogbegin_) when 
building the pflotran_rxn executable. This doesn't prevent libpflotran_rxn.a 
from being built, nor does it prevent Alquimia from working properly with PFlotran.*

The instructions below assume that you are on a UNIX or UNIX-like system, 
and you have set the environment variable PFLOTRAN_DIR to the top of your 
PFlotran source directory.

::

    cd $PFLOTRAN_DIR/src/pflotran
    make pflotran_rxn

To build PFlotran on Windows, see the instructions 
`here <https://bitbucket.org/pflotran/pflotran-dev/wiki/Installation/Windows_with_Visual_Studio>`_.

CrunchFlow engine
=================

The CrunchFlow geochemistry engine is located in a special "alquimia" branch
of the crunchtope repository on bitbucket. Currently, you need to be a 
collaborator to access this repository, but steps are being taken to release 
an open-source version.

When you have the alquimia branch of the repository located at $CRUNCHFLOW_DIR, 
you can build the geochemistry reaction library by typing

::

    cd $CRUNCHFLOW_DIR
    make libcrunchchem.a

At this time, building CrunchFlow's geochemistry engine on Windows is not 
supported.

Alquimia interface
==================

When you have built all the desired chemistry engines, you can build the 
Alquimia interface. On UNIX and UNIX-like systems, you can use the following 
command, which assumes you have set ALQUIMIA_DIR to the top of your Alquimia 
source tree. Note that you will need to create a build tree from which to 
invoke CMake.

:: 

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

**NOTE**: you can omit either of the engines if you aren't building them both. 
If you don't specify any chemistry engine, Alquimia will halt and remind you 
that building it without an engine is pointless. So, for example, to build 
Alquimia with an install of PFlotran at $PFLOTRAN_DIR using MPI compilers, 
in Debug mode, to be installed in /usr/local:

:: 

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

If you are using a geochemical engine that requires PETSc, and you want to 
specify the exact locations of its headers, and the method for linking against 
PETSc's libraries, you can specify these with the -DTPL_PETSC_INCLUDE_DIRS=<list of dirs> and 
-DTPL_PETSC_LDFLAGS=<link flags> arguments. Normally, these options are 
omitted and Alquimia automatically detects PETSc's location using the PETSC_DIR
and PETSC_ARCH environment variables.

Testing
-------

To run Alquimia's suite of tests from your build directory, just type

::

    make test

See the CMakeLists.txt file for other available build options, including
optimization level, shared/static libraries, build prefix, etc. Alquimia 
supports all xSDK-compliant build options, which can be passed to CMake 
when configuring your build.

Installation
------------

You can install the Alquimia library and the demo drivers into your desired 
location, type

::

    make install

This will install libraries into ${CMAKE_INSTALL_PREFIX}/lib, headers into 
${CMAKE_INSTALL_PREFIX}/include/alquimia, and the demo drivers into 
${CMAKE_INSTALL_PREFIX}/bin. To run some basic sanity checks on these installed
drivers, you can type

::

    make test_install

This will run a few benchmark tests to make sure that the executables have been 
properly linked and installed.
