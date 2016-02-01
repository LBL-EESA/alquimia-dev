Legal
-----

"Alquimia Copyright (c) 2013-2015, The Regents of the University of
California, through Lawrence Berkeley National Laboratory (subject to
receipt of any required approvals from the U.S. Dept. of Energy).  All
rights reserved.

If you have questions about your rights to use or distribute this
software, please contact Berkeley Lab's Technology Transfer and
Intellectual Property Management at TTD@lbl.gov referring to "Alquimia
(LBNL Ref. 2013-119)."

 

NOTICE.  This software was developed under funding from the
U.S. Department of Energy.  As such, the U.S. Government has been
granted for itself and others acting on its behalf a paid-up,
nonexclusive, irrevocable, worldwide license in the Software to
reproduce, prepare derivative works, and perform publicly and display
publicly.  Beginning five (5) years after the date permission to
assert copyright is obtained from the U.S. Department of Energy, and
subject to any subsequent five (5) year renewals, the U.S. Government
is granted for itself and others acting on its behalf a paid-up,
nonexclusive, irrevocable, worldwide license in the Software to
reproduce, prepare derivative works, distribute copies to the public,
perform publicly and display publicly, and to permit others to do so.


Description
-----------

Alquimia is a biogeochemistry API and wrapper library, originally developed as 
part of the DOE ASCEM project, and is an interoperable component of the 
Department of Energy's [IDEAS](https://ideas-productivity.org/) software 
productivity project.

The aim is to provide a unified interface to existing "geochemistry
engines" such as
[CrunchFlow](http://www.csteefel.com/CrunchFlowIntroduction.html) or
[PFLOTRAN](https://bitbucket.org/pflotran/pflotran-dev), allowing
subsurface flow and transport simulators to access a range of
functionality.

It is not an implementation of a biogeochemistry reaction library, and
does not do any geochemical calculations.

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
[here](https://bitbucket.org/pflotran/pflotran-dev/wiki/Installation/Windows_with_Visual_Studio).

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
in Debug mode:

:: 

    cd $ALQUIMIA_DIR
    mkdir build ; cd build
    cmake .. \
      -DCMAKE_C_COMPILER=`which mpicc` \
      -DCMAKE_CXX_COMPILER=`which mpicxx` \
      -DCMAKE_Fortran_COMPILER=`which mpif90` \
      -DCMAKE_BUILD_TYPE=Debug \
      -DXSDK_WITH_PFLOTRAN=ON \
      -DTPL_PFLOTRAN_LIBRARIES=$PFLOTRAN_DIR/src/pflotran/libpflotranchem.a \
      -DTPL_PFLOTRAN_INCLUDE_DIRS=$PFLOTRAN_DIR/src/pflotran
    make 

*Windows instructions go here.*

Testing
-------

To run Alquimia's suite of tests from your build directory, just type

::

    make test

See the CMakeLists.txt file for other available build options, including
optimization level, shared/static libraries, build prefix, etc. Alquimia 
supports all xSDK-compliant build options, which can be passed to CMake 
when configuring your build.

