Legal
-----

"Alquimia Copyright (c) 2013, The Regents of the University of
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

Alquimia is an biogeochemistry API and wrapper library being developed
as part of the [ASCEM](http://esd.lbl.gov/research/projects/ascem/)
project.

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
capable version of Make installed as well.

To build Alquimia with support for PFlotran's chemistry engine, you must have 
PETSc 3.5.x installed, with the PETSC_DIR and PETSC_ARCH environment variables 
set. 

Currently, Alquimia only works with a particular version of PFlotran: 
hash 611092f80ddb from the pflotran-dev repository at 
https://bitbucket.org/pflotran/pflotran-dev/wiki/Home

The instructions below assume that you have set the environment variable 
PFLOTRAN_DIR to the top of your PFlotran source directory.

::

    cd ${PFLOTRAN_DIR}/src/pflotran
    make pflotran_rxn

Finally, you can build Alquimia using the following instructions, which assume
you have set ALQUIMIA_DIR to the top of your Alquimia source tree. Note that 
you will need to create a build tree from which to invoke CMake.

::

    cd ${ALQUIMIA_DIR}
    mkdir build ; cd build
    cmake .. \
      -DCMAKE_C_COMPILER=<C compiler> \
      -DCMAKE_Fortran_COMPILER=<Fortran compiler> \
      -DXSDK_WITH_PFLOTRAN=ON \
      -DTPL_PFLOTRAN_LIBRARIES=$PFLOTRAN_DIR/src/pflotran/libpflotranchem.a \
      -DTPL_PFLOTRAN_INCLUDE_DIRS=$PFLOTRAN_DIR/src/pflotran
    make 

NOTE: The \*PFLOTRAN\* arguments can be omitted for builds without the PFlotran
chemistry engine. You also don't need PETSc in this case.

To run Alquimia's suite of tests from your build directory, just type

::

    make test

See the CMakeLists.txt file for other available build options, including
optimization level, shared/static libraries, build prefix, etc. Alquimia 
supports all xSDK-compliant build options, which can be passed to CMake 
when configuring your build.

