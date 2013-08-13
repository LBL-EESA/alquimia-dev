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
[PFloTran](https://bitbucket.org/pflotran/pflotran-dev), allowing
subsurface flow and transport simulators to access a range of
functionality.

It is not an implementation of a biogeochemistry reaction library, and
does not do any geochemical calculations.


Building
--------

To build alquimia, you must have petsc installed, with the PETSC_DIR
and PETSC_ARCH environment variables set.

::

    cd ${ALQUIMIA_DIR}/src
    CC=openmpicc CXX=openmpicxx FC=openmpif90 PFLOTRAN_DIR=${PFLOTRAN_DIR} \
        make all

    cd ${ALQUIMIA_DIR}/examples
    ../src/drivers/batch_chem -d -i calcite-kinetics-vf-pc.cfg



