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

To build alquimia, you must have petsc installed, with the PETSC_DIR
and PETSC_ARCH environment variables set.

::

    cd ${ALQUIMIA_DIR}/src
    CC=openmpicc CXX=openmpicxx FC=openmpif90 PFLOTRAN_DIR=${PFLOTRAN_DIR} \
        make all

    cd ${ALQUIMIA_DIR}/examples
    ../src/drivers/batch_chem -d -i calcite-kinetics-vf-pc.cfg



