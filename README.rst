Alquimia is an biogeochemistry API and wrapper library being developed
as part of the [ASCEM](http://esd.lbl.gov/research/projects/ascem/)
project.

The aim is to provide a unified interface to existing "geochemistry
engines" such as
[CrunchFlow](http://www.csteefel.com/CrunchFlowIntroduction.html) or
[PFloTran](https://bitbucket.org/pflotran/pflotran-dev), allowing
subsurface flow and transport simulators to access a range of
functionality.

It is not an implementation of a biogeochemistry engine, and does not
do any geochemical calculations.

Makefile based buildsystem
--------------------------

::

    cd ${ALQUIMIA_DIR}/src
    CC=openmpicc CXX=openmpicxx FC=openmpif90 PFLOTRAN_DIR=${PFLOTRAN_DIR} \
        make all

    cd ${ALQUIMIA_DIR}/examples
    ../src/drivers/batch_chem -d -i calcite-kinetics-vf-pc.cfg



Waf based buildsystem
---------------------

To build alquimia using the waf based buildsystem run something like the following:

::

    cd ${ALQUIMIA_DIR}/
    export PFLOTRAN_DIR=/path/to/pflotran-dev
    CC=openmpicc CXX=openmpicxx  FC=openmpif90 \
        ./waf --prefix ./install --pflotran=${PFLOTRAN_DIR} configure build install
    cd ${ALQUIMIA_DIR}/examples
    ../install/bin/batch_chem -d -i calcite-kinetics-vf-pc.cfg

To clobber/clean up the build and start fresh:

::
    CC=openmpicc CXX=openmpicxx  FC=openmpif90 ./waf --prefix ./install --pflotran=${PFLOTRAN_DIR} uninstall distclean

NOTE: uninstall will fail if you modify the wscript files to add new
install files. Run waf uninstall first!
