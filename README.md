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


To build alquimia using the waf based buildsystem run something like the following:

    cd ${ALQUIMIA_DIR}/
    export PFLOTRAN_DIR=/path/to/pflotran-dev
    CC=/opt/local/bin/openmpicc CXX=/opt/local/bin/openmpicxx  FC=/opt/local/bin/openmpif90 ./waf --prefix ./install --pflotran=${PFLOTRAN_DIR} distclean configure build install
    cd ${ALQUIMIA_DIR}/examples
    ../install/bin/batch_chem -i batch-test-1.cfg -d
