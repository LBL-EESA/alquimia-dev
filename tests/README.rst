=============
Test problems
=============

Automated test driver
---------------------

Since the batch chem drivers are not doing anything except calling the
engine reaction stepping many times, the numerical results should be
identical (to the precision of the text file output) at every time
step.

The automated test driver compares the pflotran observation problem
with a similarly structured alquimia batch_chem output file.

To run a all tests in a particular config file:
::

    python ./batch-compare-alquimia-pflotran.py -c alquimia-pflotran-tests.cfg -p ./pflotran -a ../src/drivers/batch_chem


To run a single test by name:
::

    python ./batch-compare-alquimia-pflotran.py -c alquimia-pflotran-tests.cfg -p ./pflotran -a ../src/drivers/batch_chem -t calcite-volume-fractions-pflotran-constraint





Calcite kinetics, short time, and pflotran native constraints
-------------------------------------------------------------

Test the alquimia interface with mineral dissolution/precipitation for
a short simulation time (no changes to calcite volume fraction). The
initial condition is specified by name only, causing alquimia look for
a pflotran native constraint with that name.

**Status: fail** 

Status Notes: Looking at the time series, there are occasional differences in the sixth decimal.

::

  ../install/bin/batch_chem -d -i calcite-short-pc.cfg

  ./pflotran -input_prefix calcite-short


NOTE: this is the same pflotran input file as
calcite-kinetics-short-dc!

Calcite kinetics, short time, and driver generated constraints
--------------------------------------------------------------

Test the alquimia interface with mineral dissolution/precipitation for
a short simulation time (no changes to calcite volume fraction). The
initial condition is fully specified by the driver and processed by
pflotran.

**Status: fail**

Status Notes: Looking at the time series, there are occasional differences in the sixth decimal.

::

  ../install/bin/batch_chem -d -i calcite-short-dc.cfg

  ./pflotran -input_prefix calcite-short


NOTE: this is the same pflotran input file as
calcite-kinetics-short-pc!


Calcite kinetics w/ volume fraction updates and pflotran native constraints
---------------------------------------------------------------------------

Test the alquimia interface with mineral dissolution/precipitation for
a long simulation time, so that the mineral volume fractions are
updated during reaction stepping. Uses pflotran native constraints.

**Status: fails**

Status Notes: slight initial numerical differences in rates accumulate error?

::

  ../install/bin/batch_chem -d -i calcite-vf-pc.cfg

  ./pflotran -input_prefix calcite-vf


