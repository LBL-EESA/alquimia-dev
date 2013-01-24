=============
Test problems
=============


Calcite kinetics, short time, and pflotran native constraints
-------------------------------------------------------------

Test the alquimia interface with mineral dissolution/precipitation for
a short simulation time. The initial condition is specified by name
only, causing alquimia look for a pflotran native constraint with that
name.

**Status: fails** total component concentrations are correct, but pH values are different.

.. sourcecode:: bash
  ../install/bin/batch_chem -d -i calcite-kinetics-short-pc.cfg

  ./pflotran -pflotranin calcite-kinetics-short.in -output_prefix calcite-kinetics-short-pc

NOTE: this is the same pflotran input file as
calcite-kinetics-short-dc!

Calcite kinetics, short time, and driver generated constraints
--------------------------------------------------------------

Test the alquimia interface with mineral dissolution/precipitation for
a short simulation time. The initial condition is fully specified by
the driver and processed by pflotran.

**Status: fails** (driver driven constraints not implemented)

.. sourcecode:: bash
  ../install/bin/batch_chem -d -i calcite-kinetics-short-dc.cfg

  ./pflotran -pflotranin calcite-kinetics-short.in -output_prefix calcite-kinetics-short-dc

NOTE: this is the same pflotran input file as
calcite-kinetics-short-pc!


Calcite kinetics w/ volume fraction updates and pflotran native constraints
---------------------------------------------------------------------------

Test the alquimia interface with mineral dissolution/precipitation for
a long simulation time, so that the mineral volume fractions are
updated during reaction stepping. Uses pflotran native constraints.

**Status: fails** (volume fractions not updating)

.. sourcecode:: bash
  ../install/bin/batch_chem -d -i calcite-kinetics-vf-pc.cfg

  ./pflotran -pflotranin calcite-kinetics-vf.in -output_prefix calcite-kinetics-vf


