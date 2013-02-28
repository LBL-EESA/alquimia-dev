Batch Test problems
===================

Since the batch chem drivers are not doing anything except calling the
engine's reaction stepping many times and writing the results to a
file, the numerical results should be identical (to the precision of
the text file output) at every time step.

Automated test driver
---------------------

The automated test driver compares the pflotran observation data for
the problem with a similarly structured alquimia batch_chem output
file. The number of columns in each file must be the same. Variable
names should be the same, but the order doesn't matter. So for
example, if "Total H+ [M]", "Total Sorbed A [mol/m^3]" or "Calcite
Rate [mol/sec]" is on one file, it must be in the other.

To run all tests in a particular config file:

.. code-block:: bash

    python ./batch-compare-alquimia-pflotran.py -c alquimia-pflotran-tests.cfg -p ./pflotran -a ../src/drivers/batch_chem


To run a single test by name:

.. code-block:: bash

    python ./batch-compare-alquimia-pflotran.py -c alquimia-pflotran-tests.cfg -p ./pflotran -a ../src/drivers/batch_chem -t calcite-volume-fractions-pflotran-constraint



General Aqueous Reactions
-------------------------

The general aqueous reaction problems are the simplest possible tests
of the alquimia interface. They track only changes in the aqeous
concentrations due to aqueous reactions. There is no mineral or sorbed
phase interaction.

General reaction, pflotran native constraint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Test the alquimia interface with pflotran's general aqueous reaction, using pflotran supplied constraint.

**Status: pass** 

Status Notes: 

.. code-block:: bash

  ../install/bin/batch_chem -d -i general-reaction-pc.cfg

  ./pflotran -input_prefix calcite-short


NOTE: this is the same pflotran input file as
general-reaction-ac!

General reaction, alquimia supplied constraint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Test the alquimia interface with pflotran's general aqueous reaction, using alquimia supplied constraint.

**Status: pass** 

Status Notes: 

.. code-block:: bash

  ../install/bin/batch_chem -d -i general-reaction-ac.cfg

  ./pflotran -input_prefix calcite-short


NOTE: this is the same pflotran input file as
general-reaction-pc!



Mineral Dissolution and Precipitation
-------------------------------------


Calcite, short time steps, and pflotran native constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Test the alquimia interface with mineral dissolution/precipitation for
a short simulation time and short time steps (no changes to calcite
volume fraction). The initial condition is specified by name only,
causing alquimia look for a pflotran native constraint with that name.

**Status: fail** 

Status Notes: Looking at the time series, there are occasional differences in the sixth decimal.

.. code-block:: bash

  ../install/bin/batch_chem -d -i calcite-short-pc.cfg

  ./pflotran -input_prefix calcite-short


NOTE: this is the same pflotran input file as
calcite-short-ac!

Calcite, short time steps, and alquimia generated constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Test the alquimia interface with mineral dissolution/precipitation for
a short simulation time (no changes to calcite volume fraction). The
initial condition is fully specified by the driver and processed by
pflotran.

**Status: fail**

Status Notes: Looking at the time series, there are occasional differences in the sixth decimal.

.. code-block:: bash

  ../install/bin/batch_chem -d -i calcite-short-dc.cfg

  ./pflotran -input_prefix calcite-short


NOTE: this is the same pflotran input file as
calcite-short-pc!


Volume fraction updates and pflotran native constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Test the alquimia interface with mineral dissolution/precipitation for
a long simulation time, so that the mineral volume fractions are
updated during reaction stepping. Uses pflotran native constraints.

**Status: fails**

Status Notes: slight initial numerical differences in rates accumulate error?

.. code-block:: bash

  ../install/bin/batch_chem -d -i calcite-vf-pc.cfg

  ./pflotran -input_prefix calcite-vf


Equilibrium Sorption Isotherms
------------------------------

Test the alquimia interface with equilibrium sorption isotherms, Kd,
langmuir, and freundlich. Because these are equilibrium isotherms, the
initial equilibrium solution obtained when processing the geochemical
constraint should not change during reaction stepping.


PFloTran native constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Status: fails**

Status Notes: final significant figure is off for the "D" species.

.. code-block:: bash

    ../src/drivers/batch_chem -d -i isotherms-pc.cfg
    ./pflotran -input_prefix isotherms

NOTE: this is the same pflotran input file as isotherms-ac.

Alquimia supplied constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Status: fails**

Status Notes: final significant figure is off for the "D" species.

.. code-block:: bash

    ../src/drivers/batch_chem -d -i isotherms-ac.cfg
    ./pflotran -input_prefix isotherms

NOTE: this is the same pflotran input file as isotherms-pc.

Equilibrium Ion exchange
------------------------

Test the alquimia interface with equilibrium ion exchange. Because
these are equilibrium reactions, the initial equilibrium solution
obtained when processing the geochemical constraint should not change
during reaction stepping.


PFloTran supplied constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Status: pass**

Status Notes: 

.. code-block:: bash

    ../src/drivers/batch_chem -d -i ion-exchange-valocchi-pc.cfg
    ./pflotran -input_prefix ion-exchange-valocchi

NOTE: this is the same pflotran input file as ion-exchange-valocchi-pc.

Alquimia supplied constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Status: pass**

Status Notes: 

.. code-block:: bash

    ../src/drivers/batch_chem -d -i ion-exchange-valocchi-ac.cfg
    ./pflotran -input_prefix ion-exchange-valocchi

NOTE: this is the same pflotran input file as ion-exchange-valocchi-pc.


Equilibrium Surface Complexation
--------------------------------

Test the alquimia interface with equilibrium surface complexation for
two surface sites. Because these are equilibrium reactions, the
initial equilibrium solution obtained when processing the geochemical
constraint should not change during reaction stepping.


PFloTran supplied constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Status: fail**

Status Notes: numerical differences in trailing significant figures

.. code-block:: bash

    ../src/drivers/batch_chem -d -i ion-exchange-valocchi-pc.cfg
    ./pflotran -input_prefix ion-exchange-valocchi

NOTE: this is the same pflotran input file as ion-exchange-valocchi-pc.

Alquimia supplied constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Status: fail**

Status Notes: numerical differences in trailing significant figures

.. code-block:: bash

    ../src/drivers/batch_chem -d -i ion-exchange-valocchi-ac.cfg
    ./pflotran -input_prefix ion-exchange-valocchi

NOTE: this is the same pflotran input file as ion-exchange-valocchi-pc.

