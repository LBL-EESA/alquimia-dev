..
   Alquimia Copyright (c) 2013-2016, The Regents of the University of California, 
   through Lawrence Berkeley National Laboratory (subject to receipt of any 
   required approvals from the U.S. Dept. of Energy).  All rights reserved.
   
   Alquimia is available under a BSD license. See LICENSE.txt for more
   information.
   
   If you have questions about your rights to use or distribute this software, 
   please contact Berkeley Lab's Technology Transfer and Intellectual Property 
   Management at TTD@lbl.gov referring to Alquimia (LBNL Ref. 2013-119).
   
   NOTICE.  This software was developed under funding from the U.S. Department 
   of Energy.  As such, the U.S. Government has been granted for itself and 
   others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide 
   license in the Software to reproduce, prepare derivative works, and perform 
   publicly and display publicly.  Beginning five (5) years after the date 
   permission to assert copyright is obtained from the U.S. Department of Energy, 
   and subject to any subsequent five (5) year renewals, the U.S. Government is 
   granted for itself and others acting on its behalf a paid-up, nonexclusive, 
   irrevocable, worldwide license in the Software to reproduce, prepare derivative
   works, distribute copies to the public, perform publicly and display publicly, 
   and to permit others to do so.
   
   Authors: Benjamin Andre <bandre@lbl.gov>
    

Batch Test problems
===================

Since the batch chem drivers are not doing anything except calling the
engine's reaction stepping many times and writing the results to a
file, the numerical results should be identical (to the precision of
the text file output) at every time step.

Automated test driver
---------------------

The automated test driver compares the PFloTran observation data for
the problem with a similarly structured alquimia batch_chem output
file. The number of columns in each file must be the same. Variable
names should be the same, but the order doesn't matter. So for
example, if "Total H+ [M]", "Total Sorbed A [mol/m^3]" or "Calcite
Rate [mol/sec]" is on one file, it must be in the other.

To run all tests in a particular config file:

.. code-block:: bash

    python ./batch-compare-alquimia-pflotran.py -c alquimia-pflotran-tests.cfg -p ./pflotran -a ./batch_chem


To run a single test by name:

.. code-block:: bash

    python ./batch-compare-alquimia-pflotran.py -c alquimia-pflotran-tests.cfg -p ./pflotran -a ./batch_chem -t calcite-volume-fractions-pflotran-constraint



General Aqueous Reactions
-------------------------

The general aqueous reaction problems are the simplest possible tests
of the alquimia interface. They track only changes in the aqueous
concentrations due to aqueous reactions. There is no mineral or sorbed
phase interaction.

General reaction, PFloTran native constraint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Test the alquimia interface with PFloTran's general aqueous reaction, using PFloTran supplied constraint.

**Status: pass** 

Status Notes: 

.. code-block:: bash

  ./batch_chem -d -i general-reaction-pc.cfg

  ./pflotran -input_prefix general-reaction


NOTE: this is the same PFloTran input file as
general-reaction-ac!

General reaction, alquimia supplied constraint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Test the alquimia interface with PFloTran's general aqueous reaction, using alquimia supplied constraint.

**Status: pass** 

Status Notes: 

.. code-block:: bash

  ./batch_chem -d -i general-reaction-ac.cfg

  ./pflotran -input_prefix general-reaction


NOTE: this is the same PFloTran input file as
general-reaction-pc!



Mineral Dissolution and Precipitation
-------------------------------------


Calcite, short time steps, and PFloTran native constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Test the alquimia interface with mineral dissolution/precipitation for
a short simulation time and short time steps (no changes to calcite
volume fraction). The initial condition is specified by name only,
causing alquimia look for a PFloTran native constraint with that name.

**Status: fail** 

Status Notes: Looking at the time series, there are occasional differences in the sixth decimal.

.. code-block:: bash

  ./batch_chem -d -i calcite-short-pc.cfg

  ./pflotran -input_prefix calcite-short


NOTE: this is the same PFloTran input file as
calcite-short-ac!

Calcite, short time steps, and alquimia generated constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Test the alquimia interface with mineral dissolution/precipitation for
a short simulation time (no changes to calcite volume fraction). The
initial condition is fully specified by the driver and processed by
PFloTran.

**Status: fail**

Status Notes: Looking at the time series, there are occasional differences in the sixth decimal.

.. code-block:: bash

  ./batch_chem -d -i calcite-short-dc.cfg

  ./pflotran -input_prefix calcite-short


NOTE: this is the same PFloTran input file as
calcite-short-pc!


Volume fraction updates and PFloTran native constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Test the alquimia interface with mineral dissolution/precipitation for
a long simulation time, so that the mineral volume fractions are
updated during reaction stepping. Uses PFloTran native constraints.

**Status: fails**

Status Notes: slight initial numerical differences in rates accumulate error?

.. code-block:: bash

  ./batch_chem -d -i calcite-vf-pc.cfg

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

    ./batch_chem -d -i isotherms-pc.cfg
    ./pflotran -input_prefix isotherms

NOTE: this is the same PFloTran input file as isotherms-ac.

Alquimia supplied constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Status: fails**

Status Notes: final significant figure is off for the "D" species.

.. code-block:: bash

    ./batch_chem -d -i isotherms-ac.cfg
    ./pflotran -input_prefix isotherms

NOTE: this is the same PFloTran input file as isotherms-pc.

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

    ./batch_chem -d -i ion-exchange-valocchi-pc.cfg
    ./pflotran -input_prefix ion-exchange-valocchi

NOTE: this is the same PFloTran input file as ion-exchange-valocchi-pc.

Alquimia supplied constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Status: pass**

Status Notes: 

.. code-block:: bash

    ./batch_chem -d -i ion-exchange-valocchi-ac.cfg
    ./pflotran -input_prefix ion-exchange-valocchi

NOTE: this is the same PFloTran input file as ion-exchange-valocchi-pc.


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

    ./batch_chem -d -i ion-exchange-valocchi-pc.cfg
    ./pflotran -input_prefix ion-exchange-valocchi

NOTE: this is the same PFloTran input file as ion-exchange-valocchi-pc.

Alquimia supplied constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Status: fail**

Status Notes: numerical differences in trailing significant figures

.. code-block:: bash

    ./batch_chem -d -i ion-exchange-valocchi-ac.cfg
    ./pflotran -input_prefix ion-exchange-valocchi

NOTE: this is the same PFloTran input file as ion-exchange-valocchi-pc.

