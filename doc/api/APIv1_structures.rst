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
   
   Authors: Benjamin Andre <bandre@lbl.gov>, Sergi Molins <smolins@lbl.gov>, 
   Jeff Johnson <jnjohnson@lbl.gov>


Alquimia Data Transfer Containers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These are the containers that alquimia will used to pass data between
the driver code and engine. The **driver** is responsible for
allocating and freeing the memory used in the data transfer
containers!


**Implementation notes**:
 
* In order to pass data between C/C++ and FORTRAN, all containers
  should be **implemented** as structures or derived types of plain old
  data (integers, doubles, characters, bool/logical) and pointers to
  arrays.
* Structures may be nested.
* The ordering and size of the data contained in the structures is
  important. It must be the same on both sides of the C/FORTRAN
  interface.
* Check your units.


The implementations of these structures are contained in
``alquimia_containers.h`` and ``alquimia_containers.F90`` for C/C++
and FORTRAN respectively.

Struct: Alquimia Vectors
========================

All arrays used in Alquimia will be AlquimiaVectors, grouping the size
and data into a single structure. AlquimiaVectorTYPE is denoted in the
following tables as "vector<type, size>".

+--------------+---------------+
| **variable** | **storage**   |
+==============+===============+
| size         | int           |
+--------------+---------------+
| data         | array<type>   |
+--------------+---------------+


Stand alone variables
=====================

+--------------+--------------+----------------------------------+
| **variable** | **storage**  |**use**                           |
+==============+==============+==================================+
| engine_state | void pointer |pointer to a struct/derived type  |
|              |              |containing all the persistent     |
|              |              |internal state data for the       |
|              |              |chemistry engine.                 |
+--------------+--------------+----------------------------------+

NOTE: The engine state variable should be used by the engine to store
any spatially independent, persistent data. All spatially or mesh
dependent must be stored in the auxiliary data structure managed by
the driver. This allows an AMR based driver to refine or interpolate
as needed.

Struct: Alquimia Sizes
======================

The memory requirements and size of the various arrays that are being
passed through the interface.

+-------------------------+-------------+---------------------------------------------------------+
| **variable**            | **storage** | **use**                                                 |
+=========================+=============+=========================================================+
| num_primary             | int         | N_p, number of primary species                          |
+-------------------------+-------------+---------------------------------------------------------+
| num_sorbed              | int         | N_sorb, number of sorbed species [3]_                   |
+-------------------------+-------------+---------------------------------------------------------+
| num_minerals            | int         | N_m, number of minerals                                 |
+-------------------------+-------------+---------------------------------------------------------+
| num_aqueous_complexes   | int         |       N_s, number of secondary aqueous complexes        |
+-------------------------+-------------+---------------------------------------------------------+
| num_aqueous_kinetics    | int         | N_k, number of kinetic aqueous reactions                |
+-------------------------+-------------+---------------------------------------------------------+
| num_surface_sites       | int         | N_ss, number of surface sites                           |
+-------------------------+-------------+---------------------------------------------------------+
| num_ion_exchange_sites  | int         | N_ix, number of ion exchange sites                      |
+-------------------------+-------------+---------------------------------------------------------+
| num_isotherm_species    | int         | N_is, number of species involved in isotherm reactions  |
+-------------------------+-------------+---------------------------------------------------------+
| num_aux_integers        | int         | N_ai, number of auxiliary integers                      |
+-------------------------+-------------+---------------------------------------------------------+
| num_aux_doubles         | int         | N_ad, number of auxiliary doubles                       |
+-------------------------+-------------+---------------------------------------------------------+

.. [3] Consider the number of sorbed species to be a flag requiring memory allocation for the immobile component of the primary species, generally assert (N_sorb == N_p || N_sorb == 0).



Struct: Alquimia State
======================

Storage for spatially and temporally varying "state" data. Read/write (chemistry may change these values).

+-----------------------------------+------------------------+-----------------------------+
| **variable**                      |      **storage**       |        **units**            |
+===================================+========================+=============================+
| water_density                     | double                 | [kg/m^3]                    |
+-----------------------------------+------------------------+-----------------------------+
| porosity                          | double                 | [m^3 pore space / m^3 bulk] |
+-----------------------------------+------------------------+-----------------------------+
| temperature                       | double                 | [deg C]                     |
+-----------------------------------+------------------------+-----------------------------+
| aqueous_pressure                  | double                 | [Pa]                        |
+-----------------------------------+------------------------+-----------------------------+
| total_mobile                      | vector<double, N_p>    | [molarity: moles/L]         |
+-----------------------------------+------------------------+-----------------------------+
| total_immobile                    | vector<double, N_sorb> | [moles/m^3 bulk]            |
+-----------------------------------+------------------------+-----------------------------+
| mineral_volume_fractions          | vector<double, N_m>    | [-]                         |
+-----------------------------------+------------------------+-----------------------------+
| mineral_specific_surface area     | vector<double, N_m>    | [m^2 mineral / m^3 bulk]    |
+-----------------------------------+------------------------+-----------------------------+
| surface_site_density              | vector<double, N_ss>   | [moles / m^3 bulk]          |
+-----------------------------------+------------------------+-----------------------------+
| cation_exchange_capacity          | vector<double, N_ix>   | [moles / m^3 bulk]          |
+-----------------------------------+------------------------+-----------------------------+


Struct: Alquimia Properties
===========================

Storage for spatially variable "parameters", not changing in time. Read only (chemistry may not change).

+---------------------------+-----------------------+-------------------------------+
| **variable**              |      **storage**      | **units**                     |
+===========================+=======================+===============================+
| volume                    | double                | [m^3]                         |
+---------------------------+-----------------------+-------------------------------+
| saturation                | double                | [m^3 liquid / m^3 pore space] |
+---------------------------+-----------------------+-------------------------------+
| Isotherm Kd               | vector<double, N_is>  | [kg H20 / m^3 bulk]           |
+---------------------------+-----------------------+-------------------------------+
| Freundlich N              | vector<double, N_is>  | [-]                           |
+---------------------------+-----------------------+-------------------------------+
| Langmuir b                | vector<double, N_is>  | [-]                           |
+---------------------------+-----------------------+-------------------------------+
| mineral_rate_cnst         | vector<double, N_m>   | [mol /m^2-sec]                |
+---------------------------+-----------------------+-------------------------------+
| aqueous_kinetic_rate_cnst | vector<double, N_k>   | [sec^-1] [5]_                 |
+---------------------------+-----------------------+-------------------------------+

.. [5] Units will vary with the reaction model in the engine; sec^-1 corresponds to a first order rate dependence.

Struct: Alquimia Auxiliary Data
===============================

Spatially and temporally variable internal state data for the
chemistry engine. The driver code **MUST** store these on the mesh
(cell by cell basis). They must be passed to the engine at the start
of each reaction step call and stored after the call.  This data
**MUST** be written to checkpoint/restart files. The driver does not
need to do anything else with them, they do not need to be
transported. Persistent data that is not mesh dependent should be
stored in the engine state variable.

This structure is intended for things like free ion concentrations,
primary and secondary activity coefficients, surface complex free site
concentration [1]_, ion exchange reference cation concentration [1]_,
etc. The engine is responsible for packing and unpacking this data as
needed.

+----------------+-----------------------+------------+
| **variable**   | **storage**           | **units**  |
+================+=======================+============+
| aux_ints       | vector<int, N_ai>     | [-]        |
+----------------+-----------------------+------------+
| aux_double     | vector<double, N_ad>  | [-]        |
+----------------+-----------------------+------------+


.. [1] PFloTran internal variable that must be stored



Struct: Alquimia Engine Status
==============================

Return the status of the geochemistry engine after the last
operation.

+--------------------------+-------------+
| **variable**             | **storage** |
+==========================+=============+
| error                    | int         |
+--------------------------+-------------+
| message                  | string      |
+--------------------------+-------------+
| converged                | bool        |
+--------------------------+-------------+
| num_rhs_evaluations      | int         |
+--------------------------+-------------+
| num_jacobian_evaluations | int         |
+--------------------------+-------------+
| num_newton_iterations    | int         |
+--------------------------+-------------+

* Every alquimia function call should be followed by a check
  of the error status. 

* Convergence failure is a normal part of numerical computing, **NOT**
  an error.

* error messages in the message string should spell out the source of
  the error as much as possible. Developer errors should be
  distinguished from user errors if possible. Use something like
  "DEV_ERROR:" or "INPUT_ERROR:" at the start of the string.


Struct: Alquimia Engine Functionality
=====================================

Information about the functionality supported by the geochemistry
engine. This is **not** necessarily a hard coded list. For example,
the engine may support temperature dependent chemistry for a
particular problem only if the user supplied database contains the
appropriate data.

+-------------------------+---------------------+-------------------------------------------+
| **variable**            | **storage**         |**comment**                                |
+=========================+=====================+===========================================+
| thread safe             | bool                |tells the client whether it can create     |
|                         |                     |multiple copies of the chemistry engine on |
|                         |                     |the same processor and farm out work using |
|                         |                     |OpenMP or something similar. Only valid if |
|                         |                     |the engine doesn't have global variables.  |
+-------------------------+---------------------+-------------------------------------------+
| temperature dependent   | bool                |Engine supports temperature dependent      |
|                         |                     |chemistry                                  |
+-------------------------+---------------------+-------------------------------------------+
| pressure dependent      | bool                |Engine supports pressure dependent         |
|                         |                     |chemistry                                  |
+-------------------------+---------------------+-------------------------------------------+
| porosity updates        | bool                |Engine supports porosity updates due to    |
|                         |                     |mineral dissolution/precipitation, biomass |
|                         |                     |clogging, etc.                             |
+-------------------------+---------------------+-------------------------------------------+
| operator splitting      | bool                |Engine supports operator splitting reaction|
|                         |                     |stepping                                   |
+-------------------------+---------------------+-------------------------------------------+
| global implicit         | bool                |Engine supports global implicit reaction   |
|                         |                     |stepping                                   |
+-------------------------+---------------------+-------------------------------------------+
| base index              | int                 |base index for vectors passed between the  |
|                         |                     |driver and engine i.e. if the engine is C  |
|                         |                     |based, base index = 0, if FORTRAN based,   |
|                         |                     |base index = 1                             |
+-------------------------+---------------------+-------------------------------------------+

Struct: Alquimia Problem Meta Data
==================================

Problem specific meta data, e.g. primary species and mineral
names. Species are in the order that the chemistry engine expects to
receive data.

+------------------------+----------------------+--------------------------------------------+
| **variable**           | **storage**          | **comment**                                |
+========================+======================+============================================+
| primary_names          | vector<string, N_p>  | names of the primary species               |
+------------------------+----------------------+--------------------------------------------+
| positivity             | vector<int, N_p>     | positivity of the primary species (1 or 0) |
+------------------------+----------------------+--------------------------------------------+
| mineral_names  | vector<string, N_m>  | names of the minerals              |
+------------------------+----------------------+--------------------------------------------+
| surface_site_names     | vector<string, N_ss> | names of the surface sites                 |
+------------------------+----------------------+--------------------------------------------+
| ion_exchange_names     | vector<string, N_ix> | names of the ion exchange sites            |
+------------------------+----------------------+--------------------------------------------+
| isotherm_species_names | vector<string, N_is> | names of the primary species involved in   |
|                        |                      | isotherm reactions                         |
+------------------------+----------------------+--------------------------------------------+
| aqueous_kinetic_names  | vector<string, N_k>  | names of the kinetic aqueous reactions     |
+------------------------+----------------------+--------------------------------------------+

The positivity array is the same size as primary_names, and its ith entry 
contains 1 if the ith primary species must be positive, 0 if it has no 
such positivity constraint.

.. _AlquimiaAuxiliaryOutputData:

Struct: Alquimia Auxiliary Output Data
======================================

Additional data that the user may request be written to the output
files. The engine ignores any value passed in with these arrays and
over writes it with the current value. If the driver does not want
data in a particular array, it should set the size to zero.

+----------------------------------+------------------------+------------------------+
|       **variable**               |        **type**        |       **units**        |
+==================================+========================+========================+
| pH                               |         double         | [-]                    |
+----------------------------------+------------------------+------------------------+
| aqueous_kinetic_rate             |  vector<double, N_k>   | [mol / src / m^3]      |
+----------------------------------+------------------------+------------------------+
| mineral_saturation_index         |  vector<double, N_m>   | [-]                    |
+----------------------------------+------------------------+------------------------+
| mineral_reaction_rate            |  vector<double, N_m>   | [mol / sec / m^3 bulk] |
+----------------------------------+------------------------+------------------------+
| primary_free_ion_concentration   |  vector<double, N_p>   | [molality: mol/kg H2O] |
+----------------------------------+------------------------+------------------------+
| primary_activity_coeff           |  vector<double, N_p>   | [-]                    |
+----------------------------------+------------------------+------------------------+
| secondary_free_ion_concentration |  vector<double, N_s>   | [molality: mol/kg H2O] |
+----------------------------------+------------------------+------------------------+
| secondary_activity_coeff         |  vector<double, N_s>   | [-]                    |
+----------------------------------+------------------------+------------------------+


TODO(bja): to keep things simple, we just write out all the mineral
data. If the driver only wants a subset, then they can grab the ones
they want using the name-index mapping provided by the problem meta
data.... 

TODO(bja): this is only considering kinetic minerals. User may want
reference minerals as well....

Struct: Alquimia Geochemical Condition
======================================

Geochemical Condition is a structure containing a name string and a
vector of geochemical constraints. There must be one constraint for
each primary species and each kinetic mineral.

+---------------------+---------------------------------+
|    **variable**     |            **type**             |
+=====================+=================================+
|        name         |             string              |
+---------------------+---------------------------------+
| aqueous_constraints | vector<aqueous_constraint, N_p> |
+---------------------+---------------------------------+
| mineral_constraints | vector<mineral_constraint, N_m> |
+---------------------+---------------------------------+


Struct: Alquimia Aqueous Constraint
===================================

An aqueous geochemical constraint is a structure with the following fields:

+--------------------+----------+
| **variable**       | **type** |
+====================+==========+
| primary species    | string   |
+--------------------+----------+
| constraint type    | string   |
+--------------------+----------+
| associated species | string   |
+--------------------+----------+
| value              | double   |
+--------------------+----------+

"Associated species" is the name of the mineral or gas associated with
that constraint, e.g. Ca++ is constrained by equilibrium with the
mineral calcite or HCO3- is constrained by equilibrium with CO2 gas.

Types of constraints supported:

* total_aqueous
* total_sorb
* free
* mineral
* gas
* pH
* charge

These are named alquimia string constants, :ref:`AlquimiaStrings`.

If an engine does not support a particular type of constraint, it
should report an error.

The units for a constraint value depend on the constraint type, and
should agree with the units defined above, e.g. total_aqueous should
agree with total_mobile from AlquimiaState, free ion concentration
should agree with free ion units from AquimiaAuxOutput.

If a constraint type does not require a supplied value, e.g. charge,
then the user/driver should supply either a safe initial guess (1.0e-9 for
example) that the engine can use, or a very small non-zero value
(1.0e-20). The engine may use this or chose to ignore it.

Struct: Alquimia Mineral Constraint
===================================

A mineral geochemical constraint is a structure with the following fields:

+---------------------+----------+---------------------------+
| **variable**        | **type** |         **units**         |
+=====================+==========+===========================+
| mineral_name        | string   | [-]                       |
+---------------------+----------+---------------------------+
| volume_fraction     | double   | [-]                       |
+---------------------+----------+---------------------------+
|specific_surface_area| double   | [m^2 mineral / m^3 bulk]  |
+---------------------+----------+---------------------------+
