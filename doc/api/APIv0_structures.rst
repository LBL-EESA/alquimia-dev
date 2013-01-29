Alquimia Data Transfer Containers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These are the containers that will be used to pass data between the
transport code and alquimia. The alquimia may pass these on to the
chemistry engine directly, or may change units, container type etc and
pass the munged data on to the engine.

**Implementation note**: In order to pass data between C/C++ and
Fortran, all containers should be **implemented** as structs or
derived types of plain old data (integers, doubles, characters,
bool/logical) and pointers to arrays. Structs may be nested.



Struct: Alquimia Vectors
========================

All arrays used in Alquimia will be AlquimiaVectors, grouping the size
and data into a single struct. AlquimiaVectorTYPE is denoted in the
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


Struct: Alquimia Sizes
======================

The size of the various arrays that are being passed through the interface.

+-------------------------+-------------+--------------------------------------------+
| **variable**            | **storage** | **use**                                    |
+=========================+=============+============================================+
| num_primary             | int         | N_p, number of primary species             |
+-------------------------+-------------+--------------------------------------------+
| num_kinetic_minerals    | int         | N_m, number of kinetic minerals            |
+-------------------------+-------------+--------------------------------------------+
| num_aqueous_complexes   | int         | N_s, number of secondary aqueous complexes |
+-------------------------+-------------+--------------------------------------------+
| num_surface_sites       | int         | N_ss, number of surface sites              |
+-------------------------+-------------+--------------------------------------------+
| num_ion_exchange_sites  | int         | N_ix, number of ion exchange sites         |
+-------------------------+-------------+--------------------------------------------+

Struct: Alquimia State
======================

Storage for spatially and temporally varying "state" data. Read/write (chemistry may change these values).

+-----------------------------------+----------------------+-------------------------+
| **variable**                      |     **storange**     |        **units**        |
+===================================+======================+=========================+
| saturation                        |        double        |           [-]           |
+-----------------------------------+----------------------+-------------------------+
| porosity                          |        double        |           [-]           |
+-----------------------------------+----------------------+-------------------------+
| temperature                       |        double        |           [?]           |
+-----------------------------------+----------------------+-------------------------+
| aqueous_pressure                  |        double        |           [-]           |
+-----------------------------------+----------------------+-------------------------+
| total_primary                     | vector<double, N_p>  |       [molarity]        |
+-----------------------------------+----------------------+-------------------------+
| total_sorbed                      | vector<double, N_p>  |    [moles/m^3 bulk]     |
+-----------------------------------+----------------------+-------------------------+
| free_ion                          | vector<double, N_p>  |       [molality]        |
+-----------------------------------+----------------------+-------------------------+
| mineral_volume_fractions          | vector<double, N_m>  |           [-]           |
+-----------------------------------+----------------------+-------------------------+
| mineral_specific_surface area     | vector<double, N_m>  |[m^2 mineral / m^3 bulk] |
+-----------------------------------+----------------------+-------------------------+
| cation_exchange_capacity          | vector<double, N_ix> |           [-]           |
+-----------------------------------+----------------------+-------------------------+
| surface_site_density              | vector<double, N_ss> |      [moles / m^3]      |
+-----------------------------------+----------------------+-------------------------+


Possible changes: 
* move free ion into alquimia auxiliary data?
* rename total_primary => total mobile?
* rename total_sorbed => total immobile?

Struct: Alquimia Material Properties
====================================

Storage for spatially variable "parameters", not changing in time. Read only (chemistry may not change).

+--------------+-----------------------+------------+
| **variable** |      **storage**      | **units**  |
+==============+=======================+============+
| volume       |        double         |   [m^3]    |
+--------------+-----------------------+------------+
| Isotherm Kd  |  vector<double, N_p>  | [?]        |
+--------------+-----------------------+------------+
| Freundlich N |  vector<double, N_p>  | [?]        |
+--------------+-----------------------+------------+
| Langmuir b   |  vector<double, N_p>  | [?]        |
+--------------+-----------------------+------------+

Struct: Alquimia Auxiliary Data
===============================

Spatially and temporally variable internal state data for the
chemistry engine. The driver code **MUST** store these on the mesh
(cell by cell basis). They must be passed to the engine at the start
of each reaction step call and stored after the call.  This data
**MUST** be written to checkpoint/restart files. The driver does not
need to do anything else with them, they do not need to be
transported.

+-------------------------------------+-----------------------+--------------------------+
| **variable**                        | **storage**           | **units**                |
+=====================================+=======================+==========================+
| primary activity coeff              |  vector<double, N_p>  | [-]                      |
+-------------------------------------+-----------------------+--------------------------+
| secondary activity coeff            |  vector<double, N_s>  | [-]                      |
+-------------------------------------+-----------------------+--------------------------+
| ion_exchange_ref_cation_conc [1]_   | vector<double, N_ix>  | [?]                      |
+-------------------------------------+-----------------------+--------------------------+
| surface_complex_free_site_conc [1]_ | vector<double, N_ss>  | [moles sites / m^3 bulk] |
+-------------------------------------+-----------------------+--------------------------+

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
| message                  |   string    |
+--------------------------+-------------+
| converged                | boo         |
+--------------------------+-------------+
| num_rhs_evaluations      | int         |
+--------------------------+-------------+
| num_jacobian_evaluations | int         |
+--------------------------+-------------+
| num_newton_iterations    | int         |
+--------------------------+-------------+

* Every alquimia function call should be followed by a check
  of the error status. 

* Convergence failure is **NOT** an error.

* error messages in the message string should spell out the source of
  the error as much as possible. Developer errors should be
  distinguished from user errors if possible. Use something like
  "DEV_ERROR:" at the start of the string.


Struct: Alquimia Meta Data
==========================

Other information exchanged between the engine and client

+-------------------------+---------------------+-------------------------------------------+
| **variable**            | **storage**         |**comment**                                |
+=========================+=====================+===========================================+
| thread safe             | bool                |tells the client whether it can create     |
|                         |                     |multiple copies of the chemistry engine on |
|                         |                     |the same processor and farm out work using |
|                         |                     |openmp or something similar. Only valid if |
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
|                         |                     |based, base index = 0, if fortran based,   |
|                         |                     |base index = 1                             |
+-------------------------+---------------------+-------------------------------------------+
| primary indices         | vector<int, N_p>    |indices of the named primaries             |
+-------------------------+---------------------+-------------------------------------------+
| primary names           | vector<string, N_p> |names of the primary species               |
+-------------------------+---------------------+-------------------------------------------+
| kinetic mineral indices | vector<string, N_m> |indices of the kinetic minerals            |
+-------------------------+---------------------+-------------------------------------------+
| kinetic mineral names   | vector<string, N_m> |names of the kinetic minerals              |
+-------------------------+---------------------+-------------------------------------------+

Struct: Alquimia Auxiliary Output Data
======================================

Additional data that the user may request be written to the output
files. The engine ignores any value passed in with these arrays and
over writes it with the current value.

+--------------------------+------------------------+-----------+
|       **variable**       |        **type**        | **units** |
+==========================+========================+===========+
| pH                       |         double         | [-]       |
+--------------------------+------------------------+-----------+
| mineral_saturation_index |  vector<double, N_m>   | [-]       |
+--------------------------+------------------------+-----------+
|  mineral_reaction_rate   |  vector<double, N_m>   | [?]       |
+--------------------------+------------------------+-----------+

Struct: Alquimia Geochemical Condition
======================================

Geochemical Condition is a struct containing a name string and a vector of geochemical constraints. There must be one constraint for each primary species.

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

An aqueous geochemical constraint is a struct with the following fields:

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

Types of constraints supported:
* mineral
* gas
* pH
* charge

"Associated species" is the name of the mineral or gas associated with
that constraint, e.g. Ca++ is constrained by equilibrium with the
mineral calcite or HCO3- is constrained by equilibrium with CO2 gas.

Struct: Alquimia Mineral Constraint
===================================

A mineral geochemical constraint is a struct with the following fields:

+---------------------+----------+-----------+
| **variable**        | **type** | **units** |
+=====================+==========+===========+
| mineral_name        | string   | [-]       |
+---------------------+----------+-----------+
| volume_fraction     | double   | [-]       |
+---------------------+----------+-----------+
|specific_surface_area| double   | [?]       |
+---------------------+----------+-----------+

