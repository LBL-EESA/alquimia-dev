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


Alquimia Functions
==================

Details of the engine independent API and wrapper interface.

Alquimia is intended for mixed language environments: a C/C++ based
driver can call a FORTRAN based engine or vice versa. The common
language for mixed language programming is pure C. All function
parameters must be either:

* plain old data (POD), i.e. int, float, bool, char, c pointers.

* structures containing POD. Nested structures are OK. The ordering
  and size of data contained in the structures is important. It must
  be the same for C and FORTRAN. All structures must defined on the
  :doc:`structures <APIv1_structures>` page.

* strings should be c style strings, i.e. null terminated arrays of characters.

* FORTRAN expects all function parameters to be pass by value. All C
  function parameters must be pointers.

* Because of name mangling, interfaces contained in C++ classes or
  FORTRAN modules must also provide wrapper functions.

All function interfaces are described in a language independent
pseudo-code. Implementations should be described in header files
``XXX_alquimia_interface.h`` or module files
``XXX_alquimia_interface.F90`` and compiled into ``libalquimia_c.a``
or ``libalquimia_fortran.a`` respectively.

NOTE: the **driver** is responsible for allocating the memory used in the data transfer containers!

Alquimia: Setup
~~~~~~~~~~~~~~~

Read data files/structures, initialize memory, basis management,
including reading database, swapping basis, etc.

.. code-block:: none

    void AlquimiaSetup(
        in: engine_native_input_filename <string>,
        in: hands_off <boolean>,
        output: engine_internal_state <void pointer>,
        output: sizes <struct: Alquimia Sizes>,
	output: functionality <struct: Alquimia Engine Functionality>
	output: status <struct: Alquimia Status>)

The engine's native input file contains input that is used to describe 
the chemical properties of the system of interest. 

The hands_off parameter determines whether input should be read directly 
from the native input file, or whether it will be provided by the caller to 
Alquimia. If hands_off is true, the data in AlquimiaProperties will be ignored, 
and only state variables will be copied between the caller and the 
engine. If false, the engine will use all of the data supplied by the caller 
in the Alquimia containers.

The sizes structure contains the number of degrees of freedom for the
mobile phase (aqueous components) and immobile phases (sorbed,
biomass), auxiliary memory requirements. The driver should use this
information for allocating memory and input file validation.

One of the output variables from Setup is a void pointer to the
engines internal representation of the reaction network. **The actual
structure will be different for different engines. The client driver
should store this and return it on all subsequent alquimia function
calls, but should not use it for anything else. Do NOT rely on any
particular form of the engine's internal state.** If the engine is
thread safe (no internal global variables), then OpenMP threading can
be used by calling setup once for each thread and storing multiple
copies of the internal state.

NOTE: PETScInitialize() or MPI_Init() must be called prior to calling
AlquimiaSetup()! PFloTran requires PETSc to be present (headers for
variable types, e.g. PetscReal, and library for logging). For now we
are assuming that it can operate on MPI_COMM_WORLD, but if this turns
out not to be a valid assumption, we will need to pass the specific
communicator to PFloTran through the setup interface.

Alquimia: Shutdown
~~~~~~~~~~~~~~~~~~

Shutdown the engine. Destroys all internal objects, frees manually
allocated memory, and any other shutdown the engine may need.

.. code-block:: none

    void AlquimiaShutdown(
        in/out: engine_internal_state <void pointer>,
	output: status <struct: Alquimia Status>)


Alquimia: Get Problem Meta Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Additional information the engine may need about the specific problem
being solved. For example, the names of the primary species and
minerals, and their required order.

.. code-block:: none

    void AlquimiaGetEngineMetaData(
        in/out: engine_internal_state <void pointer>,
        output: problem_meta_data <struct: Alquimia Problem Meta Data>,
	output: status <struct: Alquimia Status>)



Alquimia: Geochemical Condition Processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Geochemical condition processing for boundary/initial conditions. Called once for each IC/BC, or once for each grid cell for heterogeneous conditions.

.. code-block:: none

    void AlquimiaProcessCondition(
        in/out engine_internal_state <void pointer>,
        input: condition <struct: Alquimia Geochemical Condition>,
        input: properties <struct: Alquimia Properties>,
        in/out: state <struct: Alquimia State>,
        in/out: aux_data <struct: Alquimia Auxiliary Data>,
        output: status <struct: Alquimia Status>)


If the name field of the condition structure is specified and the constraint list is empty, then the engine will check for a condition with the same name in its native input file format. 

Alquimia: Operator Splitting Reaction Step
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Take one reaction step in operator split mode with the specified delta t.

.. code-block:: none

    void AlquimiaReactionStepOperatorSplit(
        in/out: engine_internal_state <void pointer>,
        input: delta_t <double>,
        input: properties <struct: Alquimia Properties>,
        in/out: state <struct: Alquimia State>,
        in/out: aux_data <struct: Alquimia Auxiliary Data>
        out: status <struct: Alquimia Status>)


Alquimia: Get Auxiliary Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Access to user selected geochemical data for output, i.e. pH, mineral SI, reaction rates.

NOTE: as currently implemented in batch mode, this **MUST** be done
immediately after each reaction step call....

:ref:`AlquimiaAuxiliaryOutputData` contains a series of arrays for
different data types. If the driver does not want a particular set of
data, it should set the array size to zero. The engine should use the
value contained in aux_output to determine how much data to write.

.. code-block:: none

    void AlquimiaGetAuxiliaryOutput(
        in/out: engine_internal_state <void pointer>,
        input: state <struct: Alquimia State>,
        input: aux_data <struct: Alquimia Auxiliary Data>,
        output: aux_output <struct: Alquimia Auxiliary Output Data>,
        output: states <struct: Alquimia Status>)


Alquimia: Global Implicit Reaction Step
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is not implemented in version 1.0.0 of the alquimia API. It can be added as a backward compatible feature in version 1.x. 

NOTE: need to keep track of whether driver and engine are using row-major or column-major ordering....

Return the function evaluation and Jacobian information for a GI step

.. code-block:: none

    void AlquimiaReactionStepGlobalImplicit(....)

