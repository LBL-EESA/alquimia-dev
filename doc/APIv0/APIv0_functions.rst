Details of the engine independent API and wrapper interface.

All function parameters must be either:

* plain old data (POD), i.e. int, float, bool, char, c pointers.

* structures containing POD. All structures must defined on the [[APIv0_Structures | structures]] section.

* strings should be c style strings, i.e. null terminated arrays of characters.

== **Alquimia: Setup** ==
Read data files/structures, initialize memory, basis management, including reading database, swapping basis, etc. 

One of the output variables from Setup is a void pointer to the engines internal representation of the reaction network. //The actual struct will be different for different engines. The client driver should store this and return it on all subsequent alquimia function calls, but should not use it for anything else. Do NOT rely on any particular form of the engine's internal state.// If the engine is thread safe (no internal global variables), then OpenMP threading can be used by calling setup once for each thread and storing multiple copies of the internal state.

{{{
void AlquimiaSetup(
    in/out: engine_native_input_file <string>,
    input: mpi_com <?>,
    output: engine_internal_state <void pointer >,
    output: sizes <struct: Alquimia Sizes>)
}}}

The sizes structure contains the number of degrees of freedom for the mobile phase (aqueous components) and immobile phases (sorbed, biomass), auxiliary memory requirements. The driver should use this information for allocating memory and input file validation.

NOTE: PETScInitialize() or MPI_Init() must be called prior to calling AlquimiaSetup()! PFloTran requires mpi/petsc to be present (TODO: document why). For now we are assuming that it can operate on MPI_COMM_WORLD, but if this turns out not to be a valid assumption, we will need to pass the specific communicator to PFloTran through the setup interface. 

== **Alquimia: Shutdown** ==
Shutdown the engine. Destroys all internal objects, frees manually allocated memory, and any other shutdown the engine may need.
{{{
void AlquimiaShutdown(
    in/out: engine_internal_state <void pointer>)
}}}

== **Alquimia: Get Engine Meta Data** ==
Additional information about the library that may be needed by the engine. Is it tread_safe (can the engine create multiple copies using openmp). Does it support temperature dependent chemistry? etc
{{{
void AlquimiaGetEngineMetaData(
    in/out: engine_internal_state <void pointer>,
    output: engine_meta_data <struct: Alquimia Meta Data>)
}}}


== **Alquimia: Geochemical Condition Processing** ==
Geochemical condition processing for boundary/initial conditions. Called once for each IC/BC, or once for each grid cell for heterogeneous conditions.
{{{
void AlquimiaProcessCondition(
    in/out engine_internal_state <void pointer>,
    input: condition <struct: Alquimia Geochemical Condition>,
    input: material_properties <struct: Alquimia Material Properties>,
    in/out: state <struct: Alquimia State>,
    in/out: aux_data <struct: Alquimia Auxiliary Data>,
    output: status <struct: Alquimia Status>)
}}}

If the name field of the condition structure is specified and the constraint list is empty, then the engine will check for a condition with the same name in its native input file format. 

== **Alquimia: Operator Splitting Reaction Step** ==
Take one (or more?) reaction steps in operator split mode
{{{
void AlquimiaReactionStepOperatorSplit(
    in/out: engine_internal_state <void pointer>,
    input: delta_t <double>,
    input: material_properties <struct: Alquimia Material Properties>,
    in/out: state <struct: Alquimia State>,
    in/out: aux_data <struct: Alquimia Auxiliary Data>
    out: status <struct: Alquimia Status>)
}}}

== **Alquimia: Get Auxiliary Output** ==
Access to user selected geochemical data for output, i.e. pH, mineral SI, reaction rates.

NOTE: as currently implemented in batch mode, this **MUST** be done after each reaction step call....

{{{
void AlquimiaGetAuxiliaryOutput(
    in/out: engine_internal_state <void pointer>,
    input: state <struct: Alquimia State>,
    input: aux_data <struct: Alquimia Auxiliary Data>,
    output: aux_output <struct: Alquimia Auxiliary Output Data>,
    output: states <struct: Alquimia Status>)
}}}

== **Alquimia: Global Implicit Reaction Step** ==

NOTE: This will not be implemented in version 1.0.0 of the alquimia API. By not specifying it in version 1.0, it can be added as a backward compatible feature in version 1.x. 

NOTE: need to keep track of whether driver and engine are using row-major or column-major ordering....

Return the function evaluation and jacobian information for a GI step
{{{
void AlquimiaReactionStepGlobalImplicit(....)
}}}
