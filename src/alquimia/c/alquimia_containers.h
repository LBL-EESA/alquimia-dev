/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_CONTAINERS_H_
#define ALQUIMIA_CONTAINERS_H_

/*******************************************************************************
 **
 ** C implementation of the alquimia containers.
 **
 ** These are passed directly into the fortran routines. The
 ** signatures must match exactly with the fortran side of things.
 **
 ******************************************************************************/

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  struct AlquimiaSizes {
    int num_primary;
    int num_kinetic_minerals;
    int num_aqueous_complexes;
    int num_surface_sites;
    int num_ion_exchange_sites;
  };
  
  struct AlquimiaState {
    double water_density;  // [kg/m^3]
    double saturation;  // [-]
    double porosity;  // [-]
    double temperature;  // [celsius]
    double aqueous_pressure; // [Pa]
    double* total_primary;  // [molarity]
    double* total_sorbed;  // [moles/m^3 bulk]
    double* free_ion;  // [molality]
    double* mineral_volume_fraction;  // [-]
    double* mineral_specific_surface_area; // [m^2 mineral/m^3 bulk]
    double* cation_exchange_capacity;  // [moles/m^3 bulk]
    double* surface_site_density;  // [moles/m^3 bulk]
  };
  
  struct AlquimiaMaterialProperties {
    double volume;  // [m^3]
    double* isotherm_kd;  // [?]
    double* freundlich_n; // [?]
    double* langmuir_b;  // [?]
  };
  
  struct AlquimiaAuxiliaryData {
    double* primary_activity_coeff;  // [-]
    double* secondary_activity_coeff;  // [-]
    double* ion_exchange_ref_cation_conc;  // [?]
    double* surface_complex_free_site_conc;  // [?]
  };
  
  struct AlquimiaEngineStatus {
    int error;
    char* message;
    bool converged;
    int num_rhs_evaluations;
    int num_jacobian_evaluations;
    int num_newton_iterations;
  };
  
  struct AlquimiaMetaData {
    bool thread_safe;
    bool temperature_dependent;
    bool pressure_dependent;
    bool porosity_update;
    bool operator_splitting;
    bool global_implicit;
    int index_base;
    int* primary_indices;
    char** primary_names;
    //char** auxiliary_output_names;
  };
  
  struct AlquimiaGeochemicalConstraint {
    char* primary_species;
    char* constraint_type;
    char* associated_species;
    double value;
  };
  
  /* A geochemical condition is an array of geochemical constraints */
  /* How is this going to work in the C/Fortran interface? */
  struct AlquimiaGeochemicalCondition {
    char* name;
    int num_constraints;
    struct AlquimiaGeochemicalConstraint* constraints;
  };

  struct AlquimiaGeochemicalConditionList {
    int num_conditions;
    struct AlquimiaGeochemicalCondition* conditions;
  };
  
  struct AlquimiaOutputData {
    double pH;
    double* mineral_saturation_index;  // [-]
    double* mineral_reaction_rate;  // [?]
  };

  struct AlquimiaInterface {
    /* read data files/structures, initialize memory, basis management
       (includes reading database, swapping basis, etc.) */
    void (*Setup)(
        const char* input_filename,
        void* pft_engine_state,
        struct AlquimiaSizes* sizes,
        struct AlquimiaEngineStatus* status);

    /* gracefully shutdown the engine, cleanup memory */
    void (*Shutdown)(
      void* pft_engine_state,
      struct AlquimiaEngineStatus* status);

    /* constrain processing for boundary/initial constraints. Called
       once for each IC/BC. */
    void (*ProcessCondition)(
        void* pft_engine_state,
        struct AlquimiaGeochemicalCondition* condition,
        struct AlquimiaMaterialProperties* material_props,
        struct AlquimiaState* state,
        struct AlquimiaAuxiliaryData* aux_data,
        struct AlquimiaEngineStatus* status);

    /* take one (or more?) reaction steps in operator split mode */
    void (*ReactionStepOperatorSplit)(
        void* pft_engine_state,
        double* delta_t,
        struct AlquimiaMaterialProperties* material_props,
        struct AlquimiaState* state,
        struct AlquimiaAuxiliaryData* aux_data,
        struct AlquimiaEngineStatus* status);
    
    /* Access to user selected geochemical data for output, i.e. pH, 
       mineral SI, reaction rates */
    void (*GetAuxiliaryOutput)(
        void* pft_engine_state,
        struct AlquimiaEngineStatus* status);
    
    void (*GetEngineMetaData)(
        void* pft_engine_state,
        struct AlquimiaSizes* sizes,
        struct AlquimiaMetaData* meta_data,
        struct AlquimiaEngineStatus* status);
    
    void (*GetPrimaryNameFromIndex)(
        void* pft_engine_state,
        int* primary_index,
        char* primary_name,
        struct AlquimiaEngineStatus* status);

    /* internal representation of the chemistry engine's state */
    void* engine_state;

  };



#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  // ALQUIMIA_CONTAINERS_H_
