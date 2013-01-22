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
    int num_sorbed;
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
    int size_total_primary;
    double* total_primary;  // [molarity]
    int size_total_sorbed;
    double* total_sorbed;  // [moles/m^3 bulk]
    int size_free_ion;
    double* free_ion;  // [molality]
    int size_mineral_volume_fraction;
    double* mineral_volume_fraction;  // [-]
    int size_mineral_specific_surface_area;
    double* mineral_specific_surface_area; // [m^2 mineral/m^3 bulk]
    int size_cation_exchange_capacity;
    double* cation_exchange_capacity;  // [moles/m^3 bulk]
    int size_surface_site_density;
    double* surface_site_density;  // [moles/m^3 bulk]
  };
  
  struct AlquimiaMaterialProperties {
    double volume;  // [m^3]
    int size_isotherm_kd;
    double* isotherm_kd;  // [?]
    int size_freundlich_n;
    double* freundlich_n; // [?]
    int size_langmuir_b;
    double* langmuir_b;  // [?]
  };
  
  struct AlquimiaAuxiliaryData {
    int size_primary_activity_coeff;
    double* primary_activity_coeff;  // [-]
    int size_secondary_activity_coeff;
    double* secondary_activity_coeff;  // [-]
    int size_ion_exchange_ref_cation_conc;
    double* ion_exchange_ref_cation_conc;  // [?]
    int size_surface_complex_free_site_conc;
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
    int size_primary;
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


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  // ALQUIMIA_CONTAINERS_H_
