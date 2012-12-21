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

#define ALQUIMIA_MAX_STRING_LENGTH 256

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

struct AlquimiaSizes_C {
  int num_primary;
  int num_kinetic_minerals;
  int num_aqueous_complexes;
  int num_surface_sites;
  int num_ion_exchange_sites;
};

struct AlquimiaState_C {
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

struct AlquimiaMaterialProperties_C {
  double volume;  // [m^3]
  double* isotherm_kd;  // [?]
  double* freundlich_n; // [?]
  double* langmuir_b;  // [?]
};

struct AlquimiaAuxiliaryData_C {
  double* primary_activity_coeff;  // [-]
  double* secondary_activity_coeff;  // [-]
  double* ion_exchange_ref_cation_conc;  // [?]
  double* surface_complex_free_site_conc;  // [?]
};

struct AlquimiaEngineStatus_C {
  int num_rhs_evaluations;
  int num_jacobian_evaluations;
  int num_newton_iterations;
  bool converged;
};

struct AlquimiaMetaData_C {
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

struct AlquimiaGeochemicalConstraint_C {
  char* primary_species;
  char* constraint_type;
  char* associated_species;
  double value;
};

/* A geochemical condition is an array of geochemical constraints */
/* How is this going to work in the C/Fortran interface? */
struct AlquimiaGeochemicalCondition_C {
  char* name;
  int num_constraints;
  struct AlquimiaGeochemicalConstraint_C* constraints;
};

struct AlquimiaGeochemicalConditionList_C {
  int num_conditions;
  struct AlquimiaGeochemicalCondition_C* conditions;
};

struct AlquimiaOutputData_C {
  double pH;
  double* mineral_saturation_index;  // [-]
  double* mineral_reaction_rate;  // [?]
};

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  // ALQUIMIA_CONTAINERS_H_
