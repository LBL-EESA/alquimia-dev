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

#ifdef WINDOWS
#include "xstdbool.h"
#else
#include <stdbool.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  struct AlquimiaVectorDouble {
    int size;
    double* data;
  };

  struct AlquimiaVectorInt {
    int size;
    int* data;
  };

  struct AlquimiaVectorString {
    /* NOTE: this is a vector of strings */
    int size;
    char** data;
  };

  struct AlquimiaSizes {
    int num_primary;
    int num_sorbed;
    int num_kinetic_minerals;
    int num_aqueous_complexes;
    int num_surface_sites;
    int num_ion_exchange_sites;
    int num_isotherm_species;
    int num_aux_integers;
    int num_aux_doubles;
  };
  
  struct AlquimiaState {
    double water_density;  /* [kg/m^3] */
    double saturation;  /* [-] */
    double porosity;  /* [-] */
    double temperature;  /* [celsius] */
    double aqueous_pressure; /* [Pa] */
    struct AlquimiaVectorDouble total_mobile;  /* [molarity] */
    struct AlquimiaVectorDouble total_immobile;  /* [moles/m^3 bulk] */
    struct AlquimiaVectorDouble mineral_volume_fraction;  /* [-] */
    struct AlquimiaVectorDouble mineral_specific_surface_area; /* [m^2 mineral/m^3 bulk] */
    struct AlquimiaVectorDouble surface_site_density;  /* [moles/m^3 bulk] */
    struct AlquimiaVectorDouble cation_exchange_capacity;  /* [moles/m^3 bulk] */
  };
  
  struct AlquimiaMaterialProperties {
    double volume;  /* [m^3] */
    struct AlquimiaVectorDouble isotherm_kd;  /* [?] */
    struct AlquimiaVectorDouble freundlich_n; /* [?] */
    struct AlquimiaVectorDouble langmuir_b;  /* [?] */
  };
  
  struct AlquimiaAuxiliaryData {
    struct AlquimiaVectorInt aux_ints;  /* [-] */
    struct AlquimiaVectorDouble aux_doubles;  /* [-] */
  };
  
  struct AlquimiaEngineStatus {
    int error;
    char* message;
    bool converged;
    int num_rhs_evaluations;
    int num_jacobian_evaluations;
    int num_newton_iterations;
  };
  
  struct AlquimiaEngineFunctionality {
    bool thread_safe;
    bool temperature_dependent;
    bool pressure_dependent;
    bool porosity_update;
    bool operator_splitting;
    bool global_implicit;
    int index_base;
  };
  
  struct AlquimiaProblemMetaData {
    struct AlquimiaVectorString primary_names;
    struct AlquimiaVectorString mineral_names;
    struct AlquimiaVectorString surface_site_names;
    struct AlquimiaVectorString ion_exchange_names;
    struct AlquimiaVectorString isotherm_species_names;
    /*char** auxiliary_output_names;*/
  };
  
  struct AlquimiaAuxiliaryOutputData {
    double pH;
    struct AlquimiaVectorDouble mineral_saturation_index;  /* [-] */
    struct AlquimiaVectorDouble mineral_reaction_rate;  /* [?] */
  };

  /* 
  ** Geochemical Conditions
  */

  struct AlquimiaAqueousConstraint {
    char* primary_species_name;
    char* constraint_type;
    char* associated_species;
    double value;
  };

  struct AlquimiaAqueousConstraintVector {
    int size;
    struct AlquimiaAqueousConstraint* data;
  };

  struct AlquimiaMineralConstraint {
    char* mineral_name;
    double volume_fraction;
    double specific_surface_area;
  };
  
  struct AlquimiaMineralConstraintVector {
    int size;
    struct AlquimiaMineralConstraint* data;
  };

  /* A geochemical condition is an array of aqueous and mineral geochemical constraints */
  struct AlquimiaGeochemicalCondition {
    char* name;
    struct AlquimiaAqueousConstraintVector aqueous_constraints;
    struct AlquimiaMineralConstraintVector mineral_constraints;
  };

  struct AlquimiaGeochemicalConditionVector {
    int size;
    struct AlquimiaGeochemicalCondition* data;
  };
  
#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* ALQUIMIA_CONTAINERS_H_ */
