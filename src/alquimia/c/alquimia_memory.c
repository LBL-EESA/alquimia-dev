/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*******************************************************************************
 **
 **  Alquimia C memory utilities to handle memory management
 **
 *******************************************************************************/

#include "alquimia_memory.h"

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "alquimia_containers.h"

/*******************************************************************************
 **
 **  State
 **
 *******************************************************************************/

void AllocateAlquimiaState(const struct AlquimiaSizes_C* sizes,
                           struct AlquimiaState_C* state) {
  state->total_primary = NULL;
  state->total_sorbed = NULL;
  state->free_ion = NULL;
  state->mineral_volume_fraction = NULL;
  state->mineral_specific_surface_area = NULL;
  state->cation_exchange_capacity = NULL;
  state->surface_site_density = NULL;
  //state-> = NULL;

  if (sizes->num_primary > 0) {
    state->total_primary = (double*) calloc(sizes->num_primary, sizeof(double));
    if (NULL == state->total_primary) {
      // TODO(bja): error handling
    }
    state->free_ion = (double*) calloc(sizes->num_primary, sizeof(double));
    if (NULL == state->free_ion) {
      // TODO(bja): error handling
    }
  }

  if (sizes->num_kinetic_minerals > 0) {
    state->mineral_volume_fraction = (double*) calloc(
        sizes->num_kinetic_minerals, sizeof(double));
    if (NULL == state->mineral_volume_fraction) {
      // TODO(bja): error handling
    }
    state->mineral_specific_surface_area = (double*) calloc(
        sizes->num_kinetic_minerals, sizeof(double));
    if (NULL == state->mineral_specific_surface_area) {
      // TODO(bja): error handling
    }
  }
}  // end AllocateAlquimiaState()

void FreeAlquimiaState(struct AlquimiaState_C* state) {
  if (state != NULL) {
    free(state->total_primary);
    free(state->free_ion);
    free(state->mineral_volume_fraction);
    free(state->mineral_specific_surface_area);
  }
  state->total_primary = NULL;
  state->free_ion = NULL;
  state->mineral_volume_fraction = NULL;
  state->mineral_specific_surface_area = NULL;
}  // end FreeAlquimiaState()

/*******************************************************************************
 **
 **  Auxiliary Data
 **
 *******************************************************************************/

void AllocateAlquimiaAuxiliaryData(const struct AlquimiaSizes_C* sizes,
                                   struct AlquimiaAuxiliaryData_C* aux_data) {
  aux_data->primary_activity_coeff = NULL;
  aux_data->secondary_activity_coeff = NULL;
  aux_data->ion_exchange_ref_cation_conc = NULL;
  aux_data->surface_complex_free_site_conc = NULL;
  //aux_data-> = NULL;

  if (sizes->num_primary > 0) {
    aux_data->primary_activity_coeff = (double*) calloc(sizes->num_primary,
                                                        sizeof(double));
    if (NULL == aux_data->primary_activity_coeff) {
      // TODO(bja): error handling
    }
  }

  if (sizes->num_aqueous_complexes > 0) {
    aux_data->secondary_activity_coeff = (double*) calloc(sizes->num_aqueous_complexes,
                                                          sizeof(double));
    if (NULL == aux_data->secondary_activity_coeff) {
      // TODO(bja): error handling
    }
  }

  if (sizes->num_ion_exchange_sites > 0) {
    aux_data->ion_exchange_ref_cation_conc = (double*) calloc(
        sizes->num_ion_exchange_sites, sizeof(double));
    if (NULL == aux_data->ion_exchange_ref_cation_conc) {
      // TODO(bja): error handling
    }
  }
  if (sizes->num_surface_sites > 0) {
    aux_data->surface_complex_free_site_conc = (double*) calloc(
        sizes->num_surface_sites, sizeof(double));
    if (NULL == aux_data->surface_complex_free_site_conc) {
      // TODO(bja): error handling
    }
  }
}  // end AllocateAlquimiaAuxiliaryData()

void FreeAlquimiaAuxiliaryData(struct AlquimiaAuxiliaryData_C* aux_data) {
  if (aux_data != NULL) {
    free(aux_data->primary_activity_coeff);
    free(aux_data->secondary_activity_coeff);
    free(aux_data->ion_exchange_ref_cation_conc);
    free(aux_data->surface_complex_free_site_conc);
  }
  aux_data->primary_activity_coeff = NULL;
  aux_data->secondary_activity_coeff = NULL;
  aux_data->ion_exchange_ref_cation_conc = NULL;
  aux_data->surface_complex_free_site_conc = NULL;
}  // end FreeAlquimiaAuxiliaryData()

/*******************************************************************************
 **
 **  Material Properties
 **
 *******************************************************************************/

void AllocateAlquimiaMaterialProperties(
    const struct AlquimiaSizes_C* sizes,
    struct AlquimiaMaterialProperties_C* material_props) {
  /* NOTE(bja) : need to be smarter about how we allocate memory for
     isotherms. (1) Only allocate if isotherms are used in chemistry, and
     (2) only allocate for the primary species that are being sorbed. */
  material_props->isotherm_kd = NULL;
  material_props->freundlich_n = NULL;
  material_props->langmuir_b = NULL;
  //material_props-> = NULL;

  if (sizes->num_primary > 0) {
    material_props->isotherm_kd = (double*) calloc(sizes->num_primary,
                                                   sizeof(double));
    if (NULL == material_props->isotherm_kd) {
      // TODO(bja): error handling
    }
    material_props->freundlich_n = (double*) calloc(sizes->num_primary,
                                                    sizeof(double));
    if (NULL == material_props->freundlich_n) {
      // TODO(bja): error handling
    }
    material_props->langmuir_b = (double*) calloc(sizes->num_primary,
                                                  sizeof(double));
    if (NULL == material_props->langmuir_b) {
      // TODO(bja): error handling
    }
  }

}  // end AllocateAlquimiaMaterialProperties()

void FreeAlquimiaMaterialProperties(
    struct AlquimiaMaterialProperties_C* material_props) {
  if (material_props != NULL) {
    free(material_props->isotherm_kd);
    free(material_props->freundlich_n);
    free(material_props->langmuir_b);
  }
  material_props->isotherm_kd = NULL;
  material_props->freundlich_n = NULL;
  material_props->langmuir_b = NULL;
}  // end FreeAlquimiaMaterialProperties()

/*******************************************************************************
 **
 **  Meta Data
 **
 *******************************************************************************/

void AllocateAlquimiaMetaData(const struct AlquimiaSizes_C* sizes,
                              struct AlquimiaMetaData_C* meta_data) {
  int i;

  meta_data->primary_indices = NULL;
  meta_data->primary_names = NULL;
  //meta_data-> = NULL;

  if (sizes->num_primary > 0) {
    meta_data->primary_indices = (int*) calloc(sizes->num_primary, sizeof(int));
    if (NULL == meta_data->primary_indices) {
      // TODO(bja): error handling
    }
    meta_data->primary_names = (char**) calloc(sizes->num_primary, sizeof(char*));
    if (NULL == meta_data->primary_names) {
      // TODO(bja): error handling
    }
    for (i = 0; i < sizes->num_primary; ++i) {
      meta_data->primary_names[i] = (char*) calloc(ALQUIMIA_MAX_STRING_LENGTH, sizeof(char));
      if (NULL == meta_data->primary_names[i]) {
        // TODO(bja): error handling
      }
    }
  }
}  // end AllocateAlquimiaMetaData()

void FreeAlquimiaMetaData(const struct AlquimiaSizes_C* sizes,
                          struct AlquimiaMetaData_C* meta_data) {
  int i;
  if (meta_data != NULL) {
    free(meta_data->primary_indices);
    for (i = 0; i < sizes->num_primary; ++i) {
      free(meta_data->primary_names[i]);
    }
    free(meta_data->primary_names);
  }
  meta_data->primary_indices = NULL;
  meta_data->primary_names = NULL;

  // NOTE(bja): meta_data->internal_state is not our memory, but is allocated by
  // the engine, we rely on them to clean it up!

}  // end FreeAlquimiaMetaData()



/*******************************************************************************
 **
 **  Geochemical conditions/constraints
 **
 *******************************************************************************/

void AllocateAlquimiaGeochemicalConditionList(
    const int num_conditions,
    struct AlquimiaGeochemicalConditionList_C* condition_list) {
  // NOTE: we are only allocating pointers to N conditions here, not
  // the actual conditions themselves.
  fprintf(stdout, " AllocateAlquimiaGeochemicalConditionList() : %d\n",
          num_conditions);
  condition_list->num_conditions = num_conditions;

  condition_list->conditions = NULL;
  if (condition_list->num_conditions > 0) {
    condition_list->conditions = (struct AlquimiaGeochemicalCondition_C*)
        calloc(condition_list->num_conditions, 
               sizeof(struct AlquimiaGeochemicalCondition_C));
  }
}  // end AllocateAlquimiaGeochemicalConditionList()

void AllocateAlquimiaGeochemicalCondition(
    const char* name, const int num_constraints,
    struct AlquimiaGeochemicalCondition_C* condition) {
  // NOTE: we are only allocating pointers to N constraints here, not
  // the actual condstraints themselves.
  condition->num_constraints = num_constraints;

  condition->name = (char*) calloc(ALQUIMIA_MAX_STRING_LENGTH, sizeof(char));
  int max_copy_length = ALQUIMIA_MAX_STRING_LENGTH;
  if (strlen(name) < ALQUIMIA_MAX_STRING_LENGTH) {
    max_copy_length = strlen(name);
  }
  strncpy(condition->name, name, max_copy_length);

  condition->constraints = NULL;
  if (condition->num_constraints > 0) {
    condition->constraints = (struct AlquimiaGeochemicalConstraint_C*)
        calloc(condition->num_constraints, 
               sizeof(struct AlquimiaGeochemicalConstraint_C));
  }
}  // end AllocateAlquimiaGeochemicalCondition()

void AllocateAlquimiaGeochemicalConstraint(
    struct AlquimiaGeochemicalConstraint_C* constraint){
  constraint->primary_species =
      (char*) calloc(ALQUIMIA_MAX_STRING_LENGTH, sizeof(char));
  constraint->constraint_type =
      (char*) calloc(ALQUIMIA_MAX_STRING_LENGTH, sizeof(char));
  constraint->associated_species =
      (char*) calloc(ALQUIMIA_MAX_STRING_LENGTH, sizeof(char));
  constraint->value = 0.0;
}  // end AllocateAlquimiaGeochemicalConstraint()

void FreeAlquimiaGeochemicalConditionList(
    struct AlquimiaGeochemicalConditionList_C* condition_list) {
  int i;
  for (i = 0; i < condition_list->num_conditions; ++i) {
    FreeAlquimiaGeochemicalCondition(&(condition_list->conditions[i]));
  }
  free(condition_list->conditions);
  condition_list->conditions = NULL;
}  // end FreeAlquimiaGeochemicalConditionList()

void FreeAlquimiaGeochemicalCondition(
    struct AlquimiaGeochemicalCondition_C* condition) {
  int i;
  for (i = 0; i < condition->num_constraints; ++i) {
    FreeAlquimiaGeochemicalConstraint(&(condition->constraints[i]));
  }
  free(condition->constraints);
  condition->constraints = NULL;
}  // end FreeAlquimiaGeochemicalCondition()

void FreeAlquimiaGeochemicalConstraint(
    struct AlquimiaGeochemicalConstraint_C* constraint) {
  free(constraint->primary_species);
  constraint->primary_species = NULL;
  free(constraint->constraint_type);
  constraint->constraint_type = NULL;
  free(constraint->associated_species);
  constraint->associated_species = NULL;
}  // end FreeAlquimiaGeochemicalCondition()

