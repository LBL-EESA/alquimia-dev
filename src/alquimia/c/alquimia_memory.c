/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*******************************************************************************
 **
 **  Alquimia C memory utilities to handle memory management
 **
 **  Notes:
 **
 **  - calloc/malloc always return NULL pointers if they fail, so
 **    there is no need to pre-assign NULL for the pointers we are
 **    allocating here. For pointers being assigned elsewhere, we still
 **    assign NULL just to be safe. In some places, the malloc calls
 **    are inside if blocks, so we pre-initialize to NULL to keep life
 **    simple.
 **
 ** - free just releases the memory, it does not change the value
 **   of the pointer. After free, the pointer is no longer valid, so
 **   we set it to NULL.
 **
 *******************************************************************************/

#include "alquimia_memory.h"

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "alquimia_constants.h"
#include "alquimia_interface.h"
#include "alquimia_containers.h"

/*******************************************************************************
 **
 **  Alquimia Interface
 **
 *******************************************************************************/

void AllocateAlquimiaInterface(struct AlquimiaInterface* interface) {
  interface->Setup = NULL;
  interface->Shutdown = NULL;
  interface->ProcessCondition = NULL;
  interface->ReactionStepOperatorSplit = NULL;
  interface->GetAuxiliaryOutput = NULL;
  interface->GetEngineMetaData = NULL;
  interface->GetPrimaryNameFromIndex = NULL;
  interface->engine_state = (void*) calloc(1, sizeof(void*));
}  // end AllocateAlquimiaInterface()

void FreeAlquimiaInterface(struct AlquimiaInterface* interface) {
  if (interface != NULL) {
    free(interface->engine_state);
    interface->engine_state = NULL;
  }
}  // end FreeAlquimiaInterface()

/*******************************************************************************
 **
 **  State
 **
 *******************************************************************************/

void AllocateAlquimiaState(const struct AlquimiaSizes* sizes,
                           struct AlquimiaState* state) {
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

void FreeAlquimiaState(struct AlquimiaState* state) {
  if (state != NULL) {
    free(state->total_primary);
    state->total_primary = NULL;

    free(state->free_ion);
    state->free_ion = NULL;

    free(state->mineral_volume_fraction);
    state->mineral_volume_fraction = NULL;

    free(state->mineral_specific_surface_area);
    state->mineral_specific_surface_area = NULL;
  }
}  // end FreeAlquimiaState()

/*******************************************************************************
 **
 **  Auxiliary Data
 **
 *******************************************************************************/

void AllocateAlquimiaAuxiliaryData(const struct AlquimiaSizes* sizes,
                                   struct AlquimiaAuxiliaryData* aux_data) {
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

void FreeAlquimiaAuxiliaryData(struct AlquimiaAuxiliaryData* aux_data) {
  if (aux_data != NULL) {
    free(aux_data->primary_activity_coeff);
    aux_data->primary_activity_coeff = NULL;

    free(aux_data->secondary_activity_coeff);
    aux_data->secondary_activity_coeff = NULL;

    free(aux_data->ion_exchange_ref_cation_conc);
    aux_data->ion_exchange_ref_cation_conc = NULL;

    free(aux_data->surface_complex_free_site_conc);
    aux_data->surface_complex_free_site_conc = NULL;
  }
}  // end FreeAlquimiaAuxiliaryData()

/*******************************************************************************
 **
 **  Material Properties
 **
 *******************************************************************************/

void AllocateAlquimiaMaterialProperties(
    const struct AlquimiaSizes* sizes,
    struct AlquimiaMaterialProperties* material_props) {
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
    struct AlquimiaMaterialProperties* material_props) {
  if (material_props != NULL) {
    free(material_props->isotherm_kd);
    material_props->isotherm_kd = NULL;

    free(material_props->freundlich_n);
    material_props->freundlich_n = NULL;

    free(material_props->langmuir_b);
    material_props->langmuir_b = NULL;
  }
}  // end FreeAlquimiaMaterialProperties()

/*******************************************************************************
 **
 **  Meta Data
 **
 *******************************************************************************/

void AllocateAlquimiaMetaData(const struct AlquimiaSizes* sizes,
                              struct AlquimiaMetaData* meta_data) {
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

void FreeAlquimiaMetaData(const struct AlquimiaSizes* sizes,
                          struct AlquimiaMetaData* meta_data) {
  int i;
  if (meta_data != NULL) {
    free(meta_data->primary_indices);
    meta_data->primary_indices = NULL;
    for (i = 0; i < sizes->num_primary; ++i) {
      free(meta_data->primary_names[i]);
    }
    free(meta_data->primary_names);
    meta_data->primary_names = NULL;
  }

}  // end FreeAlquimiaMetaData()

/*******************************************************************************
 **
 **  Engine Status
 **
 *******************************************************************************/

void AllocateAlquimiaEngineStatus(struct AlquimiaEngineStatus* status) {

  status->message = (char*) calloc(ALQUIMIA_MAX_STRING_LENGTH, sizeof(char));
  if (NULL == status->message) {
    // TODO(bja): error handling
  }
}  // end AllocateAlquimiaEngineStatus()

void FreeAlquimiaEngineStatus(struct AlquimiaEngineStatus* status) {
  if (status != NULL) {
    free(status->message);
  }
  status->message = NULL;

}  // end FreeAlquimiaEngineStatus()



/*******************************************************************************
 **
 **  Geochemical conditions/constraints
 **
 *******************************************************************************/

void AllocateAlquimiaGeochemicalConditionList(
    const int num_conditions,
    struct AlquimiaGeochemicalConditionList* condition_list) {
  // NOTE: we are only allocating pointers to N conditions here, not
  // the actual conditions themselves.
  fprintf(stdout, " AllocateAlquimiaGeochemicalConditionList() : %d\n",
          num_conditions);
  condition_list->num_conditions = num_conditions;

  condition_list->conditions = NULL;
  if (condition_list->num_conditions > 0) {
    condition_list->conditions = (struct AlquimiaGeochemicalCondition*)
        calloc(condition_list->num_conditions, 
               sizeof(struct AlquimiaGeochemicalCondition));
  }
}  // end AllocateAlquimiaGeochemicalConditionList()

void AllocateAlquimiaGeochemicalCondition(
    const char* name, const int num_constraints,
    struct AlquimiaGeochemicalCondition* condition) {
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
    condition->constraints = (struct AlquimiaGeochemicalConstraint*)
        calloc(condition->num_constraints, 
               sizeof(struct AlquimiaGeochemicalConstraint));
  }
}  // end AllocateAlquimiaGeochemicalCondition()

void AllocateAlquimiaGeochemicalConstraint(
    struct AlquimiaGeochemicalConstraint* constraint){
  constraint->primary_species =
      (char*) calloc(ALQUIMIA_MAX_STRING_LENGTH, sizeof(char));
  constraint->constraint_type =
      (char*) calloc(ALQUIMIA_MAX_STRING_LENGTH, sizeof(char));
  constraint->associated_species =
      (char*) calloc(ALQUIMIA_MAX_STRING_LENGTH, sizeof(char));
  constraint->value = 0.0;
}  // end AllocateAlquimiaGeochemicalConstraint()

void FreeAlquimiaGeochemicalConditionList(
    struct AlquimiaGeochemicalConditionList* condition_list) {

  for (int i = 0; i < condition_list->num_conditions; ++i) {
    FreeAlquimiaGeochemicalCondition(&(condition_list->conditions[i]));
  }
  free(condition_list->conditions);
  condition_list->conditions = NULL;
}  // end FreeAlquimiaGeochemicalConditionList()

void FreeAlquimiaGeochemicalCondition(
    struct AlquimiaGeochemicalCondition* condition) {
  if (condition != NULL) {
    free(condition->name);
    condition->name = NULL;
    for (int i = 0; i < condition->num_constraints; ++i) {
      FreeAlquimiaGeochemicalConstraint(&(condition->constraints[i]));
    }
    free(condition->constraints);
    condition->constraints = NULL;
  }
}  // end FreeAlquimiaGeochemicalCondition()

void FreeAlquimiaGeochemicalConstraint(
    struct AlquimiaGeochemicalConstraint* constraint) {
  free(constraint->primary_species);
  constraint->primary_species = NULL;
  free(constraint->constraint_type);
  constraint->constraint_type = NULL;
  free(constraint->associated_species);
  constraint->associated_species = NULL;
}  // end FreeAlquimiaGeochemicalCondition()

