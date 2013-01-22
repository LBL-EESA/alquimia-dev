/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*******************************************************************************
 **
 **  Alquimia C memory utilities to handle memory management
 **
 **  Notes:
 **
 **  - calloc/malloc always return NULL pointers if they fail, so
 **    there is no need to pre-assign NULL for the pointers we are
 **    allocating here. For zero size or zero members, the returned
 **    pointer should be NULL or something that can be freed....
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
#include <assert.h>

#include "alquimia_interface.h"
#include "alquimia_constants.h"
#include "alquimia_containers.h"

void AllocateDoubleArray(const int size, int* array_size, double** array) {
  if (size > 0) {
    *array_size = size;
    *array = (double*) calloc(size, sizeof(double));
    assert(NULL != *array);
  } else {
    *array_size = 0;
    *array = NULL;
  }
}  // end AllocateDoubleArray()

void AllocateIntArray(const int size, int* array_size, int** array) {
  if (size > 0) {
    *array_size = size;
    *array = (int*) calloc(size, sizeof(int));
    assert(NULL != *array);
  } else {
    *array_size = 0;
    *array = NULL;
  }
}  // end AllocateIntArray()

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
  AllocateDoubleArray(sizes->num_primary, 
                      &(state->size_total_primary), &(state->total_primary));
  assert(state->total_primary != NULL);

  AllocateDoubleArray(sizes->num_primary, 
                      &(state->size_free_ion), &(state->free_ion));
  assert(state->free_ion != NULL);

  AllocateDoubleArray(sizes->num_sorbed,
                      &(state->size_total_sorbed), &(state->total_sorbed));

  AllocateDoubleArray(sizes->num_surface_sites,
                      &(state->size_surface_site_density),
                      &(state->surface_site_density));

  AllocateDoubleArray(sizes->num_ion_exchange_sites,
                      &(state->size_cation_exchange_capacity),
                      &(state->cation_exchange_capacity));

  AllocateDoubleArray(sizes->num_kinetic_minerals,
                      &(state->size_mineral_volume_fraction),
                      &(state->mineral_volume_fraction));

  AllocateDoubleArray(sizes->num_kinetic_minerals,
                      &(state->size_mineral_specific_surface_area),
                      &(state->mineral_specific_surface_area));
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

    free(state->cation_exchange_capacity);
    state->cation_exchange_capacity = NULL;

    free(state->surface_site_density);
    state->surface_site_density = NULL;
  }
}  // end FreeAlquimiaState()

/*******************************************************************************
 **
 **  Auxiliary Data
 **
 *******************************************************************************/

void AllocateAlquimiaAuxiliaryData(const struct AlquimiaSizes* sizes,
                                   struct AlquimiaAuxiliaryData* aux_data) {
  AllocateDoubleArray(sizes->num_primary,
                      &(aux_data->size_primary_activity_coeff),
                      &(aux_data->primary_activity_coeff));

  AllocateDoubleArray(sizes->num_aqueous_complexes,
                      &(aux_data->size_secondary_activity_coeff),
                      &(aux_data->secondary_activity_coeff));

  AllocateDoubleArray(sizes->num_ion_exchange_sites,
                      &(aux_data->size_ion_exchange_ref_cation_conc),
                      &(aux_data->ion_exchange_ref_cation_conc));

  AllocateDoubleArray(sizes->num_surface_sites,
                      &(aux_data->size_surface_complex_free_site_conc),
                      &(aux_data->surface_complex_free_site_conc));

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
  AllocateDoubleArray(sizes->num_primary,
                      &(material_props->size_isotherm_kd),
                      &(material_props->isotherm_kd));
  AllocateDoubleArray(sizes->num_primary,
                      &(material_props->size_freundlich_n),
                      &(material_props->freundlich_n));
  AllocateDoubleArray(sizes->num_primary,
                      &(material_props->size_langmuir_b),
                      &(material_props->langmuir_b));

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

  if (sizes->num_primary > 0) {
    AllocateIntArray(sizes->num_primary,
                     &(meta_data->size_primary), &(meta_data->primary_indices));
    assert(meta_data->primary_indices != NULL);

    meta_data->primary_names = (char**) calloc(meta_data->size_primary, sizeof(char*));
    assert(NULL != meta_data->primary_names);

    for (int i = 0; i < meta_data->size_primary; ++i) {
      meta_data->primary_names[i] = (char*) calloc(ALQUIMIA_MAX_STRING_LENGTH,
                                                   sizeof(char));
      assert(NULL != meta_data->primary_names[i]);
    }
  }
}  // end AllocateAlquimiaMetaData()

void FreeAlquimiaMetaData(struct AlquimiaMetaData* meta_data) {

  if (meta_data != NULL) {
    free(meta_data->primary_indices);
    meta_data->primary_indices = NULL;
    for (int i = 0; i < meta_data->size_primary; ++i) {
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


/*******************************************************************************
 **
 **  Data convenience struct
 **
 *******************************************************************************/
void AllocateAlquimiaData(struct AlquimiaData* data) {
    AllocateAlquimiaState(&data->sizes, &data->state);
    AllocateAlquimiaMaterialProperties(&data->sizes, &data->material_properties);
    AllocateAlquimiaAuxiliaryData(&data->sizes, &data->aux_data);
    AllocateAlquimiaMetaData(&data->sizes, &data->meta_data);

}  // end AllocateAlquimiaData()


void FreeAlquimiaData(struct AlquimiaData* data) {
  FreeAlquimiaState(&data->state);
  FreeAlquimiaMaterialProperties(&data->material_properties);
  FreeAlquimiaAuxiliaryData(&data->aux_data);
  FreeAlquimiaMetaData(&data->meta_data);
}  // end FreeAlquimiaData()
