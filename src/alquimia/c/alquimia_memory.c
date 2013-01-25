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

/*******************************************************************************
 **
 **  Alquimia Vectors
 **
 *******************************************************************************/
void AllocateAlquimiaVectorDouble(const int size, struct AlquimiaVectorDouble* vector) {
  if (size > 0) {
    vector->size = size;
    vector->data = (double*) calloc(size, sizeof(double));
    assert(NULL != vector->data);
  } else {
    vector->size = 0;
    vector->data = NULL;
  }
}  // end AllocateAlquimiaVectorDouble()

void FreeAlquimiaVectorDouble(struct AlquimiaVectorDouble* vector) {
  if (vector != NULL) {
    free(vector->data);
    vector->data = NULL;
    vector->size = 0;
  }
}  // end FreeAlquimiaVectorDouble()

void AllocateAlquimiaVectorInt(const int size, struct AlquimiaVectorInt* vector) {
  if (size > 0) {
    vector->size = size;
    vector->data = (int*) calloc(size, sizeof(int));
    assert(NULL != vector->data);
  } else {
    vector->size = 0;
    vector->data = NULL;
  }
}  // end AllocateAlquimiaVectorInt()
 
void FreeAlquimiaVectorInt(struct AlquimiaVectorInt* vector) {
  if (vector != NULL) {
    free(vector->data);
    vector->data = NULL;
    vector->size = 0;
  }
}  // end FreeAlquimiaVectorInt()

void AllocateAlquimiaVectorString(const int size, struct AlquimiaVectorString* vector) {
  if (size > 0) {
    vector->size = size;
    vector->data = (char**) calloc(size, sizeof(char*));
    assert(NULL != vector->data);

    for (int i = 0; i < vector->size; ++i) {
      vector->data[i] = (char*) calloc(kAlquimiaMaxStringLength, sizeof(char));
      assert(NULL != vector->data[i]);
    }
  } else {
    vector->size = 0;
    vector->data = NULL;
  }
}  // end AllocateAlquimiaVectorString()
 
void FreeAlquimiaVectorString(struct AlquimiaVectorString* vector) {
  if (vector != NULL) {
    for (int i = 0; i < vector->size; ++i) {
      free(vector->data[i]);
    }
    free(vector->data);
    vector->data = NULL;
    vector->size = 0;
  }
}  // end FreeAlquimiaVectorString()

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
  AllocateAlquimiaVectorDouble(sizes->num_primary, &(state->total_primary));
  assert(state->total_primary.data != NULL);

  AllocateAlquimiaVectorDouble(sizes->num_primary, &(state->free_ion));
  assert(state->free_ion.data != NULL);

  AllocateAlquimiaVectorDouble(sizes->num_sorbed, &(state->total_sorbed));

  AllocateAlquimiaVectorDouble(sizes->num_surface_sites,
                               &(state->surface_site_density));

  AllocateAlquimiaVectorDouble(sizes->num_ion_exchange_sites,
                               &(state->cation_exchange_capacity));

  AllocateAlquimiaVectorDouble(sizes->num_kinetic_minerals,
                               &(state->mineral_volume_fraction));

  AllocateAlquimiaVectorDouble(sizes->num_kinetic_minerals,
                               &(state->mineral_specific_surface_area));
}  // end AllocateAlquimiaState()

void FreeAlquimiaState(struct AlquimiaState* state) {
  if (state != NULL) {
    FreeAlquimiaVectorDouble(&(state->total_primary));
    FreeAlquimiaVectorDouble(&(state->free_ion));
    FreeAlquimiaVectorDouble(&(state->mineral_volume_fraction));
    FreeAlquimiaVectorDouble(&(state->mineral_specific_surface_area));
    FreeAlquimiaVectorDouble(&(state->cation_exchange_capacity));
    FreeAlquimiaVectorDouble(&(state->surface_site_density));
  }
}  // end FreeAlquimiaState()

/*******************************************************************************
 **
 **  Auxiliary Data
 **
 *******************************************************************************/

void AllocateAlquimiaAuxiliaryData(const struct AlquimiaSizes* sizes,
                                   struct AlquimiaAuxiliaryData* aux_data) {
  AllocateAlquimiaVectorDouble(sizes->num_primary,
                               &(aux_data->primary_activity_coeff));

  AllocateAlquimiaVectorDouble(sizes->num_aqueous_complexes,
                               &(aux_data->secondary_activity_coeff));

  AllocateAlquimiaVectorDouble(sizes->num_ion_exchange_sites,
                               &(aux_data->ion_exchange_ref_cation_conc));

  AllocateAlquimiaVectorDouble(sizes->num_surface_sites,
                               &(aux_data->surface_complex_free_site_conc));

}  // end AllocateAlquimiaAuxiliaryData()

void FreeAlquimiaAuxiliaryData(struct AlquimiaAuxiliaryData* aux_data) {
  if (aux_data != NULL) {
    FreeAlquimiaVectorDouble(&(aux_data->primary_activity_coeff));
    FreeAlquimiaVectorDouble(&(aux_data->secondary_activity_coeff));
    FreeAlquimiaVectorDouble(&(aux_data->ion_exchange_ref_cation_conc));
    FreeAlquimiaVectorDouble(&(aux_data->surface_complex_free_site_conc));
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
  AllocateAlquimiaVectorDouble(sizes->num_primary,
                               &(material_props->isotherm_kd));
  AllocateAlquimiaVectorDouble(sizes->num_primary,
                               &(material_props->freundlich_n));
  AllocateAlquimiaVectorDouble(sizes->num_primary,
                               &(material_props->langmuir_b));

}  // end AllocateAlquimiaMaterialProperties()

void FreeAlquimiaMaterialProperties(
    struct AlquimiaMaterialProperties* material_props) {
  if (material_props != NULL) {
    FreeAlquimiaVectorDouble(&(material_props->isotherm_kd));
    FreeAlquimiaVectorDouble(&(material_props->freundlich_n));
    FreeAlquimiaVectorDouble(&(material_props->langmuir_b));
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
    AllocateAlquimiaVectorInt(sizes->num_primary, &(meta_data->primary_indices));
    assert(meta_data->primary_indices.data != NULL);

    AllocateAlquimiaVectorString(sizes->num_primary, &(meta_data->primary_names));
    assert(meta_data->primary_names.data != NULL);
  }

  if (sizes->num_kinetic_minerals > 0) {
    AllocateAlquimiaVectorInt(sizes->num_kinetic_minerals,
                              &(meta_data->mineral_indices));
    assert(meta_data->mineral_indices.data != NULL);

    AllocateAlquimiaVectorString(sizes->num_kinetic_minerals, &(meta_data->mineral_names));
    assert(meta_data->mineral_names.data != NULL);
  }
}  // end AllocateAlquimiaMetaData()

void FreeAlquimiaMetaData(struct AlquimiaMetaData* meta_data) {

  if (meta_data != NULL) {
    FreeAlquimiaVectorInt(&(meta_data->primary_indices));
    FreeAlquimiaVectorString(&(meta_data->primary_names));
    FreeAlquimiaVectorInt(&(meta_data->mineral_indices));
    FreeAlquimiaVectorString(&(meta_data->mineral_names));
  }
}  // end FreeAlquimiaMetaData()

/*******************************************************************************
 **
 **  Auxiliary Output Data
 **
 *******************************************************************************/

void AllocateAlquimiaAuxiliaryOutputData(
    const struct AlquimiaSizes* sizes,
    struct AlquimiaAuxiliaryOutputData* aux_output) {
  AllocateAlquimiaVectorDouble(sizes->num_kinetic_minerals,
                               &(aux_output->mineral_saturation_index));

  AllocateAlquimiaVectorDouble(sizes->num_kinetic_minerals,
                               &(aux_output->mineral_reaction_rate));


}  // end AllocateAlquimiaAuxiliaryOutputData()

void FreeAlquimiaAuxiliaryOutputData(
    struct AlquimiaAuxiliaryOutputData* aux_output) {
  if (aux_output != NULL) {
    FreeAlquimiaVectorDouble(&(aux_output->mineral_saturation_index));
    FreeAlquimiaVectorDouble(&(aux_output->mineral_reaction_rate));
  }
}  // end FreeAlquimiaAuxiliaryOutputData()

/*******************************************************************************
 **
 **  Engine Status
 **
 *******************************************************************************/

void AllocateAlquimiaEngineStatus(struct AlquimiaEngineStatus* status) {

  status->message = (char*) calloc(kAlquimiaMaxStringLength, sizeof(char));
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

  condition->name = (char*) calloc(kAlquimiaMaxStringLength, sizeof(char));
  int max_copy_length = kAlquimiaMaxStringLength;
  if (strlen(name) < kAlquimiaMaxStringLength) {
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
      (char*) calloc(kAlquimiaMaxStringLength, sizeof(char));
  constraint->constraint_type =
      (char*) calloc(kAlquimiaMaxStringLength, sizeof(char));
  constraint->associated_species =
      (char*) calloc(kAlquimiaMaxStringLength, sizeof(char));
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
    AllocateAlquimiaAuxiliaryOutputData(&data->sizes, &data->aux_output);
}  // end AllocateAlquimiaData()


void FreeAlquimiaData(struct AlquimiaData* data) {
  FreeAlquimiaState(&data->state);
  FreeAlquimiaMaterialProperties(&data->material_properties);
  FreeAlquimiaAuxiliaryData(&data->aux_data);
  FreeAlquimiaMetaData(&data->meta_data);
  FreeAlquimiaAuxiliaryOutputData(&data->aux_output);
}  // end FreeAlquimiaData()
