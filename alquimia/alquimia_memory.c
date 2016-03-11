/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
** Alquimia Copyright (c) 2013-2016, The Regents of the University of California, 
** through Lawrence Berkeley National Laboratory (subject to receipt of any 
** required approvals from the U.S. Dept. of Energy).  All rights reserved.
** 
** Alquimia is available under a BSD license. See LICENSE.txt for more
** information.
**
** If you have questions about your rights to use or distribute this software, 
** please contact Berkeley Lab's Technology Transfer and Intellectual Property 
** Management at TTD@lbl.gov referring to Alquimia (LBNL Ref. 2013-119).
** 
** NOTICE.  This software was developed under funding from the U.S. Department 
** of Energy.  As such, the U.S. Government has been granted for itself and 
** others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide 
** license in the Software to reproduce, prepare derivative works, and perform 
** publicly and display publicly.  Beginning five (5) years after the date 
** permission to assert copyright is obtained from the U.S. Department of Energy, 
** and subject to any subsequent five (5) year renewals, the U.S. Government is 
** granted for itself and others acting on its behalf a paid-up, nonexclusive, 
** irrevocable, worldwide license in the Software to reproduce, prepare derivative
** works, distribute copies to the public, perform publicly and display publicly, 
** and to permit others to do so.
** 
** Authors: Benjamin Andre <bandre@lbl.gov>
*/


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

#include "alquimia/alquimia_memory.h"
#include "alquimia/alquimia_interface.h"
#include "alquimia/alquimia_constants.h"
#include "alquimia/alquimia_containers.h"

// Returns the nearest power of 2 greater than or equal to n, or 0 if n == 0.
static inline int nearest_power_of_2(int n)
{
  if (n == 0) return 0;
  int twop = 1;
  while (twop < n)
    twop *= 2;
  return twop;
}

/*******************************************************************************
 **
 **  Alquimia Vectors
 **
 *******************************************************************************/
void AllocateAlquimiaVectorDouble(const int size, AlquimiaVectorDouble* vector) {
  if (size > 0) {
    vector->size = size;
    vector->capacity = nearest_power_of_2(size);
    vector->data = (double*) calloc((size_t)vector->capacity, sizeof(double));
    ALQUIMIA_ASSERT(NULL != vector->data);
  } else {
    vector->size = 0;
    vector->capacity = 0;
    vector->data = NULL;
  }
}  /* end AllocateAlquimiaVectorDouble() */

void FreeAlquimiaVectorDouble(AlquimiaVectorDouble* vector) {
  if (vector != NULL) {
    free(vector->data);
    vector->data = NULL;
    vector->size = 0;
    vector->capacity = 0;
  }
}  /* end FreeAlquimiaVectorDouble() */

void AllocateAlquimiaVectorInt(const int size, AlquimiaVectorInt* vector) {
  if (size > 0) {
    vector->size = size;
    vector->capacity = nearest_power_of_2(size);
    vector->data = (int*) calloc((size_t)vector->capacity, sizeof(int));
    ALQUIMIA_ASSERT(NULL != vector->data);
  } else {
    vector->size = 0;
    vector->capacity = 0;
    vector->data = NULL;
  }
}  /* end AllocateAlquimiaVectorInt() */
 
void FreeAlquimiaVectorInt(AlquimiaVectorInt* vector) {
  if (vector != NULL) {
    free(vector->data);
    vector->data = NULL;
    vector->size = 0;
    vector->capacity = 0;
  }
}  /* end FreeAlquimiaVectorInt() */

void AllocateAlquimiaVectorString(const int size, AlquimiaVectorString* vector) {
  int i;
  if (size > 0) {
    vector->size = size;
    vector->capacity = nearest_power_of_2(size);
    vector->data = (char**) calloc((size_t)vector->capacity, sizeof(char*));
    ALQUIMIA_ASSERT(NULL != vector->data);
    for (i = 0; i < vector->size; ++i) {
      vector->data[i] = (char*) calloc((size_t)kAlquimiaMaxStringLength, sizeof(char));
      ALQUIMIA_ASSERT(NULL != vector->data[i]);
    }
  } else {
    vector->size = 0;
    vector->capacity = 0;
    vector->data = NULL;
  }
}  /* end AllocateAlquimiaVectorString() */
 
void FreeAlquimiaVectorString(AlquimiaVectorString* vector) {
  int i;
  if (vector != NULL) {
    for (i = 0; i < vector->size; ++i) {
      free(vector->data[i]);
    }
    free(vector->data);
    vector->data = NULL;
    vector->size = 0;
    vector->capacity = 0;
  }
}  /* end FreeAlquimiaVectorString() */

/*******************************************************************************
 **
 **  State
 **
 *******************************************************************************/

void AllocateAlquimiaState(const AlquimiaSizes* const sizes,
                           AlquimiaState* state) {
  AllocateAlquimiaVectorDouble(sizes->num_primary, &(state->total_mobile));
  ALQUIMIA_ASSERT(state->total_mobile.data != NULL);

  AllocateAlquimiaVectorDouble(sizes->num_sorbed, &(state->total_immobile));

  AllocateAlquimiaVectorDouble(sizes->num_surface_sites,
                               &(state->surface_site_density));

  AllocateAlquimiaVectorDouble(sizes->num_ion_exchange_sites,
                               &(state->cation_exchange_capacity));

  AllocateAlquimiaVectorDouble(sizes->num_minerals,
                               &(state->mineral_volume_fraction));

  AllocateAlquimiaVectorDouble(sizes->num_minerals,
                               &(state->mineral_specific_surface_area));
}  /* end AllocateAlquimiaState() */

void FreeAlquimiaState(AlquimiaState* state) {
  if (state != NULL) {
    FreeAlquimiaVectorDouble(&(state->total_mobile));
    FreeAlquimiaVectorDouble(&(state->total_immobile));
    FreeAlquimiaVectorDouble(&(state->mineral_volume_fraction));
    FreeAlquimiaVectorDouble(&(state->mineral_specific_surface_area));
    FreeAlquimiaVectorDouble(&(state->cation_exchange_capacity));
    FreeAlquimiaVectorDouble(&(state->surface_site_density));
  }
}  /* end FreeAlquimiaState() */

/*******************************************************************************
 **
 **  Auxiliary Data
 **
 *******************************************************************************/

void AllocateAlquimiaAuxiliaryData(const AlquimiaSizes* const sizes,
                                   AlquimiaAuxiliaryData* aux_data) {
  AllocateAlquimiaVectorInt(sizes->num_aux_integers,
                            &(aux_data->aux_ints));

  AllocateAlquimiaVectorDouble(sizes->num_aux_doubles,
                               &(aux_data->aux_doubles));

}  /* end AllocateAlquimiaAuxiliaryData() */

void FreeAlquimiaAuxiliaryData(AlquimiaAuxiliaryData* aux_data) {
  if (aux_data != NULL) {
    FreeAlquimiaVectorInt(&(aux_data->aux_ints));
    FreeAlquimiaVectorDouble(&(aux_data->aux_doubles));
  }
}  /* end FreeAlquimiaAuxiliaryData() */

/*******************************************************************************
 **
 **  Properties
 **
 *******************************************************************************/

void AllocateAlquimiaProperties(const AlquimiaSizes* const sizes,
                                AlquimiaProperties* props) {
  AllocateAlquimiaVectorDouble(sizes->num_isotherm_species,
                               &(props->isotherm_kd));
  AllocateAlquimiaVectorDouble(sizes->num_isotherm_species,
                               &(props->freundlich_n));
  AllocateAlquimiaVectorDouble(sizes->num_isotherm_species,
                               &(props->langmuir_b));
  AllocateAlquimiaVectorDouble(sizes->num_minerals,
                               &(props->mineral_rate_cnst));
  AllocateAlquimiaVectorDouble(sizes->num_aqueous_kinetics,
                               &(props->aqueous_kinetic_rate_cnst));

}  /* end AllocateAlquimiaProperties() */

void FreeAlquimiaProperties(AlquimiaProperties* props) {
  if (props != NULL) {
    FreeAlquimiaVectorDouble(&(props->isotherm_kd));
    FreeAlquimiaVectorDouble(&(props->freundlich_n));
    FreeAlquimiaVectorDouble(&(props->langmuir_b));
    FreeAlquimiaVectorDouble(&(props->mineral_rate_cnst));
    FreeAlquimiaVectorDouble(&(props->aqueous_kinetic_rate_cnst));
  }
}  /* end FreeAlquimiaProperties() */

/*******************************************************************************
 **
 **  Problem Meta Data
 **
 *******************************************************************************/

void AllocateAlquimiaProblemMetaData(const AlquimiaSizes* const sizes,
                                     AlquimiaProblemMetaData* meta_data) {

  AllocateAlquimiaVectorString(sizes->num_primary, &(meta_data->primary_names));
  ALQUIMIA_ASSERT(meta_data->primary_names.data != NULL);

  AllocateAlquimiaVectorInt(sizes->num_primary, &(meta_data->positivity));
  memset(meta_data->positivity.data, 0, sizeof(int) * sizes->num_primary);

  AllocateAlquimiaVectorString(sizes->num_minerals,
                               &(meta_data->mineral_names));

  AllocateAlquimiaVectorString(sizes->num_surface_sites,
                               &(meta_data->surface_site_names));

  AllocateAlquimiaVectorString(sizes->num_ion_exchange_sites,
                               &(meta_data->ion_exchange_names));

  AllocateAlquimiaVectorString(sizes->num_isotherm_species,
                               &(meta_data->isotherm_species_names));

  AllocateAlquimiaVectorString(sizes->num_aqueous_kinetics,
                               &(meta_data->aqueous_kinetic_names));

}  /* end AllocateAlquimiaProblemMetaData() */

void FreeAlquimiaProblemMetaData(AlquimiaProblemMetaData* meta_data) {

  if (meta_data != NULL) {
    FreeAlquimiaVectorString(&(meta_data->primary_names));
    FreeAlquimiaVectorInt(&(meta_data->positivity));
    FreeAlquimiaVectorString(&(meta_data->mineral_names));
    FreeAlquimiaVectorString(&(meta_data->surface_site_names));
    FreeAlquimiaVectorString(&(meta_data->ion_exchange_names));
    FreeAlquimiaVectorString(&(meta_data->isotherm_species_names));
    FreeAlquimiaVectorString(&(meta_data->aqueous_kinetic_names));
  }
}  /* end FreeAlquimiaProblemMetaData() */

/*******************************************************************************
 **
 **  Auxiliary Output Data
 **
 *******************************************************************************/

void AllocateAlquimiaAuxiliaryOutputData(const AlquimiaSizes* const sizes,
                                         AlquimiaAuxiliaryOutputData* aux_output) {
  aux_output->pH = -999.9;
  AllocateAlquimiaVectorDouble(sizes->num_minerals,
                               &(aux_output->mineral_saturation_index));

  AllocateAlquimiaVectorDouble(sizes->num_aqueous_kinetics,
                               &(aux_output->aqueous_kinetic_rate));

  AllocateAlquimiaVectorDouble(sizes->num_minerals,
                               &(aux_output->mineral_reaction_rate));

  AllocateAlquimiaVectorDouble(sizes->num_primary,
                               &(aux_output->primary_free_ion_concentration));
  AllocateAlquimiaVectorDouble(sizes->num_primary,
                               &(aux_output->primary_activity_coeff));

  AllocateAlquimiaVectorDouble(sizes->num_aqueous_complexes,
                               &(aux_output->secondary_free_ion_concentration));
  AllocateAlquimiaVectorDouble(sizes->num_aqueous_complexes,
                               &(aux_output->secondary_activity_coeff));

}  /* end AllocateAlquimiaAuxiliaryOutputData() */

void FreeAlquimiaAuxiliaryOutputData(AlquimiaAuxiliaryOutputData* aux_output) {
  if (aux_output != NULL) {
    FreeAlquimiaVectorDouble(&(aux_output->aqueous_kinetic_rate));
    FreeAlquimiaVectorDouble(&(aux_output->mineral_saturation_index));
    FreeAlquimiaVectorDouble(&(aux_output->mineral_reaction_rate));
    FreeAlquimiaVectorDouble(&(aux_output->primary_free_ion_concentration));
    FreeAlquimiaVectorDouble(&(aux_output->primary_activity_coeff));
    FreeAlquimiaVectorDouble(&(aux_output->secondary_free_ion_concentration));
    FreeAlquimiaVectorDouble(&(aux_output->secondary_activity_coeff));
  }
}  /* end FreeAlquimiaAuxiliaryOutputData() */

/*******************************************************************************
 **
 **  Engine Status
 **
 *******************************************************************************/

void AllocateAlquimiaEngineStatus(AlquimiaEngineStatus* status) {

  status->message = (char*) calloc((size_t)kAlquimiaMaxStringLength, sizeof(char));
  if (NULL == status->message) {
    /* TODO(bja): error handling */
  }
}  /* end AllocateAlquimiaEngineStatus() */

void FreeAlquimiaEngineStatus(AlquimiaEngineStatus* status) {
  if (status != NULL) {
    free(status->message);
  }
  status->message = NULL;

}  /* end FreeAlquimiaEngineStatus() */



/*******************************************************************************
 **
 **  Geochemical conditions/constraints
 **
 *******************************************************************************/

void AllocateAlquimiaGeochemicalConditionVector(const int num_conditions,
                                                AlquimiaGeochemicalConditionVector* condition_list) {
  /* NOTE: we are only allocating pointers to N conditions here, not
     the actual conditions themselves. */
  fprintf(stdout, " AllocateAlquimiaGeochemicalConditionList() : %d\n",
          num_conditions);
  condition_list->size = num_conditions;
  condition_list->capacity = nearest_power_of_2(num_conditions);

  if (condition_list->size > 0) {
    condition_list->data = (AlquimiaGeochemicalCondition*)
        calloc((size_t)condition_list->capacity, 
               sizeof(AlquimiaGeochemicalCondition));
  }
}  /* end AllocateAlquimiaGeochemicalConditionVector() */

void AllocateAlquimiaGeochemicalCondition(const int size_name,
                                          const int num_aqueous_constraints, 
                                          const int num_mineral_constraints,
    AlquimiaGeochemicalCondition* condition) {
  /* NOTE: we are only allocating pointers to N constraints here, not
     the actual condstraints themselves. */
  if (condition != NULL) {
    /* size_name + 1 to include the null character! */
    condition->name = (char*) calloc((size_t)size_name+1, sizeof(char));
    AllocateAlquimiaAqueousConstraintVector(num_aqueous_constraints, &condition->aqueous_constraints);
    AllocateAlquimiaMineralConstraintVector(num_mineral_constraints, &condition->mineral_constraints);
  }
}  /* end AllocateAlquimiaGeochemicalCondition() */

void AllocateAlquimiaAqueousConstraint(AlquimiaAqueousConstraint* constraint) {
  constraint->primary_species_name =
      (char*) calloc((size_t)kAlquimiaMaxStringLength, sizeof(char));
  constraint->constraint_type =
      (char*) calloc((size_t)kAlquimiaMaxStringLength, sizeof(char));
  constraint->associated_species =
      (char*) calloc((size_t)kAlquimiaMaxStringLength, sizeof(char));
  constraint->value = 0.0;
}  /* end AllocateAlquimiaAqueousConstraint() */

void AllocateAlquimiaAqueousConstraintVector(int num_constraints,
                                             AlquimiaAqueousConstraintVector* constraint_list) {
  constraint_list->size = num_constraints;
  constraint_list->capacity = nearest_power_of_2(num_constraints);
  if (constraint_list->size > 0) {
    constraint_list->data = (AlquimiaAqueousConstraint*)
      calloc((size_t)constraint_list->capacity, 
             sizeof(AlquimiaAqueousConstraint));
  }
  else
    constraint_list->data = NULL;
}

void AllocateAlquimiaMineralConstraint(AlquimiaMineralConstraint* constraint) {
  constraint->mineral_name =
      (char*) calloc((size_t)kAlquimiaMaxStringLength, sizeof(char));
  constraint->volume_fraction = -1.0;
  constraint->specific_surface_area = -1.0;
}  /* end AllocateAlquimiaMineralConstraint() */

void AllocateAlquimiaMineralConstraintVector(int num_constraints,
                                             AlquimiaMineralConstraintVector* constraint_list) {
  constraint_list->size = num_constraints;
  constraint_list->capacity = nearest_power_of_2(num_constraints);
  if (constraint_list->size > 0) {
    constraint_list->data = (AlquimiaMineralConstraint*)
      calloc((size_t)constraint_list->capacity, 
             sizeof(AlquimiaMineralConstraint));
  }
  else
    constraint_list->data = NULL;
}

void FreeAlquimiaGeochemicalConditionVector(AlquimiaGeochemicalConditionVector* condition_list) {
  int i;
  if (condition_list != NULL) {
    for (i = 0; i < condition_list->size; ++i) {
      FreeAlquimiaGeochemicalCondition(&(condition_list->data[i]));
    }
    if (condition_list->data != NULL) {
      free(condition_list->data);
      condition_list->data = NULL;
    }
    condition_list->size = 0;
    condition_list->capacity = 0;
  }
}  /* end FreeAlquimiaGeochemicalConditionList() */

void FreeAlquimiaGeochemicalCondition(AlquimiaGeochemicalCondition* condition) {
  if (condition != NULL) {
    if (condition->name != NULL) {
      free(condition->name);
      condition->name = NULL;
    }
    FreeAlquimiaAqueousConstraintVector(&(condition->aqueous_constraints));
    FreeAlquimiaMineralConstraintVector(&(condition->mineral_constraints));
  }
}  /* end FreeAlquimiaGeochemicalCondition() */

void FreeAlquimiaAqueousConstraintVector(AlquimiaAqueousConstraintVector* vector) {
  int i;
  if (vector != NULL) {
    for (i = 0; i < vector->size; ++i) {
      FreeAlquimiaAqueousConstraint(&vector->data[i]);
    }
    if (vector->data != NULL) {
      free(vector->data);
      vector->data = NULL;
    }
    vector->size = 0;
    vector->capacity = 0;
  }
}  /* end FreeAlquimiaAqueousConstraintVector() */

void FreeAlquimiaAqueousConstraint(AlquimiaAqueousConstraint* constraint) {
  free(constraint->primary_species_name);
  constraint->primary_species_name = NULL;
  free(constraint->constraint_type);
  constraint->constraint_type = NULL;
  free(constraint->associated_species);
  constraint->associated_species = NULL;
}  /* end FreeAlquimiaAqueousConstraint() */

void FreeAlquimiaMineralConstraintVector(AlquimiaMineralConstraintVector* vector) {
  int i;
  if (vector != NULL) {
    for (i = 0; i < vector->size; ++i) {
      FreeAlquimiaMineralConstraint(&vector->data[i]);
    }
    free(vector->data);
    vector->data = NULL;
    vector->size = 0;
    vector->capacity = 0;
  }
}  /* end FreeAlquimiaMineralConstraintVector() */

void FreeAlquimiaMineralConstraint(AlquimiaMineralConstraint* constraint) {
  free(constraint->mineral_name);
  constraint->mineral_name = NULL;
}  /* end FreeAlquimiaMineralConstraint() */


/*******************************************************************************
 **
 **  Data convenience struct
 **
 *******************************************************************************/
void AllocateAlquimiaData(AlquimiaData* data) {
    AllocateAlquimiaState(&data->sizes, &data->state);
    AllocateAlquimiaProperties(&data->sizes, &data->properties);
    AllocateAlquimiaAuxiliaryData(&data->sizes, &data->aux_data);
    AllocateAlquimiaProblemMetaData(&data->sizes, &data->meta_data);
    AllocateAlquimiaAuxiliaryOutputData(&data->sizes, &data->aux_output);
}  /* end AllocateAlquimiaData() */


void FreeAlquimiaData(AlquimiaData* data) {
  FreeAlquimiaState(&data->state);
  FreeAlquimiaProperties(&data->properties);
  FreeAlquimiaAuxiliaryData(&data->aux_data);
  FreeAlquimiaProblemMetaData(&data->meta_data);
  FreeAlquimiaAuxiliaryOutputData(&data->aux_output);
}  /* end FreeAlquimiaData() */
