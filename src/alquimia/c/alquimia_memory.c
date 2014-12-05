/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
** Alquimia Copyright (c) 2013, The Regents of the University of California, 
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
    vector->data = (double*) calloc((unsigned int)size, sizeof(double));
    assert(NULL != vector->data);
  } else {
    vector->size = 0;
    vector->data = NULL;
  }
}  /* end AllocateAlquimiaVectorDouble() */

void FreeAlquimiaVectorDouble(struct AlquimiaVectorDouble* vector) {
  if (vector != NULL) {
    free(vector->data);
    vector->data = NULL;
    vector->size = 0;
  }
}  /* end FreeAlquimiaVectorDouble() */

void AllocateAlquimiaVectorInt(const int size, struct AlquimiaVectorInt* vector) {
  if (size > 0) {
    vector->size = size;
    vector->data = (int*) calloc((unsigned int)size, sizeof(int));
    assert(NULL != vector->data);
  } else {
    vector->size = 0;
    vector->data = NULL;
  }
}  /* end AllocateAlquimiaVectorInt() */
 
void FreeAlquimiaVectorInt(struct AlquimiaVectorInt* vector) {
  if (vector != NULL) {
    free(vector->data);
    vector->data = NULL;
    vector->size = 0;
  }
}  /* end FreeAlquimiaVectorInt() */

void AllocateAlquimiaVectorString(const int size, struct AlquimiaVectorString* vector) {
  int i;
  if (size > 0) {
    vector->size = size;
    vector->data = (char**) calloc((unsigned int)size, sizeof(char*));
    assert(NULL != vector->data);
    for (i = 0; i < vector->size; ++i) {
      vector->data[i] = (char*) calloc((unsigned int)kAlquimiaMaxStringLength, sizeof(char));
      assert(NULL != vector->data[i]);
    }
  } else {
    vector->size = 0;
    vector->data = NULL;
  }
}  /* end AllocateAlquimiaVectorString() */
 
void FreeAlquimiaVectorString(struct AlquimiaVectorString* vector) {
  int i;
  if (vector != NULL) {
    for (i = 0; i < vector->size; ++i) {
      free(vector->data[i]);
    }
    free(vector->data);
    vector->data = NULL;
    vector->size = 0;
  }
}  /* end FreeAlquimiaVectorString() */

/*******************************************************************************
 **
 **  State
 **
 *******************************************************************************/

void AllocateAlquimiaState(const struct AlquimiaSizes* const sizes,
                           struct AlquimiaState* state) {
  AllocateAlquimiaVectorDouble(sizes->num_primary, &(state->total_mobile));
  assert(state->total_mobile.data != NULL);

  AllocateAlquimiaVectorDouble(sizes->num_sorbed, &(state->total_immobile));

  AllocateAlquimiaVectorDouble(sizes->num_surface_sites,
                               &(state->surface_site_density));

  AllocateAlquimiaVectorDouble(sizes->num_ion_exchange_sites,
                               &(state->cation_exchange_capacity));

  AllocateAlquimiaVectorDouble(sizes->num_kinetic_minerals,
                               &(state->mineral_volume_fraction));

  AllocateAlquimiaVectorDouble(sizes->num_kinetic_minerals,
                               &(state->mineral_specific_surface_area));
  
  AllocateAlquimiaVectorDouble(sizes->num_total_gases,
                               &(state->total_gas));
  
  AllocateAlquimiaVectorDouble(sizes->num_gas_species,
                               &(state->gas_concentration));
}  /* end AllocateAlquimiaState() */

void FreeAlquimiaState(struct AlquimiaState* state) {
  if (state != NULL) {
    FreeAlquimiaVectorDouble(&(state->total_mobile));
    FreeAlquimiaVectorDouble(&(state->total_immobile));
    FreeAlquimiaVectorDouble(&(state->mineral_volume_fraction));
    FreeAlquimiaVectorDouble(&(state->mineral_specific_surface_area));
    FreeAlquimiaVectorDouble(&(state->cation_exchange_capacity));
    FreeAlquimiaVectorDouble(&(state->surface_site_density));
    FreeAlquimiaVectorDouble(&(state->total_gas));
    FreeAlquimiaVectorDouble(&(state->gas_concentration));
  }
}  /* end FreeAlquimiaState() */

/*******************************************************************************
 **
 **  Auxiliary Data
 **
 *******************************************************************************/

void AllocateAlquimiaAuxiliaryData(const struct AlquimiaSizes* const sizes,
                                   struct AlquimiaAuxiliaryData* aux_data) {
  AllocateAlquimiaVectorInt(sizes->num_aux_integers,
                            &(aux_data->aux_ints));

  AllocateAlquimiaVectorDouble(sizes->num_aux_doubles,
                               &(aux_data->aux_doubles));

}  /* end AllocateAlquimiaAuxiliaryData() */

void FreeAlquimiaAuxiliaryData(struct AlquimiaAuxiliaryData* aux_data) {
  if (aux_data != NULL) {
    FreeAlquimiaVectorInt(&(aux_data->aux_ints));
    FreeAlquimiaVectorDouble(&(aux_data->aux_doubles));
  }
}  /* end FreeAlquimiaAuxiliaryData() */

/*******************************************************************************
 **
 **  Material Properties
 **
 *******************************************************************************/

void AllocateAlquimiaMaterialProperties(
    const struct AlquimiaSizes* const sizes,
    struct AlquimiaMaterialProperties* material_props) {
  AllocateAlquimiaVectorDouble(sizes->num_isotherm_species,
                               &(material_props->isotherm_kd));
  AllocateAlquimiaVectorDouble(sizes->num_isotherm_species,
                               &(material_props->freundlich_n));
  AllocateAlquimiaVectorDouble(sizes->num_isotherm_species,
                               &(material_props->langmuir_b));

}  /* end AllocateAlquimiaMaterialProperties() */

void FreeAlquimiaMaterialProperties(
    struct AlquimiaMaterialProperties* material_props) {
  if (material_props != NULL) {
    FreeAlquimiaVectorDouble(&(material_props->isotherm_kd));
    FreeAlquimiaVectorDouble(&(material_props->freundlich_n));
    FreeAlquimiaVectorDouble(&(material_props->langmuir_b));
  }
}  /* end FreeAlquimiaMaterialProperties() */

/*******************************************************************************
 **
 **  Problem Meta Data
 **
 *******************************************************************************/

void AllocateAlquimiaProblemMetaData(const struct AlquimiaSizes* const sizes,
                                     struct AlquimiaProblemMetaData* meta_data) {

  AllocateAlquimiaVectorString(sizes->num_primary, &(meta_data->primary_names));
  assert(meta_data->primary_names.data != NULL);

  AllocateAlquimiaVectorString(sizes->num_kinetic_minerals,
                               &(meta_data->mineral_names));

  AllocateAlquimiaVectorString(sizes->num_surface_sites,
                               &(meta_data->surface_site_names));

  AllocateAlquimiaVectorString(sizes->num_ion_exchange_sites,
                               &(meta_data->ion_exchange_names));

  AllocateAlquimiaVectorString(sizes->num_isotherm_species,
                               &(meta_data->isotherm_species_names));

  AllocateAlquimiaVectorString(sizes->num_gas_species,
                               &(meta_data->gas_names));

}  /* end AllocateAlquimiaProblemMetaData() */

void FreeAlquimiaProblemMetaData(struct AlquimiaProblemMetaData* meta_data) {

  if (meta_data != NULL) {
    FreeAlquimiaVectorString(&(meta_data->primary_names));
    FreeAlquimiaVectorString(&(meta_data->mineral_names));
    FreeAlquimiaVectorString(&(meta_data->surface_site_names));
    FreeAlquimiaVectorString(&(meta_data->ion_exchange_names));
    FreeAlquimiaVectorString(&(meta_data->isotherm_species_names));
	FreeAlquimiaVectorString(&(meta_data->gas_names));
  }
}  /* end FreeAlquimiaProblemMetaData() */

/*******************************************************************************
 **
 **  Auxiliary Output Data
 **
 *******************************************************************************/

void AllocateAlquimiaAuxiliaryOutputData(
    const struct AlquimiaSizes* const sizes,
    struct AlquimiaAuxiliaryOutputData* aux_output) {
  aux_output->pH = -999.9;
  AllocateAlquimiaVectorDouble(sizes->num_kinetic_minerals,
                               &(aux_output->mineral_saturation_index));

  AllocateAlquimiaVectorDouble(sizes->num_kinetic_minerals,
                               &(aux_output->mineral_reaction_rate));

  AllocateAlquimiaVectorDouble(sizes->num_primary,
                               &(aux_output->primary_free_ion_concentration));
  AllocateAlquimiaVectorDouble(sizes->num_primary,
                               &(aux_output->primary_activity_coeff));

  AllocateAlquimiaVectorDouble(sizes->num_aqueous_complexes,
                               &(aux_output->secondary_free_ion_concentration));
  AllocateAlquimiaVectorDouble(sizes->num_aqueous_complexes,
                               &(aux_output->secondary_activity_coeff));
  
  AllocateAlquimiaVectorDouble(sizes->num_gas_species,
                               &(aux_output->gas_partial_pressure));

}  /* end AllocateAlquimiaAuxiliaryOutputData() */

void FreeAlquimiaAuxiliaryOutputData(
    struct AlquimiaAuxiliaryOutputData* aux_output) {
  if (aux_output != NULL) {
    FreeAlquimiaVectorDouble(&(aux_output->mineral_saturation_index));
    FreeAlquimiaVectorDouble(&(aux_output->mineral_reaction_rate));
    FreeAlquimiaVectorDouble(&(aux_output->primary_free_ion_concentration));
    FreeAlquimiaVectorDouble(&(aux_output->primary_activity_coeff));
    FreeAlquimiaVectorDouble(&(aux_output->secondary_free_ion_concentration));
    FreeAlquimiaVectorDouble(&(aux_output->secondary_activity_coeff));
	FreeAlquimiaVectorDouble(&(aux_output->gas_partial_pressure));
  }
}  /* end FreeAlquimiaAuxiliaryOutputData() */

/*******************************************************************************
 **
 **  Engine Status
 **
 *******************************************************************************/

void AllocateAlquimiaEngineStatus(struct AlquimiaEngineStatus* status) {

  status->message = (char*) calloc((unsigned int)kAlquimiaMaxStringLength, sizeof(char));
  if (NULL == status->message) {
    /* TODO(bja): error handling */
  }
}  /* end AllocateAlquimiaEngineStatus() */

void FreeAlquimiaEngineStatus(struct AlquimiaEngineStatus* status) {
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

void AllocateAlquimiaGeochemicalConditionVector(
    const int num_conditions,
    struct AlquimiaGeochemicalConditionVector* condition_list) {
  /* NOTE: we are only allocating pointers to N conditions here, not
     the actual conditions themselves. */
  fprintf(stdout, " AllocateAlquimiaGeochemicalConditionList() : %d\n",
          num_conditions);
  condition_list->size = num_conditions;

  if (condition_list->size > 0) {
    condition_list->data = (struct AlquimiaGeochemicalCondition*)
        calloc((unsigned int)condition_list->size, 
               sizeof(struct AlquimiaGeochemicalCondition));
  }
}  /* end AllocateAlquimiaGeochemicalConditionVector() */

void AllocateAlquimiaGeochemicalCondition(
    const int size_name,
    const int num_aqueous_constraints, const int num_mineral_constraints,
    struct AlquimiaGeochemicalCondition* condition) {
  /* NOTE: we are only allocating pointers to N constraints here, not
     the actual condstraints themselves. */
  if (condition != NULL) {
    /* size_name + 1 to include the null character! */
    condition->name = (char*) calloc((unsigned int)size_name+1, sizeof(char));

    condition->aqueous_constraints.size = num_aqueous_constraints;
    if (condition->aqueous_constraints.size > 0) {
      condition->aqueous_constraints.data = (struct AlquimiaAqueousConstraint*)
          calloc((unsigned int)condition->aqueous_constraints.size, 
                 sizeof(struct AlquimiaAqueousConstraint));
    }

    condition->mineral_constraints.size = num_mineral_constraints;
    if (condition->mineral_constraints.size > 0) {
      condition->mineral_constraints.data = (struct AlquimiaMineralConstraint*)
          calloc((unsigned int)condition->mineral_constraints.size, 
                 sizeof(struct AlquimiaMineralConstraint));
    }

  }
}  /* end AllocateAlquimiaGeochemicalCondition() */

void AllocateAlquimiaAqueousConstraint(
    struct AlquimiaAqueousConstraint* constraint) {
  constraint->primary_species_name =
      (char*) calloc((unsigned int)kAlquimiaMaxStringLength, sizeof(char));
  constraint->constraint_type =
      (char*) calloc((unsigned int)kAlquimiaMaxStringLength, sizeof(char));
  constraint->associated_species =
      (char*) calloc((unsigned int)kAlquimiaMaxStringLength, sizeof(char));
  constraint->value = 0.0;
}  /* end AllocateAlquimiaAqueousConstraint() */

void AllocateAlquimiaMineralConstraint(
    struct AlquimiaMineralConstraint* constraint) {
  constraint->mineral_name =
      (char*) calloc((unsigned int)kAlquimiaMaxStringLength, sizeof(char));
  constraint->volume_fraction = -1.0;
  constraint->specific_surface_area = -1.0;
}  /* end AllocateAlquimiaMineralConstraint() */

void FreeAlquimiaGeochemicalConditionVector(
    struct AlquimiaGeochemicalConditionVector* condition_list) {
  int i;
  if (condition_list != NULL) {
    for (i = 0; i < condition_list->size; ++i) {
      FreeAlquimiaGeochemicalCondition(&(condition_list->data[i]));
    }
    free(condition_list->data);
    condition_list->data = NULL;
    condition_list->size = 0;
  }
}  /* end FreeAlquimiaGeochemicalConditionList() */

void FreeAlquimiaGeochemicalCondition(
    struct AlquimiaGeochemicalCondition* condition) {
  if (condition != NULL) {
    free(condition->name);
    condition->name = NULL;
    FreeAlquimiaAqueousConstraintVector(&(condition->aqueous_constraints));
    FreeAlquimiaMineralConstraintVector(&(condition->mineral_constraints));
  }
}  /* end FreeAlquimiaGeochemicalCondition() */

void FreeAlquimiaAqueousConstraintVector(
    struct AlquimiaAqueousConstraintVector* vector) {
  int i;
  if (vector != NULL) {
    for (i = 0; i < vector->size; ++i) {
      FreeAlquimiaAqueousConstraint(&vector->data[i]);
    }
    free(vector->data);
    vector->data = NULL;
    vector->size = 0;
  }
}  /* end FreeAlquimiaAqueousConstraintVector() */

void FreeAlquimiaAqueousConstraint(
    struct AlquimiaAqueousConstraint* constraint) {
  free(constraint->primary_species_name);
  constraint->primary_species_name = NULL;
  free(constraint->constraint_type);
  constraint->constraint_type = NULL;
  free(constraint->associated_species);
  constraint->associated_species = NULL;
}  /* end FreeAlquimiaAqueousConstraint() */

void FreeAlquimiaMineralConstraintVector(
    struct AlquimiaMineralConstraintVector* vector) {
  int i;
  if (vector != NULL) {
    for (i = 0; i < vector->size; ++i) {
      FreeAlquimiaMineralConstraint(&vector->data[i]);
    }
    free(vector->data);
    vector->data = NULL;
    vector->size = 0;
  }
}  /* end FreeAlquimiaMineralConstraintVector() */

void FreeAlquimiaMineralConstraint(
    struct AlquimiaMineralConstraint* constraint) {
  free(constraint->mineral_name);
  constraint->mineral_name = NULL;
}  /* end FreeAlquimiaMineralConstraint() */


/*******************************************************************************
 **
 **  Data convenience struct
 **
 *******************************************************************************/
void AllocateAlquimiaData(struct AlquimiaData* data) {
    AllocateAlquimiaState(&data->sizes, &data->state);
    AllocateAlquimiaMaterialProperties(&data->sizes, &data->material_properties);
    AllocateAlquimiaAuxiliaryData(&data->sizes, &data->aux_data);
    AllocateAlquimiaProblemMetaData(&data->sizes, &data->meta_data);
    AllocateAlquimiaAuxiliaryOutputData(&data->sizes, &data->aux_output);
}  /* end AllocateAlquimiaData() */


void FreeAlquimiaData(struct AlquimiaData* data) {
  FreeAlquimiaState(&data->state);
  FreeAlquimiaMaterialProperties(&data->material_properties);
  FreeAlquimiaAuxiliaryData(&data->aux_data);
  FreeAlquimiaProblemMetaData(&data->meta_data);
  FreeAlquimiaAuxiliaryOutputData(&data->aux_output);
}  /* end FreeAlquimiaData() */
