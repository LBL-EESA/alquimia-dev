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


//------------------------------------------------------------------------
//
//   C utilities for working with alquimia data structures
//
//------------------------------------------------------------------------

#include "alquimia/alquimia_util.h"
#include "alquimia/alquimia_memory.h"
#include "alquimia/alquimia_interface.h"
#include "alquimia/alquimia_constants.h"

//------------------------------------------------------------------------
//
//   Strings
//
//------------------------------------------------------------------------
char* AlquimiaStringDup(const char* str)
{
  int len = strlen(str);
  char* copy = malloc((len+1)*sizeof(char));
  strcpy(copy, str);
  return copy;
}

bool AlquimiaCaseInsensitiveStringCompare(const char* const str1,
                                          const char* const str2) {
  int i;
  bool equal = true;
  if (strlen(str1) != strlen(str2)) {
    equal = false;
  } else {
    for (i = 0; i < (int)strlen(str1); ++i) {
      if (tolower(str1[i]) != tolower(str2[i])) {
        equal = false;
        break;
      }
    }
  }
  return equal;
}  /* end AlquimiaCaseInsensitiveStringCompare() */

//------------------------------------------------------------------------
//
//   Mapping Species names - and indices
//
//------------------------------------------------------------------------
void AlquimiaFindIndexFromName(const char* const name,
                               const AlquimiaVectorString* const names,
                               int* index) {
  int i;
  *index = -1;
  for (i = 0; i < names->size; ++i) {
    if (strncmp(name, names->data[i], (unsigned int)kAlquimiaMaxStringLength) == 0) {
      *index = i;
      break;
    }
  }
}  /* end AlquimiaFindIndexFromName() */


//------------------------------------------------------------------------
//
//   Copying containers
//
//------------------------------------------------------------------------
void CopyAlquimiaVectorDouble(const AlquimiaVectorDouble* const source,
                              AlquimiaVectorDouble* destination)
{
  ResizeAlquimiaVectorDouble(destination, source->size);
  memcpy(destination->data, source->data, sizeof(double) * destination->size);
}

void CopyAlquimiaVectorInt(const AlquimiaVectorInt* const source,
                           AlquimiaVectorInt* destination)
{
  ResizeAlquimiaVectorInt(destination, source->size);
  memcpy(destination->data, source->data, sizeof(int) * destination->size);
}

void CopyAlquimiaVectorString(const AlquimiaVectorString* const source,
                              AlquimiaVectorString* destination)
{
  ResizeAlquimiaVectorString(destination, source->size);
  for (int i = 0; i < destination->size; ++i)
    destination->data[i] = AlquimiaStringDup(source->data[i]);
}

void CopyAlquimiaSizes(const AlquimiaSizes* const source, 
                       AlquimiaSizes* destination)
{
  memcpy(destination, source, sizeof(AlquimiaSizes));
}

void CopyAlquimiaProblemMetaData(const AlquimiaProblemMetaData* const source, 
                                 AlquimiaProblemMetaData* destination)
{
  CopyAlquimiaVectorString(&source->primary_names, &destination->primary_names);
  CopyAlquimiaVectorInt(&source->positivity, &destination->positivity);
  CopyAlquimiaVectorString(&source->mineral_names, &destination->mineral_names);
  CopyAlquimiaVectorString(&source->surface_site_names, &destination->surface_site_names);
  CopyAlquimiaVectorString(&source->ion_exchange_names, &destination->ion_exchange_names);
  CopyAlquimiaVectorString(&source->isotherm_species_names, &destination->isotherm_species_names);
  CopyAlquimiaVectorString(&source->aqueous_kinetic_names, &destination->aqueous_kinetic_names);
}

void CopyAlquimiaProperties(const AlquimiaProperties* const source, 
                            AlquimiaProperties* destination)
{
  destination->volume = source->volume;
  destination->saturation = source->saturation;
  CopyAlquimiaVectorDouble(&source->isotherm_kd, &destination->isotherm_kd);
  CopyAlquimiaVectorDouble(&source->freundlich_n, &destination->freundlich_n);
  CopyAlquimiaVectorDouble(&source->langmuir_b, &destination->langmuir_b);
  CopyAlquimiaVectorDouble(&source->mineral_rate_cnst, &destination->mineral_rate_cnst);
  CopyAlquimiaVectorDouble(&source->aqueous_kinetic_rate_cnst, &destination->aqueous_kinetic_rate_cnst);
}

void CopyAlquimiaEngineFunctionality(const AlquimiaEngineFunctionality* const source, 
                                     AlquimiaEngineFunctionality* destination)
{
  memcpy(destination, source, sizeof(AlquimiaEngineFunctionality));
}

void CopyAlquimiaState(const AlquimiaState* const source, 
                       AlquimiaState* destination)
{
  destination->water_density = source->water_density;
  destination->porosity = source->porosity;
  destination->temperature = source->temperature;
  destination->aqueous_pressure = source->aqueous_pressure;
  CopyAlquimiaVectorDouble(&source->total_mobile, &destination->total_mobile);
  CopyAlquimiaVectorDouble(&source->total_immobile, &destination->total_immobile);
  CopyAlquimiaVectorDouble(&source->mineral_volume_fraction, &destination->mineral_volume_fraction);
  CopyAlquimiaVectorDouble(&source->mineral_specific_surface_area, &destination->mineral_specific_surface_area);
  CopyAlquimiaVectorDouble(&source->surface_site_density, &destination->surface_site_density);
  CopyAlquimiaVectorDouble(&source->cation_exchange_capacity, &destination->cation_exchange_capacity);
}

void CopyAlquimiaAuxiliaryData(const AlquimiaAuxiliaryData* const source, 
                               AlquimiaAuxiliaryData* destination)
{
  CopyAlquimiaVectorInt(&source->aux_ints, &destination->aux_ints);
  CopyAlquimiaVectorDouble(&source->aux_doubles, &destination->aux_doubles);
}

void CopyAlquimiaAuxiliaryOutputData(const AlquimiaAuxiliaryOutputData* const source, 
                                     AlquimiaAuxiliaryOutputData* destination)
{
  destination->pH = source->pH;
  CopyAlquimiaVectorDouble(&source->aqueous_kinetic_rate, &destination->aqueous_kinetic_rate);
  CopyAlquimiaVectorDouble(&source->mineral_saturation_index, &destination->mineral_saturation_index);
  CopyAlquimiaVectorDouble(&source->mineral_reaction_rate, &destination->mineral_reaction_rate);
  CopyAlquimiaVectorDouble(&source->primary_free_ion_concentration, &destination->primary_free_ion_concentration);
  CopyAlquimiaVectorDouble(&source->primary_activity_coeff, &destination->primary_activity_coeff);
  CopyAlquimiaVectorDouble(&source->secondary_free_ion_concentration, &destination->secondary_free_ion_concentration);
  CopyAlquimiaVectorDouble(&source->secondary_activity_coeff, &destination->secondary_activity_coeff);
}

void CopyAlquimiaGeochemicalCondition(const AlquimiaGeochemicalCondition* const source, 
                                      AlquimiaGeochemicalCondition* destination)
{
  if (destination->name != NULL)
    free(destination->name);
  destination->name = AlquimiaStringDup(source->name);
  CopyAlquimiaAqueousConstraintVector(&source->aqueous_constraints, &destination->aqueous_constraints);
  CopyAlquimiaMineralConstraintVector(&source->mineral_constraints, &destination->mineral_constraints);
}

void CopyAlquimiaGeochemicalConditionVector(const AlquimiaGeochemicalConditionVector* source, 
                                            AlquimiaGeochemicalConditionVector* destination)
{
  ResizeAlquimiaGeochemicalConditionVector(destination, source->size);
  for (int i = 0; i < destination->size; ++i)
    CopyAlquimiaGeochemicalCondition(&source->data[i], &destination->data[i]);
}

void CopyAlquimiaAqueousConstraint(const AlquimiaAqueousConstraint* const source, 
                                   AlquimiaAqueousConstraint* destination)
{
  if (destination->primary_species_name != NULL)
    free(destination->primary_species_name);
  destination->primary_species_name = AlquimiaStringDup(source->primary_species_name);
  if (destination->constraint_type != NULL)
    free(destination->constraint_type);
  destination->constraint_type = AlquimiaStringDup(source->constraint_type);
  if (destination->associated_species != NULL)
    free(destination->associated_species);
  destination->associated_species = AlquimiaStringDup(source->associated_species);
  destination->value = source->value;
}

void CopyAlquimiaAqueousConstraintVector(const AlquimiaAqueousConstraintVector* const source, 
                                         AlquimiaAqueousConstraintVector* destination)
{
  ResizeAlquimiaAqueousConstraintVector(destination, source->size);
  for (int i = 0; i < destination->size; ++i)
    CopyAlquimiaAqueousConstraint(&source->data[i], &destination->data[i]);
}

void CopyAlquimiaMineralConstraint(const AlquimiaMineralConstraint* const source, 
                                   AlquimiaMineralConstraint* destination)
{
  if (destination->mineral_name != NULL)
    free(destination->mineral_name);
  destination->mineral_name = AlquimiaStringDup(source->mineral_name);
  destination->volume_fraction = source->volume_fraction;
  destination->specific_surface_area = source->specific_surface_area;
}

void CopyAlquimiaMineralConstraintVector(const AlquimiaMineralConstraintVector* const source, 
                                         AlquimiaMineralConstraintVector* destination)
{
  ResizeAlquimiaMineralConstraintVector(destination, source->size);
  for (int i = 0; i < destination->size; ++i)
    CopyAlquimiaMineralConstraint(&source->data[i], &destination->data[i]);
}

//------------------------------------------------------------------------
//
//   Resizing vectors
//
//------------------------------------------------------------------------
void ResizeAlquimiaVectorDouble(AlquimiaVectorDouble* vec, int new_size)
{
  if (vec->size == 0)
    AllocateAlquimiaVectorDouble(new_size, vec);
  else if (new_size > vec->capacity)
  {
    while (vec->capacity < new_size)
      vec->capacity *= 2;
    vec->data = realloc(vec->data, sizeof(double) * vec->capacity);
  }
  vec->size = new_size;
}

void ResizeAlquimiaVectorInt(AlquimiaVectorInt* vec, int new_size)
{
  if (vec->size == 0)
    AllocateAlquimiaVectorInt(new_size, vec);
  else if (new_size > vec->capacity)
  {
    while (vec->capacity < new_size)
      vec->capacity *= 2;
    vec->data = realloc(vec->data, sizeof(int) * vec->capacity);
  }
  vec->size = new_size;
}

void ResizeAlquimiaVectorString(AlquimiaVectorString* vec, int new_size)
{
  if (vec->size == 0)
    AllocateAlquimiaVectorString(new_size, vec);
  else if (new_size > vec->capacity)
  {
    while (vec->capacity < new_size)
      vec->capacity *= 2;
    vec->data = realloc(vec->data, sizeof(char*) * vec->capacity);
  }
  vec->size = new_size;
}

void ResizeAlquimiaGeochemicalConditionVector(AlquimiaGeochemicalConditionVector* vec, int new_size)
{
  if (vec->size == 0)
    AllocateAlquimiaGeochemicalConditionVector(new_size, vec);
  else if (new_size > vec->capacity)
  {
    while (vec->capacity < new_size)
      vec->capacity *= 2;
    vec->data = realloc(vec->data, sizeof(AlquimiaGeochemicalCondition) * vec->capacity);
  }
  vec->size = new_size;
}

void ResizeAlquimiaAqueousConstraintVector(AlquimiaAqueousConstraintVector* vec, int new_size)
{
  if (vec->size == 0)
    AllocateAlquimiaAqueousConstraintVector(new_size, vec);
  else if (new_size > vec->capacity)
  {
    while (vec->capacity < new_size)
      vec->capacity *= 2;
    vec->data = realloc(vec->data, sizeof(AlquimiaAqueousConstraint) * vec->capacity);
  }
  vec->size = new_size;
}

void ResizeAlquimiaMineralConstraintVector(AlquimiaMineralConstraintVector* vec, int new_size)
{
  if (vec->size == 0)
    AllocateAlquimiaMineralConstraintVector(new_size, vec);
  else if (new_size > vec->capacity)
  {
    while (vec->capacity < new_size)
      vec->capacity *= 2;
    vec->data = realloc(vec->data, sizeof(AlquimiaMineralConstraint) * vec->capacity);
  }
  vec->size = new_size;
}

//------------------------------------------------------------------------
//
//   Printing vectors
//
//------------------------------------------------------------------------
void PrintAlquimiaVectorDouble(const char* const name,
                               const AlquimiaVectorDouble* const vector,
                               FILE* file) {
  int i;
  fprintf(file, "    %s (%d) (%p):\n", name, vector->size, (void*)(&vector->data));
  fprintf(file, "   [ ");
  for (i = 0; i < vector->size; ++i) {
    fprintf(file, "%e, ", vector->data[i]);
  }
  fprintf(file, "]\n");
}  /* end PrintAlqumiaVectorDouble() */

void PrintAlquimiaVectorInt(const char* const name,
                            const AlquimiaVectorInt* const vector,
                            FILE* file) {
  int i;
  fprintf(file, "    %s (%d) (%p):\n", name, vector->size, (void*)(&vector->data));
  fprintf(file, "   [ ");
  for (i = 0; i < vector->size; ++i) {
    fprintf(file, "%d, ", vector->data[i]);
  }
  fprintf(file, "]\n");
}  /* end PrintAlqumiaVectorInt() */

void PrintAlquimiaVectorString(const char* const name,
                               const AlquimiaVectorString* const vector,
                               FILE* file) {
  int i;
  fprintf(file, "    %s (%d) (%p):\n", name, vector->size, (void*)(&vector->data));
  fprintf(file, "   [ ");
  for (i = 0; i < vector->size; ++i) {
    fprintf(file, "'%s', ", vector->data[i]);
  }
  fprintf(file, "]\n");
}  /* end PrintAlqumiaVectorInt() */


//------------------------------------------------------------------------
//
//   Printing containers
//
//------------------------------------------------------------------------
void PrintAlquimiaData(const AlquimiaData* const data, FILE* file) {
  fprintf(file, "- Alquimia Data ----------------------------------------\n");
  fprintf(file, "  engine_state : %p\n", data->engine_state);
  PrintAlquimiaSizes(&data->sizes, file);
  PrintAlquimiaEngineFunctionality(&data->functionality, file);
  PrintAlquimiaState(&data->state, file);
  PrintAlquimiaProperties(&data->properties, file);
  PrintAlquimiaAuxiliaryData(&data->aux_data, file);
  PrintAlquimiaProblemMetaData(&data->meta_data, file);
  PrintAlquimiaAuxiliaryOutputData(&data->aux_output, file);
  fprintf(file, "---------------------------------------- Alquimia Data -\n");
}  /* end PrintAlquimiaData() */

void PrintAlquimiaSizes(const AlquimiaSizes* const sizes, FILE* file) {
  fprintf(file, "-- Alquimia Sizes :\n");
  fprintf(file, "     num primary species : %d\n", sizes->num_primary);
  fprintf(file, "     num sorbed : %d\n", sizes->num_sorbed);
  fprintf(file, "     num minerals : %d\n", sizes->num_minerals);
  fprintf(file, "     num aqueous complexes : %d\n", sizes->num_aqueous_complexes);
  fprintf(file, "     num aqueous kinetics : %d\n", sizes->num_aqueous_kinetics);
  fprintf(file, "     num surface sites : %d\n", sizes->num_surface_sites);
  fprintf(file, "     num ion exchange sites : %d\n", sizes->num_ion_exchange_sites);
  fprintf(file, "     num auxiliary integers : %d\n", sizes->num_aux_integers);
  fprintf(file, "     num auxiliary doubles : %d\n", sizes->num_aux_doubles);
}  /* end PrintAlquimiaSizes() */

void PrintAlquimiaEngineFunctionality(const AlquimiaEngineFunctionality* const functionality, FILE* file) {

  fprintf(file, "-- Alquimia Engine Functionality :\n");
  fprintf(file, "     thread_safe : %d\n", functionality->thread_safe);
  fprintf(file, "     temperature_dependent : %d\n",
          functionality->temperature_dependent);
  fprintf(file, "     pressure_dependent : %d\n", 
          functionality->pressure_dependent);
  fprintf(file, "     porosity_update  : %d\n", functionality->porosity_update);
  fprintf(file, "     index base : %d\n", functionality->index_base);
}  /* end PrintAlquimiaEngineFunctionality() */

void PrintAlquimiaProblemMetaData(const AlquimiaProblemMetaData* const meta_data, FILE* file) {

  fprintf(file, "-- Alquimia Problem Meta Data :\n");
  PrintAlquimiaVectorString("primary names", &(meta_data->primary_names), file);
  PrintAlquimiaVectorInt("positivity names", &(meta_data->positivity), file);
  PrintAlquimiaVectorString("mineral names", &(meta_data->mineral_names), file);
  PrintAlquimiaVectorString("surface site names", &(meta_data->surface_site_names), file);
  PrintAlquimiaVectorString("ion exchange names", &(meta_data->ion_exchange_names), file);
  PrintAlquimiaVectorString("isotherm species names", &(meta_data->isotherm_species_names), file);
  PrintAlquimiaVectorString("aqueous kinetic names", &(meta_data->aqueous_kinetic_names), file);
}  /* end PrintAlquimiaProblemMetaData() */

void PrintAlquimiaProperties(const AlquimiaProperties* const mat_prop, FILE* file) {

  fprintf(file, "-- Alquimia Properties :\n");
  fprintf(file, "     volume : %f\n", mat_prop->volume);
  fprintf(file, "     saturation : %f\n", mat_prop->saturation);
  PrintAlquimiaVectorDouble("isotherm kd", &(mat_prop->isotherm_kd), file);
  PrintAlquimiaVectorDouble("freundlich n", &(mat_prop->freundlich_n), file);
  PrintAlquimiaVectorDouble("langmuir b", &(mat_prop->langmuir_b), file);
  PrintAlquimiaVectorDouble("mineral rate cnst", &(mat_prop->mineral_rate_cnst), file);
  PrintAlquimiaVectorDouble("aqueous kinetic rate cnst", &(mat_prop->aqueous_kinetic_rate_cnst), file);
}  /* end PrintAlquimiaProperties() */

void PrintAlquimiaState(const AlquimiaState* const state, FILE* file) {

  fprintf(file, "-- Alquimia State:\n");
  fprintf(file, "     water density : %f\n", state->water_density);
  fprintf(file, "     porosity : %f\n", state->porosity);
  fprintf(file, "     temperature : %f\n", state->temperature);
  fprintf(file, "     aqueous_pressure : %f\n", state->aqueous_pressure);

  PrintAlquimiaVectorDouble("total_mobile", &(state->total_mobile), file);
  PrintAlquimiaVectorDouble("total_immobile", &(state->total_immobile), file);
  PrintAlquimiaVectorDouble("kinetic minerals volume fraction",
                            &(state->mineral_volume_fraction), file);
  PrintAlquimiaVectorDouble("kinetic minerals specific surface area",
                            &(state->mineral_specific_surface_area), file);
  PrintAlquimiaVectorDouble("cation_exchange_capacity",
                            &(state->cation_exchange_capacity), file);
  PrintAlquimiaVectorDouble("surface_site_density",
                            &(state->surface_site_density), file);
}  /* end PrintAlquimiaState() */

void PrintAlquimiaAuxiliaryData(const AlquimiaAuxiliaryData* const aux_data, FILE* file) {

  fprintf(file, "-- Alquimia Auxiliary Data:\n");
  PrintAlquimiaVectorInt("auxiliary integers", &(aux_data->aux_ints), file);
  PrintAlquimiaVectorDouble("auxiliary doubles", &(aux_data->aux_doubles), file);
}  /* end PrintAlquimiaAuxiliaryData() */

void PrintAlquimiaAuxiliaryOutputData(const AlquimiaAuxiliaryOutputData* const aux_output, FILE* file) {

  fprintf(file, "-- Alquimia Auxiliary Output Data:\n");
  fprintf(file, "     pH : %f\n", aux_output->pH);

  PrintAlquimiaVectorDouble("mineral saturation index",
                            &(aux_output->mineral_saturation_index), file);
  PrintAlquimiaVectorDouble("mineral reaction rate",
                            &(aux_output->mineral_reaction_rate), file);
  PrintAlquimiaVectorDouble("primary free ion concentrations",
                            &(aux_output->primary_free_ion_concentration), file);
  PrintAlquimiaVectorDouble("primary activity coeff",
                            &(aux_output->primary_activity_coeff), file);
  PrintAlquimiaVectorDouble("secondary free ion concentrations",
                            &(aux_output->secondary_free_ion_concentration), file);
  PrintAlquimiaVectorDouble("secondary activity coeff",
                            &(aux_output->secondary_activity_coeff), file);
}  /* end PrintAlquimiaAuxiliaryOutputData() */

void PrintAlquimiaGeochemicalConditionVector(const AlquimiaGeochemicalConditionVector* const condition_list, FILE* file) {
  int i;
  fprintf(file, "- Alquimia Geochemical Condition List ------------------\n");
  for (i = 0; i < condition_list->size; ++i) {
    PrintAlquimiaGeochemicalCondition(&(condition_list->data[i]), file);
    fprintf(file, "\n");
  }
  fprintf(file, "------------------ Alquimia Geochemical Condition List -\n");
}  /*  PrintAlquimiaGeochemicalConditionVector() */

void PrintAlquimiaGeochemicalCondition(const AlquimiaGeochemicalCondition* const condition, FILE* file) {
  int i;
  fprintf(file, "-- Alquimia Geochemical Condition : %s\n", condition->name);
  for (i = 0; i < condition->aqueous_constraints.size; ++i) {
    PrintAlquimiaAqueousConstraint(&(condition->aqueous_constraints.data[i]), file);
  }
  for (i = 0; i < condition->mineral_constraints.size; ++i) {
    PrintAlquimiaMineralConstraint(&(condition->mineral_constraints.data[i]), file);
  }
  fprintf(file, "\n");
}  /*  PrintAlquimiaGeochemicalCondition() */

void PrintAlquimiaAqueousConstraint(const AlquimiaAqueousConstraint* const constraint, FILE* file) {
  fprintf(file, "--- Alquimia Aqueous Constraint : \n");
  fprintf(file, "      primary species : %s\n", constraint->primary_species_name);
  fprintf(file, "      constraint type : %s\n", constraint->constraint_type);
  fprintf(file, "      associated species : %s\n", constraint->associated_species);
  fprintf(file, "      value : %e\n", constraint->value);
}  /*  PrintAlquimiaAqueousConstraint() */

void PrintAlquimiaMineralConstraint(const AlquimiaMineralConstraint* const constraint, FILE* file) {
  fprintf(file, "--- Alquimia Mineral Constraint : \n");
  fprintf(file, "      mineral : %s\n", constraint->mineral_name);
  fprintf(file, "      volume fraction : %e\n", constraint->volume_fraction);
  fprintf(file, "      specific surface area : %e\n", constraint->specific_surface_area);
}  /*  PrintAlquimiaMineralConstraint() */
