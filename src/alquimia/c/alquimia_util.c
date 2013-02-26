/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*******************************************************************************
 **
 **  C utilities for working with alquimia data structures
 **
 *******************************************************************************/

#include "alquimia_util.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#include "alquimia_containers.h"
#include "alquimia_interface.h"
#include "alquimia_constants.h"

/*******************************************************************************
 **
 **  Strings
 **
 *******************************************************************************/
bool AlquimiaCaseInsensitiveStringCompare(const char* const str1,
                                          const char* const str2) {
  int i;
  bool equal = true;
  if (strlen(str1) != strlen(str2)) {
    equal = false;
  } else {
    for (i = 0; i < strlen(str1); ++i) {
      if (tolower(str1[i]) != tolower(str2[i])) {
        equal = false;
        break;
      }  // end if()
    }  // end for()
  } // end else()
  return equal;
}  // end AlquimiaCaseInsensitiveStringCompare()

/*******************************************************************************
 **
 **  Mapping Species names - and indices
 **
 *******************************************************************************/
void AlquimiaIndexFromName(const char* const name,
                           const struct AlquimiaVectorString* const names,
                           const struct AlquimiaVectorInt* const indices,
                           int* index) {
  int i;
  *index = -1;
  assert(names->size == indices->size);
  for (i = 0; i < names->size; ++i) {
    if (strncmp(name, names->data[i], kAlquimiaMaxStringLength)) {
      *index = i;
      break;
    }
  }
}  // end AlquimiaIndexFromName()

void AlquimiaNameFromIndex(const int index,
                           const struct AlquimiaVectorString* const names,
                           const struct AlquimiaVectorInt* const indices,
                           char* name) {
  int i;
  assert(names->size == indices->size);
  assert(index >= 0);
  assert(index <= names->size);
  for (i = 0; i < names->size; ++i) {
    if (index == indices->data[i]) {
      strncpy(name, names->data[i], kAlquimiaMaxStringLength);
      break;
    }
  }
}  // end AlquimiaNameFromIndex()


/*******************************************************************************
 **
 **  Printing Vectors
 **
 *******************************************************************************/
void PrintAlquimiaVectorDouble(const char* const name,
                               const struct AlquimiaVectorDouble* const vector) {
  int i;
  fprintf(stdout, "    %s (%d) (%p):\n", name, vector->size, &(vector->data));
  fprintf(stdout, "   [ ");
  for (i = 0; i < vector->size; ++i) {
    fprintf(stdout, "%e, ", vector->data[i]);
  }
  fprintf(stdout, "]\n");
}  // end PrintAlqumiaVectorDouble()

void PrintAlquimiaVectorInt(const char* const name,
                            const struct AlquimiaVectorInt* const vector) {
  int i;
  fprintf(stdout, "    %s (%d) (%p):\n", name, vector->size, &(vector->data));
  fprintf(stdout, "   [ ");
  for (i = 0; i < vector->size; ++i) {
    fprintf(stdout, "%d, ", vector->data[i]);
  }
  fprintf(stdout, "]\n");
}  // end PrintAlqumiaVectorInt()

void PrintAlquimiaVectorString(const char* const name,
                               const struct AlquimiaVectorString* const vector) {
  int i;
  fprintf(stdout, "    %s (%d) (%p):\n", name, vector->size, &(vector->data));
  fprintf(stdout, "   [ ");
  for (i = 0; i < vector->size; ++i) {
    fprintf(stdout, "'%s', ", vector->data[i]);
  }
  fprintf(stdout, "]\n");
}  // end PrintAlqumiaVectorInt()


/*******************************************************************************
 **
 **  Printing Containers
 **
 *******************************************************************************/
void PrintAlquimiaData(const struct AlquimiaData* const data) {
  fprintf(stdout, "- Alquimia Data ----------------------------------------\n");
  fprintf(stdout, "  engine_state : %p\n", data->engine_state);
  PrintAlquimiaSizes(&data->sizes);
  PrintAlquimiaEngineFunctionality(&data->functionality);
  PrintAlquimiaState(&data->state);
  PrintAlquimiaMaterialProperties(&data->material_properties);
  PrintAlquimiaAuxiliaryData(&data->aux_data);
  PrintAlquimiaProblemMetaData(&data->meta_data);
  PrintAlquimiaAuxiliaryOutputData(&data->aux_output);
  fprintf(stdout, "---------------------------------------- Alquimia Data -\n");
}  // end PrintAlquimiaData()

void PrintAlquimiaSizes(const struct AlquimiaSizes* const sizes) {
  fprintf(stdout, "-- Alquimia Sizes :\n");
  fprintf(stdout, "     num primary species : %d\n", sizes->num_primary);
  fprintf(stdout, "     num sorbed : %d\n", sizes->num_sorbed);
  fprintf(stdout, "     num kinetic minerals : %d\n", sizes->num_kinetic_minerals);
  fprintf(stdout, "     num aqueous complexes : %d\n", sizes->num_aqueous_complexes);
  fprintf(stdout, "     num surface sites : %d\n", sizes->num_surface_sites);
  fprintf(stdout, "     num ion exchange sites : %d\n", sizes->num_ion_exchange_sites);
  fprintf(stdout, "     num auxiliary integers : %d\n", sizes->num_aux_integers);
  fprintf(stdout, "     num auxiliary doubles : %d\n", sizes->num_aux_doubles);
}  // end PrintAlquimiaSizes()

void PrintAlquimiaEngineFunctionality(
    const struct AlquimiaEngineFunctionality* const functionality) {

  fprintf(stdout, "-- Alquimia Engine Functionality :\n");
  fprintf(stdout, "     thread_safe : %d\n", functionality->thread_safe);
  fprintf(stdout, "     temperature_dependent : %d\n",
          functionality->temperature_dependent);
  fprintf(stdout, "     pressure_dependent : %d\n", 
          functionality->pressure_dependent);
  fprintf(stdout, "     porosity_update  : %d\n", functionality->porosity_update);
  fprintf(stdout, "     index base : %d\n", functionality->index_base);
}  // end PrintAlquimiaEngineFunctionality()

void PrintAlquimiaProblemMetaData(const struct AlquimiaProblemMetaData* const meta_data) {

  fprintf(stdout, "-- Alquimia Problem Meta Data :\n");
  PrintAlquimiaVectorInt("primary indices", &(meta_data->primary_indices));
  PrintAlquimiaVectorString("primary names", &(meta_data->primary_names));
  PrintAlquimiaVectorInt("mineral indices", &(meta_data->mineral_indices));
  PrintAlquimiaVectorString("mineral names", &(meta_data->mineral_names));
  PrintAlquimiaVectorInt("surface site indices", &(meta_data->surface_site_indices));
  PrintAlquimiaVectorString("surface site names", &(meta_data->surface_site_names));
  PrintAlquimiaVectorInt("isotherm species indicies", &(meta_data->isotherm_species_indices));
}  // end PrintAlquimiaProblemMetaData()

void PrintAlquimiaMaterialProperties(const struct AlquimiaMaterialProperties* const mat_prop) {

  fprintf(stdout, "-- Alquimia Material Properties :\n");
  PrintAlquimiaVectorDouble("isotherm kd", &(mat_prop->isotherm_kd));
  PrintAlquimiaVectorDouble("freundlich n", &(mat_prop->freundlich_n));
  PrintAlquimiaVectorDouble("langmuir b", &(mat_prop->langmuir_b));
}  // end PrintAlquimiaMaterialProperties()

void PrintAlquimiaState(const struct AlquimiaState* const state) {

  fprintf(stdout, "-- Alquimia State:\n");
  fprintf(stdout, "     water density : %f\n", state->water_density);
  fprintf(stdout, "     saturation : %f\n", state->saturation);
  fprintf(stdout, "     porosity : %f\n", state->porosity);
  fprintf(stdout, "     temperature : %f\n", state->temperature);
  fprintf(stdout, "     aqueous_pressure : %f\n", state->aqueous_pressure);

  PrintAlquimiaVectorDouble("total_mobile", &(state->total_mobile));
  PrintAlquimiaVectorDouble("total_immobile", &(state->total_immobile));
  PrintAlquimiaVectorDouble("kinetic minerals volume fraction",
                            &(state->mineral_volume_fraction));
  PrintAlquimiaVectorDouble("kinetic minerals specific surface area",
                            &(state->mineral_specific_surface_area));
  PrintAlquimiaVectorDouble("cation_exchange_capacity",
                            &(state->cation_exchange_capacity));
  PrintAlquimiaVectorDouble("surface_site_density",
                            &(state->surface_site_density));
}  // end PrintAlquimiaState()

void PrintAlquimiaAuxiliaryData(const struct AlquimiaAuxiliaryData* const aux_data) {

  fprintf(stdout, "-- Alquimia Auxiliary Data:\n");
  PrintAlquimiaVectorInt("auxiliary integers", &(aux_data->aux_ints));
  PrintAlquimiaVectorDouble("auxiliary doubles", &(aux_data->aux_doubles));
}  // end PrintAlquimiaAuxiliaryData()

void PrintAlquimiaAuxiliaryOutputData(
    const struct AlquimiaAuxiliaryOutputData* const aux_output) {

  fprintf(stdout, "-- Alquimia Auxiliary Output Data:\n");
  fprintf(stdout, "     pH : %f\n", aux_output->pH);

  PrintAlquimiaVectorDouble("mineral saturation index",
                            &(aux_output->mineral_saturation_index));
  PrintAlquimiaVectorDouble("mineral reaction rate",
                            &(aux_output->mineral_reaction_rate));
}  // end PrintAlquimiaAuxiliaryOutputData()

void PrintAlquimiaGeochemicalConditionVector(
    const struct AlquimiaGeochemicalConditionVector* const condition_list) {
  int i;
  fprintf(stdout, "- Alquimia Geochemical Condition List ------------------\n");
  for (i = 0; i < condition_list->size; ++i) {
    PrintAlquimiaGeochemicalCondition(&(condition_list->data[i]));
    fprintf(stdout, "\n");
  }
  fprintf(stdout, "------------------ Alquimia Geochemical Condition List -\n");
}  //  PrintAlquimiaGeochemicalConditionVector()

void PrintAlquimiaGeochemicalCondition(
    const struct AlquimiaGeochemicalCondition* const condition) {
  int i;
  fprintf(stdout, "-- Alquimia Geochemical Condition : %s\n", condition->name);
  for (i = 0; i < condition->aqueous_constraints.size; ++i) {
    PrintAlquimiaAqueousConstraint(&(condition->aqueous_constraints.data[i]));
  }
  for (i = 0; i < condition->mineral_constraints.size; ++i) {
    PrintAlquimiaMineralConstraint(&(condition->mineral_constraints.data[i]));
  }
  fprintf(stdout, "\n");
}  //  PrintAlquimiaGeochemicalCondition()

void PrintAlquimiaAqueousConstraint(
    const struct AlquimiaAqueousConstraint* const constraint) {
  fprintf(stdout, "--- Alquimia Aqueous Constraint : \n");
  fprintf(stdout, "      primary species : %s\n", constraint->primary_species_name);
  fprintf(stdout, "      constraint type : %s\n", constraint->constraint_type);
  fprintf(stdout, "      associated species : %s\n", constraint->associated_species);
  fprintf(stdout, "      value : %e\n", constraint->value);
}  //  PrintAlquimiaAqueousConstraint()

void PrintAlquimiaMineralConstraint(
    const struct AlquimiaMineralConstraint* const constraint) {
  fprintf(stdout, "--- Alquimia Mineral Constraint : \n");
  fprintf(stdout, "      mineral : %s\n", constraint->mineral_name);
  fprintf(stdout, "      volume fraction : %e\n", constraint->volume_fraction);
  fprintf(stdout, "      specific surface area : %e\n", constraint->specific_surface_area);
}  //  PrintAlquimiaMineralConstraint()
