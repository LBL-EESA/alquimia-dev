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

#include "alquimia_containers.h"

/*******************************************************************************
 **
 **  Strings
 **
 *******************************************************************************/
bool AlquimiaCaseInsensitiveStringCompare(const char* const str1,
                                          const char* const str2) {
  bool equal = true;
  if (strlen(str1) != strlen(str2)) {
    equal = false;
  } else {
    for (size_t i = 0; i < strlen(str1); ++i) {
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
 **  Printing
 **
 *******************************************************************************/
void PrintAlquimiaDoubleArray(const char* name, const int size, const double* array) {
  fprintf(stdout, "    %s (%d) (%p):\n", name, size, array);
  fprintf(stdout, "   [ ");
  for (int i = 0; i < size; ++i) {
    fprintf(stdout, "%e, ", array[i]);
  }
  fprintf(stdout, "]\n");
}  // end PrintAlquimiaDoubleArray()

void PrintAlquimiaIntArray(const char* name, const int size, const int* array) {
  fprintf(stdout, "    %s (%d) (%p):\n", name, size, array);
  fprintf(stdout, "   [ ");
  for (int i = 0; i < size; ++i) {
    fprintf(stdout, "%d, ", array[i]);
  }
  fprintf(stdout, "]\n");
}  // end PrintAlqumiaIntArray()

void PrintAlquimiaStringArray(const char* name, const int size,
                              char** array) {
  fprintf(stdout, "    %s (%d) (%p):\n", name, size, array);
  fprintf(stdout, "   [ ");
  for (int i = 0; i < size; ++i) {
    fprintf(stdout, "'%s', ", array[i]);
  }
  fprintf(stdout, "]\n");
}  // end PrintAlquimiaStringArray()


void PrintAlquimiaSizes(const struct AlquimiaSizes* sizes) {
  fprintf(stdout, "  Alquimia Sizes :\n");
  fprintf(stdout, "    num primary species : %d\n", sizes->num_primary);
  fprintf(stdout, "    num sorbed : %d\n", sizes->num_sorbed);
  fprintf(stdout, "    num kinetic minerals : %d\n", sizes->num_kinetic_minerals);
  fprintf(stdout, "    num aqueous complexes : %d\n", sizes->num_aqueous_complexes);
  fprintf(stdout, "    num surface sites : %d\n", sizes->num_surface_sites);
  fprintf(stdout, "    num ion exchange sites : %d\n", sizes->num_ion_exchange_sites);
}  // end PrintAlquimiaSizes()

void PrintAlquimiaMetaData(const struct AlquimiaMetaData* meta_data) {

  fprintf(stdout, "  Alquimia Meta Data :\n");
  fprintf(stdout, "    thread_safe : %d\n", meta_data->thread_safe);
  fprintf(stdout, "    temperature_dependent : %d\n",
          meta_data->temperature_dependent);
  fprintf(stdout, "    pressure_dependent : %d\n", meta_data->pressure_dependent);
  fprintf(stdout, "    porosity_update  : %d\n", meta_data->porosity_update);
  fprintf(stdout, "    index base : %d\n", meta_data->index_base);
  PrintAlquimiaIntArray("primary indices", meta_data->size_primary,
                        meta_data->primary_indices);
  PrintAlquimiaStringArray("primary names", meta_data->size_primary,
                           meta_data->primary_names);
  PrintAlquimiaIntArray("mineral indices",
                        meta_data->size_minerals, meta_data->mineral_indices);
  PrintAlquimiaStringArray("mineral names",
                           meta_data->size_minerals, meta_data->mineral_names);
}  // end PrintAlquimiaMetaData()

void PrintAlquimiaState(const struct AlquimiaState* state) {

  fprintf(stdout, "  Alquimia State:\n");
  fprintf(stdout, "    water density : %f\n", state->water_density);
  fprintf(stdout, "    saturation : %f\n", state->saturation);
  fprintf(stdout, "    porosity : %f\n", state->porosity);
  fprintf(stdout, "    temperature : %f\n", state->temperature);
  fprintf(stdout, "    aqueous_pressure : %f\n", state->aqueous_pressure);

  PrintAlquimiaDoubleArray("total_primary",
                           state->size_total_primary, state->total_primary);
  PrintAlquimiaDoubleArray("free_ion", state->size_free_ion,
                           state->free_ion);
  PrintAlquimiaDoubleArray("kinetic minerals volume fraction",
                           state->size_mineral_volume_fraction,
                           state->mineral_volume_fraction);
  PrintAlquimiaDoubleArray("kinetic minerals specific surface area",
                           state->size_mineral_specific_surface_area,
                           state->mineral_specific_surface_area);
}  // end PrintAlquimiaState()

void PrintAlquimiaAuxiliaryData(const struct AlquimiaAuxiliaryData* aux_data) {

  fprintf(stdout, "  Alquimia Auxiliary Data:\n");
  PrintAlquimiaDoubleArray("primary activity coeff",
                           aux_data->size_primary_activity_coeff,
                           aux_data->primary_activity_coeff);
  PrintAlquimiaDoubleArray("secondary activity coeff",
                           aux_data->size_secondary_activity_coeff,
                           aux_data->secondary_activity_coeff);
  PrintAlquimiaDoubleArray("ion exchange ref cation conc",
                           aux_data->size_ion_exchange_ref_cation_conc,
                           aux_data->ion_exchange_ref_cation_conc);
  PrintAlquimiaDoubleArray("surface complex free site conc",
                           aux_data->size_surface_complex_free_site_conc,
                           aux_data->surface_complex_free_site_conc);
}  // end PrintAlquimiaAuxiliaryData()

void PrintAlquimiaAuxiliaryOutputData(
    const struct AlquimiaAuxiliaryOutputData* aux_output) {

  fprintf(stdout, "  Alquimia Auxiliary Output Data:\n");
  fprintf(stdout, "    pH : %f\n", aux_output->pH);

  PrintAlquimiaDoubleArray("mineral saturation index", aux_output->size_minerals,
                           aux_output->mineral_saturation_index);
  PrintAlquimiaDoubleArray("mineral reaction rate", aux_output->size_minerals,
                           aux_output->mineral_reaction_rate);
}  // end PrintAlquimiaAuxiliaryOutputData()

void PrintAlquimiaGeochemicalConditionList(
    const struct AlquimiaGeochemicalConditionList* condition_list) {

  fprintf(stdout, "Alquimia Geochemical Condition List : \n");
  for (int i = 0; i < condition_list->num_conditions; ++i) {
    PrintAlquimiaGeochemicalCondition(&(condition_list->conditions[i]));
    fprintf(stdout, "\n");
  }
}  //  PrintAlquimiaGeochemicalConditionList()

void PrintAlquimiaGeochemicalCondition(
    const struct AlquimiaGeochemicalCondition* condition) {

  fprintf(stdout, "  Alquimia Geochemical Condition : %s\n", condition->name);
  for (int i = 0; i < condition->num_constraints; ++i) {
    PrintAlquimiaGeochemicalConstraint(&(condition->constraints[i]));
    fprintf(stdout, "\n");
  }
}  //  PrintAlquimiaGeochemicalCondition()

void PrintAlquimiaGeochemicalConstraint(
    const struct AlquimiaGeochemicalConstraint* constraint) {
  fprintf(stdout, "    Alquimia Geochemical Constraint : \n");
  fprintf(stdout, "      primary species : %s\n", constraint->primary_species);
  fprintf(stdout, "      constraint type : %s\n", constraint->constraint_type);
  fprintf(stdout, "      associated species : %s\n", constraint->associated_species);
  fprintf(stdout, "      value : %e\n", constraint->value);
}  //  PrintAlquimiaGeochemicalConstraint()
