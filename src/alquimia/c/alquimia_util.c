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
 **  Printing Vectors
 **
 *******************************************************************************/
void PrintAlquimiaVectorDouble(const char* const name,
                               const struct AlquimiaVectorDouble* const vector) {
  fprintf(stdout, "    %s (%d) (%p):\n", name, vector->size, &(vector->data));
  fprintf(stdout, "   [ ");
  for (int i = 0; i < vector->size; ++i) {
    fprintf(stdout, "%e, ", vector->data[i]);
  }
  fprintf(stdout, "]\n");
}  // end PrintAlqumiaVectorDouble()

void PrintAlquimiaVectorInt(const char* const name,
                            const struct AlquimiaVectorInt* const vector) {
  fprintf(stdout, "    %s (%d) (%p):\n", name, vector->size, &(vector->data));
  fprintf(stdout, "   [ ");
  for (int i = 0; i < vector->size; ++i) {
    fprintf(stdout, "%d, ", vector->data[i]);
  }
  fprintf(stdout, "]\n");
}  // end PrintAlqumiaVectorInt()

void PrintAlquimiaVectorString(const char* const name,
                               const struct AlquimiaVectorString* const vector) {
  fprintf(stdout, "    %s (%d) (%p):\n", name, vector->size, &(vector->data));
  fprintf(stdout, "   [ ");
  for (int i = 0; i < vector->size; ++i) {
    fprintf(stdout, "'%s', ", vector->data[i]);
  }
  fprintf(stdout, "]\n");
}  // end PrintAlqumiaVectorInt()


/*******************************************************************************
 **
 **  Printing Containers
 **
 *******************************************************************************/
void PrintAlquimiaSizes(const struct AlquimiaSizes* const sizes) {
  fprintf(stdout, "  Alquimia Sizes :\n");
  fprintf(stdout, "    num primary species : %d\n", sizes->num_primary);
  fprintf(stdout, "    num sorbed : %d\n", sizes->num_sorbed);
  fprintf(stdout, "    num kinetic minerals : %d\n", sizes->num_kinetic_minerals);
  fprintf(stdout, "    num aqueous complexes : %d\n", sizes->num_aqueous_complexes);
  fprintf(stdout, "    num surface sites : %d\n", sizes->num_surface_sites);
  fprintf(stdout, "    num ion exchange sites : %d\n", sizes->num_ion_exchange_sites);
  fprintf(stdout, "    num auxiliary integers : %d\n", sizes->num_aux_integers);
  fprintf(stdout, "    num auxiliary doubles : %d\n", sizes->num_aux_doubles);
}  // end PrintAlquimiaSizes()

void PrintAlquimiaEngineFunctionality(
    const struct AlquimiaEngineFunctionality* const functionality) {

  fprintf(stdout, "  Alquimia Engine Functionality :\n");
  fprintf(stdout, "    thread_safe : %d\n", functionality->thread_safe);
  fprintf(stdout, "    temperature_dependent : %d\n",
          functionality->temperature_dependent);
  fprintf(stdout, "    pressure_dependent : %d\n", 
          functionality->pressure_dependent);
  fprintf(stdout, "    porosity_update  : %d\n", functionality->porosity_update);
  fprintf(stdout, "    index base : %d\n", functionality->index_base);
}  // end PrintAlquimiaEngineFunctionality()

void PrintAlquimiaProblemMetaData(const struct AlquimiaProblemMetaData* const meta_data) {

  fprintf(stdout, "  Alquimia Problem Meta Data :\n");
  PrintAlquimiaVectorInt("primary indices", &(meta_data->primary_indices));
  PrintAlquimiaVectorString("primary names", &(meta_data->primary_names));
  PrintAlquimiaVectorInt("mineral indices", &(meta_data->mineral_indices));
  PrintAlquimiaVectorString("mineral names", &(meta_data->mineral_names));
}  // end PrintAlquimiaProblemMetaData()

void PrintAlquimiaState(const struct AlquimiaState* const state) {

  fprintf(stdout, "  Alquimia State:\n");
  fprintf(stdout, "    water density : %f\n", state->water_density);
  fprintf(stdout, "    saturation : %f\n", state->saturation);
  fprintf(stdout, "    porosity : %f\n", state->porosity);
  fprintf(stdout, "    temperature : %f\n", state->temperature);
  fprintf(stdout, "    aqueous_pressure : %f\n", state->aqueous_pressure);

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

  fprintf(stdout, "  Alquimia Auxiliary Data:\n");
  PrintAlquimiaVectorInt("auxiliary integers", &(aux_data->aux_ints));
  PrintAlquimiaVectorDouble("auxiliary doubles", &(aux_data->aux_doubles));
}  // end PrintAlquimiaAuxiliaryData()

void PrintAlquimiaAuxiliaryOutputData(
    const struct AlquimiaAuxiliaryOutputData* const aux_output) {

  fprintf(stdout, "  Alquimia Auxiliary Output Data:\n");
  fprintf(stdout, "    pH : %f\n", aux_output->pH);

  PrintAlquimiaVectorDouble("mineral saturation index",
                            &(aux_output->mineral_saturation_index));
  PrintAlquimiaVectorDouble("mineral reaction rate",
                            &(aux_output->mineral_reaction_rate));
}  // end PrintAlquimiaAuxiliaryOutputData()

void PrintAlquimiaGeochemicalConditionVector(
    const struct AlquimiaGeochemicalConditionVector* const condition_list) {

  fprintf(stdout, "Alquimia Geochemical Condition List : \n");
  for (int i = 0; i < condition_list->size; ++i) {
    PrintAlquimiaGeochemicalCondition(&(condition_list->data[i]));
    fprintf(stdout, "\n");
  }
}  //  PrintAlquimiaGeochemicalConditionVector()

void PrintAlquimiaGeochemicalCondition(
    const struct AlquimiaGeochemicalCondition* const condition) {

  fprintf(stdout, "  Alquimia Geochemical Condition : %s\n", condition->name);
  for (int i = 0; i < condition->aqueous_constraints.size; ++i) {
    PrintAlquimiaAqueousConstraint(&(condition->aqueous_constraints.data[i]));
  }
  for (int i = 0; i < condition->mineral_constraints.size; ++i) {
    PrintAlquimiaMineralConstraint(&(condition->mineral_constraints.data[i]));
  }
  fprintf(stdout, "\n");
}  //  PrintAlquimiaGeochemicalCondition()

void PrintAlquimiaAqueousConstraint(
    const struct AlquimiaAqueousConstraint* const constraint) {
  fprintf(stdout, "    Alquimia Aqueous Constraint : \n");
  fprintf(stdout, "      primary species : %s\n", constraint->primary_species_name);
  fprintf(stdout, "      constraint type : %s\n", constraint->constraint_type);
  fprintf(stdout, "      associated species : %s\n", constraint->associated_species);
  fprintf(stdout, "      value : %e\n", constraint->value);
}  //  PrintAlquimiaAqueousConstraint()

void PrintAlquimiaMineralConstraint(
    const struct AlquimiaMineralConstraint* const constraint) {
  fprintf(stdout, "    Alquimia Mineral Constraint : \n");
  fprintf(stdout, "      mineral : %s\n", constraint->mineral_name);
  fprintf(stdout, "      volume fraction : %e\n", constraint->volume_fraction);
  fprintf(stdout, "      specific surface area : %e\n", constraint->specific_surface_area);
}  //  PrintAlquimiaMineralConstraint()
