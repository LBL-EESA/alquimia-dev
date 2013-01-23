/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*******************************************************************************
 **
 **  C utilities for working with alquimia data structures
 **
 *******************************************************************************/

#include "alquimia_util.h"

#include <stdlib.h>
#include <stdio.h>

#include "alquimia_containers.h"

/*******************************************************************************
 **
 **  Printing
 **
 *******************************************************************************/

void PrintAlquimiaSizes(const struct AlquimiaSizes* sizes) {
  fprintf(stdout, "  Alquimia Sizes :\n");
  fprintf(stdout, "    num primary species : %d\n", sizes->num_primary);
  fprintf(stdout, "    num sorbed : %d\n", sizes->num_sorbed);
  fprintf(stdout, "    num kinetic minerals : %d\n", sizes->num_kinetic_minerals);
  fprintf(stdout, "    num aqueous complexes : %d\n", sizes->num_aqueous_complexes);
  fprintf(stdout, "    num surface sites : %d\n", sizes->num_surface_sites);
  fprintf(stdout, "    num ion exchange sites : %d\n", sizes->num_ion_exchange_sites);
}  // end PrintAlquimiaSizes()

void PrintAlquimiaMetaData(const struct AlquimiaMetaData* metadata) {

  fprintf(stdout, "  Alquimia Metadata :\n");
  fprintf(stdout, "    thread_safe : %d\n", metadata->thread_safe);
  fprintf(stdout, "    temperature_dependent : %d\n",
          metadata->temperature_dependent);
  fprintf(stdout, "    pressure_dependent : %d\n", metadata->pressure_dependent);
  fprintf(stdout, "    porosity_update  : %d\n", metadata->porosity_update);
  fprintf(stdout, "    index base : %d\n", metadata->index_base);
  fprintf(stdout, "    primary indices (%d) (%p) :\n      ",
          metadata->size_primary, metadata->primary_indices);
  for (int i = 0; i < metadata->size_primary; ++i) {
    fprintf(stdout, "%d, ", metadata->primary_indices[i]);
  }
  fprintf(stdout, "\n");
  fprintf(stdout, "    primary names (%d) (%p) :\n      ",
          metadata->size_primary, metadata->primary_names);
  for (int i = 0; i < metadata->size_primary; ++i) {
    fprintf(stdout, "'%s', ", metadata->primary_names[i]);
  }
  fprintf(stdout, "\n");
  fprintf(stdout, "    mineral indices (%d) (%p) :\n      ",
          metadata->size_minerals, metadata->mineral_indices);
  for (int i = 0; i < metadata->size_minerals; ++i) {
    fprintf(stdout, "%d, ", metadata->mineral_indices[i]);
  }
  fprintf(stdout, "\n");
  fprintf(stdout, "    mineral names (%d) (%p) :\n      ",
          metadata->size_minerals, metadata->mineral_names);
  for (int i = 0; i < metadata->size_minerals; ++i) {
    fprintf(stdout, "'%s', ", metadata->mineral_names[i]);
  }
  fprintf(stdout, "\n");
}  // end PrintAlquimiaMetaData()

void PrintAlquimiaState(const struct AlquimiaState* state) {

  fprintf(stdout, "  Alquimia State:\n");
  fprintf(stdout, "    water density : %f\n", state->water_density);
  fprintf(stdout, "    saturation : %f\n", state->saturation);
  fprintf(stdout, "    porosity : %f\n", state->porosity);
  fprintf(stdout, "    temperature : %f\n", state->temperature);
  fprintf(stdout, "    aqueous_pressure : %f\n", state->aqueous_pressure);

  fprintf(stdout, "    total_primary (%d) :\n    ", state->size_total_primary);
  for (int i = 0; i < state->size_total_primary; ++i) {
    fprintf(stdout, "%e, ", state->total_primary[i]);
  }
  fprintf(stdout, "\n");

  fprintf(stdout, "    free_ion (%d) :\n    ", state->size_free_ion);
  for (int i = 0; i < state->size_free_ion; ++i) {
    fprintf(stdout, "%e, ", state->free_ion[i]);
  }
  fprintf(stdout, "\n");

  fprintf(stdout, "    kinetic minerals volume fraction (%d) :\n    ",
          state->size_mineral_volume_fraction);
  for (int i = 0; i < state->size_mineral_volume_fraction; ++i) {
    fprintf(stdout, "%e, ", state->mineral_volume_fraction[i]);
  }
  fprintf(stdout, "\n");

  fprintf(stdout, "    kinetic minerals specific surface area (%d) :\n    ",
          state->size_mineral_specific_surface_area);
  for (int i = 0; i < state->size_mineral_specific_surface_area; ++i) {
    fprintf(stdout, "%e, ", state->mineral_specific_surface_area[i]);
  }
  fprintf(stdout, "\n");
}  // end PrintAlquimiaState()

void PrintAlquimiaAuxiliaryData(const struct AlquimiaAuxiliaryData* aux_data) {

  fprintf(stdout, "  Alquimia Auxiliary Data:\n");
  fprintf(stdout, "    primary activity coeff (%d) :\n    [ ",
          aux_data->size_primary_activity_coeff);
  for (int i = 0; i < aux_data->size_primary_activity_coeff; ++i) {
    fprintf(stdout, "%e, ", aux_data->primary_activity_coeff[i]);
  }
  fprintf(stdout, "]\n");

  fprintf(stdout, "    secondary activity coeff (%d) :\n    [ ",
          aux_data->size_secondary_activity_coeff);
  for (int i = 0; i < aux_data->size_secondary_activity_coeff; ++i) {
    fprintf(stdout, "%e, ", aux_data->secondary_activity_coeff[i]);
  }
  fprintf(stdout, "]\n");

  fprintf(stdout, "    ion exchange ref cation conc (%d) :\n    [ ",
          aux_data->size_ion_exchange_ref_cation_conc);
  for (int i = 0; i < aux_data->size_ion_exchange_ref_cation_conc; ++i) {
    fprintf(stdout, "%e, ", aux_data->ion_exchange_ref_cation_conc[i]);
  }
  fprintf(stdout, "]\n");

  fprintf(stdout, "    surface complex free site conc (%d) :\n    [ ",
          aux_data->size_surface_complex_free_site_conc);
  for (int i = 0; i < aux_data->size_surface_complex_free_site_conc; ++i) {
    fprintf(stdout, "%e, ", aux_data->surface_complex_free_site_conc[i]);
  }
  fprintf(stdout, "]\n");
}  // end PrintAlquimiaAuxiliaryData()

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
