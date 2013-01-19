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
  fprintf(stdout, "    num kinetic minerals : %d\n", sizes->num_kinetic_minerals);
  fprintf(stdout, "    num aqueous complexes : %d\n", sizes->num_aqueous_complexes);
  fprintf(stdout, "    num surface sites : %d\n", sizes->num_surface_sites);
  fprintf(stdout, "    num ion exchange sites : %d\n", sizes->num_ion_exchange_sites);
}  // end PrintAlquimiaSizes()

void PrintAlquimiaMetaData(const struct AlquimiaSizes* sizes,
                           const struct AlquimiaMetaData* metadata) {
  int i;
  fprintf(stdout, "  Alquimia Metadata :\n");
  fprintf(stdout, "    thread_safe : %d\n", metadata->thread_safe);
  fprintf(stdout, "    temperature_dependent : %d\n",
          metadata->temperature_dependent);
  fprintf(stdout, "    pressure_dependent : %d\n", metadata->pressure_dependent);
  fprintf(stdout, "    porosity_update  : %d\n", metadata->porosity_update);
  fprintf(stdout, "    index base : %d\n", metadata->index_base);
  fprintf(stdout, "    primary indices (%d) (%p) :\n      ",
          sizes->num_primary, metadata->primary_indices);
  for (i = 0; i < sizes->num_primary; ++i) {
    fprintf(stdout, "%d, ", metadata->primary_indices[i]);
  }
  fprintf(stdout, "\n");
  fprintf(stdout, "    primary names (%d) (%p) :\n      ",
          sizes->num_primary, metadata->primary_names);
  for (i = 0; i < sizes->num_primary; ++i) {
    fprintf(stdout, "'%s', ", metadata->primary_names[i]);
  }
  fprintf(stdout, "\n");
}  // end PrintAlquimiaMetaData()

void PrintAlquimiaState(const struct AlquimiaSizes* sizes,
                        const struct AlquimiaState* state) {
  int i;
  fprintf(stdout, "  Alquimia State:\n");
  fprintf(stdout, "    water density : %f\n", state->water_density);
  fprintf(stdout, "    saturation : %f\n", state->saturation);
  fprintf(stdout, "    porosity : %f\n", state->porosity);
  fprintf(stdout, "    temperature : %f\n", state->temperature);
  fprintf(stdout, "    aqueous_pressure : %f\n", state->aqueous_pressure);

  fprintf(stdout, "    total_primary (%d) :\n    ", sizes->num_primary);
  for (i = 0; i < sizes->num_primary; ++i) {
    fprintf(stdout, "%e, ", state->total_primary[i]);
  }
  fprintf(stdout, "\n");

  fprintf(stdout, "    free_ion (%d) :\n    ", sizes->num_primary);
  for (i = 0; i < sizes->num_primary; ++i) {
    fprintf(stdout, "%e, ", state->free_ion[i]);
  }
  fprintf(stdout, "\n");

  fprintf(stdout, "    kinetic minerals volume fraction (%d) :\n    ",
          sizes->num_kinetic_minerals);
  for (i = 0; i < sizes->num_kinetic_minerals; ++i) {
    fprintf(stdout, "%e, ", state->mineral_volume_fraction[i]);
  }
  fprintf(stdout, "\n");

  fprintf(stdout, "    kinetic minerals specific surface area (%d) :\n    ",
          sizes->num_kinetic_minerals);
  for (i = 0; i < sizes->num_kinetic_minerals; ++i) {
    fprintf(stdout, "%e, ", state->mineral_specific_surface_area[i]);
  }
  fprintf(stdout, "\n");
}  // end PrintAlquimiaState()

void PrintAlquimiaAuxiliaryData(const struct AlquimiaSizes* sizes,
                        const struct AlquimiaAuxiliaryData* aux_data) {
  int i;
  fprintf(stdout, "  Alquimia Auxiliary Data:\n");
  fprintf(stdout, "    primary activity coeff (%d) :\n    [ ", sizes->num_primary);
  for (i = 0; i < sizes->num_primary; ++i) {
    fprintf(stdout, "%e, ", aux_data->primary_activity_coeff[i]);
  }
  fprintf(stdout, "]\n");

  fprintf(stdout, "    secondary activity coeff (%d) :\n    [ ",
          sizes->num_aqueous_complexes);
  for (i = 0; i < sizes->num_aqueous_complexes; ++i) {
    fprintf(stdout, "%e, ", aux_data->secondary_activity_coeff[i]);
  }
  fprintf(stdout, "]\n");

  fprintf(stdout, "    ion exchange ref cation conc (%d) :\n    [ ",
          sizes->num_ion_exchange_sites);
  for (i = 0; i < sizes->num_ion_exchange_sites; ++i) {
    fprintf(stdout, "%e, ", aux_data->ion_exchange_ref_cation_conc[i]);
  }
  fprintf(stdout, "]\n");

  fprintf(stdout, "    surface complex free site conc (%d) :\n    [ ",
          sizes->num_surface_sites);
  for (i = 0; i < sizes->num_surface_sites; ++i) {
    fprintf(stdout, "%e, ", aux_data->surface_complex_free_site_conc[i]);
  }
  fprintf(stdout, "]\n");
}  // end PrintAlquimiaAuxiliaryData()

void PrintAlquimiaGeochemicalConditionList(
    const struct AlquimiaGeochemicalConditionList* condition_list) {
  int i;
  fprintf(stdout, "Alquimia Geochemical Condition List : \n");
  for (i = 0; i < condition_list->num_conditions; ++i) {
    PrintAlquimiaGeochemicalCondition(&(condition_list->conditions[i]));
    fprintf(stdout, "\n");
  }
}  //  PrintAlquimiaGeochemicalConditionList()

void PrintAlquimiaGeochemicalCondition(
    const struct AlquimiaGeochemicalCondition* condition) {
  int i;
  fprintf(stdout, "  Alquimia Geochemical Condition : %s\n", condition->name);
  for (i = 0; i < condition->num_constraints; ++i) {
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
