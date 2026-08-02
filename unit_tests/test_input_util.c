/* -*-  mode: c; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
** Alquimia Copyright (c) 2013-2016, The Regents of the University of
** California, through Lawrence Berkeley National Laboratory.
**
** Alquimia is available under a BSD license. See LICENSE.txt for more
** information.
*/

#include "alquimia/alquimia_constants.h"
#include "alquimia/alquimia_memory.h"
#include "drivers/input_util.h"

static void TestGeochemicalConditionSections(void)
{
  char condition_name[kAlquimiaMaxStringLength+1];
  ALQUIMIA_ASSERT(Input_IsGeochemicalConditionSection(
      "condition:first", condition_name));
  ALQUIMIA_ASSERT(strcmp(condition_name, "first") == 0);
  ALQUIMIA_ASSERT(!Input_IsGeochemicalConditionSection(
      "notcondition:first", condition_name));
  ALQUIMIA_ASSERT(!Input_IsGeochemicalConditionSection(
      "condition:", condition_name));
}

static void TestGeochemicalConditionParsing(void)
{
  char input_file[FILENAME_MAX];
  snprintf(input_file, FILENAME_MAX, "%s/test_input_util.cfg",
           CMAKE_CURRENT_SOURCE_DIR);

  AlquimiaGeochemicalConditionVector conditions = {
      .size = 0, .capacity = 0, .data = NULL};
  Input_GetGeochemicalConditions(input_file, &conditions);

  ALQUIMIA_ASSERT(conditions.size == 2);
  const AlquimiaGeochemicalCondition* first =
      Input_FindGeochemicalCondition(&conditions, "first");
  const AlquimiaGeochemicalCondition* second =
      Input_FindGeochemicalCondition(&conditions, "second");
  ALQUIMIA_ASSERT(first != NULL);
  ALQUIMIA_ASSERT(second != NULL);
  ALQUIMIA_ASSERT(Input_FindGeochemicalCondition(
      &conditions, "missing") == NULL);

  ALQUIMIA_ASSERT(first->aqueous_constraints.size == 3);
  ALQUIMIA_ASSERT(first->mineral_constraints.size == 1);
  ALQUIMIA_ASSERT(strcmp(
      first->aqueous_constraints.data[0].primary_species_name, "H+") == 0);
  ALQUIMIA_ASSERT(strcmp(
      first->aqueous_constraints.data[0].constraint_type, "pH") == 0);
  ALQUIMIA_ASSERT(first->aqueous_constraints.data[0].value == 6.0);
  ALQUIMIA_ASSERT(strcmp(
      first->aqueous_constraints.data[1].constraint_type, "mineral") == 0);
  ALQUIMIA_ASSERT(strcmp(
      first->aqueous_constraints.data[1].associated_species, "Halite") == 0);
  ALQUIMIA_ASSERT(strcmp(
      first->mineral_constraints.data[0].mineral_name, "Halite") == 0);
  ALQUIMIA_ASSERT(
      first->mineral_constraints.data[0].volume_fraction == 1.0e-5);
  ALQUIMIA_ASSERT(
      first->mineral_constraints.data[0].specific_surface_area == 1.0);

  ALQUIMIA_ASSERT(second->aqueous_constraints.size == 2);
  ALQUIMIA_ASSERT(second->mineral_constraints.size == 0);
  ALQUIMIA_ASSERT(strcmp(
      second->aqueous_constraints.data[0].constraint_type,
      "total_aqueous") == 0);
  ALQUIMIA_ASSERT(strcmp(
      second->aqueous_constraints.data[1].constraint_type, "free") == 0);

  FreeAlquimiaGeochemicalConditionVector(&conditions);
  ALQUIMIA_ASSERT(conditions.size == 0);
  ALQUIMIA_ASSERT(conditions.capacity == 0);
  ALQUIMIA_ASSERT(conditions.data == NULL);
}

int main(int argc, char** argv)
{
  (void)argc;
  (void)argv;

  printf("Testing input utilities.\n");
  TestGeochemicalConditionSections();
  TestGeochemicalConditionParsing();
  printf("All tests passed.\n");
  return EXIT_SUCCESS;
}
