/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
** Alquimia Copyright (c) 2013-2015, The Regents of the University of California, 
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
*/

#include <string.h>
#include "input_util.h"
#include "alquimia/alquimia_util.h"

bool IsStateSection(const char* section, char* state_name)
{
  bool result = (strstr(section, "state") != NULL);
  if (result && strstr(section, "state:") != NULL)
    strcpy(state_name, &section[6]);
  else
    state_name[0] = '\0';
  return result;
}

void ParseStateInput(const char* section,
                     const char* name,
                     const char* value,
                     AlquimiaProblemMetaData* metadata,
                     AlquimiaState* state)
{
  if (strcmp(name, "density") == 0)
  {
    state->water_density = atof(value);
    if (state->water_density <= 0)
      alquimia_error("In section %s: Non-positive density: %g", state->water_density);
  }
  else if (strcmp(name, "porosity") == 0)
  {
    state->porosity = atof(value);
    if ((state->porosity <= 0) || (state->porosity > 1.0))
      alquimia_error("In section %s: Invalid porosity: %g", state->porosity);
  }
  else if (strcmp(name, "temperature") == 0)
  {
    state->temperature = atof(value);
    if (state->temperature <= -273.0) 
      alquimia_error("In section %s: Invalid temperature: %g", state->temperature);
  }
  else if (strcmp(name, "pressure") == 0)
    state->aqueous_pressure = atof(value);
}

bool IsPropertiesSection(const char* section, char* properties_name)
{
  bool result = (strstr(section, "properties") != NULL);
  if (result && strstr(section, "properties:") != NULL)
    strcpy(properties_name, &section[11]);
  else
    properties_name[0] = '\0';
  return result;
}

void ParsePropertyInput(const char* section,
                        const char* name,
                        const char* value,
                        AlquimiaProblemMetaData* metadata,
                        AlquimiaProperties* properties)
{
  if (strcmp(name, "volume") == 0)
  {
    properties->volume = atof(value);
    if (properties->volume <= 0)
      alquimia_error("In section %s: Non-positive volume: %g", properties->volume);
  }
  else if (strcmp(name, "saturation") == 0)
  {
    properties->saturation = atof(value);
    if ((properties->saturation < 0) || (properties->saturation > 1.0))
      alquimia_error("In section %s: Invalid saturation: %g", properties->saturation);
  }
  else if (strstr(name, "[") != NULL)
  {
    char* begin = (char*)&name[11];
    char* end = strstr(begin, "]");
    char id[end - begin + 1];
    strncpy(id, begin, end-begin);
    if (strstr(name, "isotherm_kd[") != NULL)
    {
      int index;
      AlquimiaFindIndexFromName(id, &metadata->isotherm_species_names, &index);
      if (index == -1)
        alquimia_error("In section %s: Invalid isotherm species identifier: %s", section, name);
      properties->isotherm_kd.data[index] = atof(value);
    }
    else if (strstr(name, "freundlich_n[") != NULL)
    {
      int index;
      AlquimiaFindIndexFromName(id, &metadata->isotherm_species_names, &index);
      if (index == -1)
        alquimia_error("In section %s: Invalid isotherm species identifier: %s", section, name);
      properties->freundlich_n.data[index] = atof(value);
    }
    else if (strstr(name, "langmuir_b[") != NULL)
    {
      int index;
      AlquimiaFindIndexFromName(id, &metadata->isotherm_species_names, &index);
      if (index == -1)
        alquimia_error("In section %s: Invalid isotherm species identifier: %s", section, name);
      properties->langmuir_b.data[index] = atof(value);
    }
    else if (strstr(name, "mineral_rate_cnst[") != NULL)
    {
      int index;
      AlquimiaFindIndexFromName(id, &metadata->mineral_names, &index);
      if (index == -1)
        alquimia_error("In section %s: Invalid mineral identifier: %s", section, name);
      properties->mineral_rate_cnst.data[index] = atof(value);
    }
    else if (strstr(name, "aqueous_kinetic_rate_cnst[") != NULL)
    {
      int index;
      AlquimiaFindIndexFromName(id, &metadata->aqueous_kinetic_names, &index);
      if (index == -1)
        alquimia_error("In section %s: Invalid aqueous kinetic identifier: %s", section, name);
      properties->aqueous_kinetic_rate_cnst.data[index] = atof(value);
    }
  }
}

bool IsGeochemicalConditionSection(const char* section, char* condition_name)
{
  bool result = (strstr(section, "condition:") != NULL);
  if (result)
    strcpy(condition_name, &section[10]);
  return result;
}

void ParseGeochemicalConditionInput(const char* name,
                                    const char* value,
                                    AlquimiaGeochemicalCondition* condition)
{
  if (strstr(value, "aqueous_constraint(") != NULL)
  {
    ResizeAlquimiaAqueousConstraintVector(&condition->aqueous_constraints, condition->aqueous_constraints.size+1);
    ParseAqueousConstraintInput(name, value, &condition->aqueous_constraints.data[condition->aqueous_constraints.size-1]);
  }
  else if (strstr(value, "mineral_constraint(") != NULL)
  {
    ResizeAlquimiaMineralConstraintVector(&condition->mineral_constraints, condition->mineral_constraints.size+1);
    ParseMineralConstraintInput(name, value, &condition->mineral_constraints.data[condition->mineral_constraints.size-1]);
  }
}

void ParseAqueousConstraintInput(const char* primary_species,
                                 const char* text,
                                 AlquimiaAqueousConstraint* constraint)
{
  char* start = strstr(text, "aqueous_constraint(");
  if (start == NULL) return;

  char* comma = strstr(&start[19], ",");
  if (comma == NULL)
    alquimia_error("Invalid aqueous_constraint entry for %s: %s", primary_species, text);

  char* last_comma = comma;
  char* comma2 = strstr(comma+1, ",");
  if (comma2 != NULL)
    last_comma = comma2;
  char* end = strstr(last_comma, ")");
  if (end == NULL)
    alquimia_error("Invalid aqueous_constraint entry for %s: %s", primary_species, text);

  // Primary species.
  strcpy(constraint->primary_species_name, primary_species);

  // Jot down the value.
  char value_str[comma - &start[19]+1];
  strncpy(value_str, &start[19], comma-&start[19]);
  constraint->value = atof(value_str);

  // Figure out the type of constraint.
  char* end_type_delim = (last_comma == comma2) ? comma2 : end;
  strncpy(constraint->constraint_type, comma+1, end_type_delim-(comma+1));

  // Is there an associated species?
  if (last_comma == comma2)
    strncpy(constraint->associated_species, comma2+1, end-(comma2+1));
}

void ParseMineralConstraintInput(const char* mineral_name,
                                 const char* text,
                                 AlquimiaMineralConstraint* constraint)
{
  char* start = strstr(text, "mineral_constraint(");
  if (start == NULL) return;

  char* comma = strstr(&start[19], ",");
  if (comma == NULL)
    alquimia_error("Invalid mineral_constraint entry for %s: %s", mineral_name, text);

  char* end = strstr(comma, ")");
  if (end == NULL)
    alquimia_error("Invalid mineral_constraint entry for %s: %s", mineral_name, text);

  // Mineral name.
  strcpy(constraint->mineral_name, mineral_name);

  // Volume fraction.
  char vf_str[comma - &start[19]+1];
  strncpy(vf_str, &start[19], comma-&start[19]);
  constraint->volume_fraction = atof(vf_str);

  // Specific surface area.
  char ssa_str[end - comma + 1];
  strncpy(ssa_str, comma+1, end-comma-1);
  constraint->specific_surface_area = atof(ssa_str);
}

