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
#include "ini.h"
#include "input_util.h"
#include "alquimia/alquimia_constants.h"
#include "alquimia/alquimia_util.h"
#include "alquimia/alquimia_memory.h"

#define MATCH(s, n) (strcmp(section, s) == 0) && (strcmp(name, n) == 0)

static inline double Min(double a, double b)
{
  return (a < b) ? a : b;
}

static inline double Max(double a, double b)
{
  return (a >= b) ? a : b;
}

// Returns true if the given value string matches a "true" boolean value.
static bool ValueIsTrue(const char* value)
{
  return ((strcmp(value, "1") == 0) || 
          AlquimiaCaseInsensitiveStringCompare(value, "true") ||
          AlquimiaCaseInsensitiveStringCompare(value, "yes") ||
          AlquimiaCaseInsensitiveStringCompare(value, "on")) ? true : false;
}

typedef struct
{
  char engine_name[FILENAME_MAX+1];
  char input_file[FILENAME_MAX+1];
  bool hands_off;
} EngineData;

static int ParseEngineData(void* user,
                           const char* section,
                           const char* name,
                           const char* value)
{
  EngineData* data = user;

  // Look through the chemistry section and get the information we need.
  if (MATCH("chemistry", "engine"))
    strncpy(data->engine_name, value, kAlquimiaMaxStringLength);
  else if (MATCH("chemistry", "input_file"))
    strncpy(data->input_file, value, kAlquimiaMaxStringLength);
  else if (MATCH("chemistry", "hands_off"))
    data->hands_off = ValueIsTrue(value);

  return 1;
}

void Input_CreateAlquimiaInterface(const char* input_file,
                                   AlquimiaInterface* engine_interface,
                                   AlquimiaSizes* engine_sizes,
                                   AlquimiaEngineFunctionality* engine_functionality,
                                   AlquimiaEngineStatus* engine_status)
{
  // Get the engine and other parameters.
  EngineData data;
  int error = ini_parse(input_file, ParseEngineData, &data);
  if (error != 0)
    alquimia_error("Input_CreateAlqumiaInterface: Error parsing input: %s", input_file);

  // Initialize the engine.
  CreateAlquimiaInterface(data.engine_name, engine_interface, engine_status);
  if (engine_status->error != 0) 
    alquimia_error("Input_CreateAlquimiaInterface: %s", engine_status->message);

  // Set it up with our parameters and retrieve functionality.
  engine_interface->Setup(data.input_file,
                          data.hands_off,
                          engine_interface,
                          engine_sizes,
                          engine_functionality,
                          engine_status);

  if (engine_status->error != 0) 
    alquimia_error("Input_CreateAlquimiaInterface: %s", engine_status->message);
}

static int ParseRegionName(void* user,
                           const char* section,
                           const char* name,
                           const char* value)
{
  AlquimiaVectorString* region_names = user;
  if (strstr(section, "region:") != NULL)
  {
    ResizeAlquimiaVectorString(region_names, region_names->size+1);
    strncpy(region_names->data[region_names->size-1], &section[7], kAlquimiaMaxStringLength);
  }
  return 1;
}

void Input_GetRegions(const char* input_file,
                      AlquimiaVectorString* region_names)
{
  int error = ini_parse(input_file, ParseRegionName, region_names);
  if (error != 0)
    alquimia_error("Input_GetRegions: Error parsing input: %s", input_file);
}

static void ParseRegionCells(const char* name,
                             const char* value, 
                             AlquimiaVectorInt* cells)
{
  ResizeAlquimiaVectorInt(cells, 0);

  char* first_colon = strstr(name, ":");
  if (first_colon == NULL)
  {
    // This region is just one cell!
    ResizeAlquimiaVectorInt(cells, 1);
    cells->data[0] = atoi(value);
  }
  else
  {
    char begin[first_colon - name + 1];
    strncpy(begin, name, first_colon - name);
    int b = atoi(begin);

    char* second_colon = strstr(name, ":");
    if (second_colon == NULL)
    {
      // This region is a slice with a stride of 1.
      int e = atoi(&first_colon[1]);
      for (int i = 0; i < e - b; ++i)
      {
        ResizeAlquimiaVectorInt(cells, cells->size+1);
        cells->data[cells->size-1] = b + i;
      }
    }
    else
    {
      // This region is a slice with a specific stride.
      char end[second_colon - first_colon];
      strncpy(end, &first_colon[1], second_colon-first_colon-1);
      int e = atoi(end), stride = atoi(&second_colon[1]);
      for (int i = 0; i < e - b; i += stride)
      {
        ResizeAlquimiaVectorInt(cells, cells->size+1);
        cells->data[cells->size-1] = b + i;
      }
    }
  }
}

typedef struct
{
  char name[FILENAME_MAX+1];
  AlquimiaProblemMetaData* metadata;
  AlquimiaState* state;
  AlquimiaProperties* properties;
  AlquimiaGeochemicalCondition* initial_condition;
  AlquimiaVectorInt* cells;
} RegionData;

static int ParseRegionData(void* user,
                           const char* section,
                           const char* name,
                           const char* value)
{
  RegionData* region_data = user;

  char state_name[kAlquimiaMaxStringLength+1];
  snprintf(state_name, kAlquimiaMaxStringLength, "state:%s", region_data->name);
  char props_name[kAlquimiaMaxStringLength+1];
  snprintf(props_name, kAlquimiaMaxStringLength, "properties:%s", region_data->name);
  char ic_name[kAlquimiaMaxStringLength+1];
  snprintf(ic_name, kAlquimiaMaxStringLength, "initial_cond:%s", region_data->name);
  char region_name[kAlquimiaMaxStringLength+1];
  snprintf(region_name, kAlquimiaMaxStringLength, "region:%s", region_data->name);

  if (strcmp(section, state_name) == 0)
    Input_ParseState(section, name, value, region_data->metadata, region_data->state);
  else if (strcmp(section, props_name) == 0)
    Input_ParseProperty(section, name, value, region_data->metadata, region_data->properties);
  else if (strcmp(section, ic_name) == 0)
    Input_ParseGeochemicalCondition(name, value, region_data->initial_condition);
  else if (strcmp(section, region_name) == 0)
    ParseRegionCells(name, value, region_data->cells);
  return 1;
}

void Input_GetRegionData(const char* input_file,
                         const char* region_name,
                         AlquimiaProblemMetaData* problem_metadata,
                         AlquimiaState* region_state,
                         AlquimiaProperties* region_properties,
                         AlquimiaGeochemicalCondition* region_initial_condition,
                         AlquimiaVectorInt* region_cells)
{
  RegionData region_data = {.metadata = problem_metadata,
                            .state = region_state,
                            .properties = region_properties,
                            .initial_condition = region_initial_condition,
                            .cells = region_cells};
  strncpy(region_data.name, region_name, kAlquimiaMaxStringLength);
  int error = ini_parse(input_file, ParseRegionData, &region_data);
  if (error != 0)
    alquimia_error("Input_GetRegionData: Error parsing input: %s", input_file);
}

static int ParseGeochemicalConditions(void* user,
                                      const char* section,
                                      const char* name,
                                      const char* value)
{
  AlquimiaGeochemicalConditionVector* conditions = user;

  char cond_name[kAlquimiaMaxStringLength+1];
  if (Input_IsGeochemicalConditionSection(section, cond_name))
  {
    int icond = 0;
    for (; icond < conditions->size; ++icond)
    {
      if (strcmp(conditions->data[icond].name, cond_name) == 0)
        break;
    }
    if (icond == conditions->size)
    {
      ResizeAlquimiaGeochemicalConditionVector(conditions, conditions->size+1);
      strncpy(conditions->data[icond].name, cond_name, kAlquimiaMaxStringLength);
    }
    AlquimiaGeochemicalCondition* condition = &conditions->data[icond];
    Input_ParseGeochemicalCondition(name, value, condition);
  }

  return 1;
}

void Input_GetGeochemicalConditions(const char* input_file,
                                    AlquimiaGeochemicalConditionVector* conditions)
{
  int error = ini_parse(input_file, ParseGeochemicalConditions, conditions);
  if (error != 0)
    alquimia_error("Input_GetGeochemicalConditions: Error parsing input: %s", input_file);
}

typedef struct
{
  char type[FILENAME_MAX+1];
  char file[FILENAME_MAX+1];
  bool verbose;
} OutputData;

static int ParseOutputParameters(void* user,
                                 const char* section,
                                 const char* name,
                                 const char* value)
{
  OutputData* output = user;
  if (MATCH("output", "verbose"))
    output->verbose = ValueIsTrue(value);
  else if (MATCH("output", "type"))
    strncpy(output->type, value, kAlquimiaMaxStringLength);
  else if (MATCH("output", "filename"))
    strncpy(output->file, value, kAlquimiaMaxStringLength);

  return 1;
}

void Input_GetOutputParameters(const char* input_file,
                               char* output_type,
                               char* output_file,
                               bool* verbose)
{
  OutputData output;
  output.type[0] = '\0';
  output.file[0] = '\0';
  int error = ini_parse(input_file, ParseOutputParameters, &output);
  if (error != 0)
    alquimia_error("Input_GetOutputParameters: Error parsing input: %s", input_file);

  if (strlen(output.type) == 0)
    strcpy(output_type, "gnuplot"); // default output type
  else
    strncpy(output_type, output.type, kAlquimiaMaxStringLength);

  char default_output_file[kAlquimiaMaxStringLength+1];
  if (strlen(output.file) == 0)
  {
    // Find the last '.' in the input filename.
    int dot = strlen(input_file)-1;
    while ((dot > 0) && (input_file[dot] != '.')) --dot;
    char suffix[16];
    suffix[0] = '\0';

    // Determine the suffix from the output type.
    if (AlquimiaCaseInsensitiveStringCompare(output.type, "gnuplot"))
      sprintf(suffix, ".gnuplot");
    else if (AlquimiaCaseInsensitiveStringCompare(output.type, "python"))
      sprintf(suffix, ".py");

    // Append the suffix.
    if (dot == 0)
      sprintf(default_output_file, "%s%s", input_file, suffix);
    else
    {
      memcpy(default_output_file, input_file, sizeof(char) * dot);
      default_output_file[dot] = '\0';
      strcat(default_output_file, suffix);
    }

    // Python modules can't have hyphens in their names, so we replace them 
    // with underscores.
    if (AlquimiaCaseInsensitiveStringCompare(output.type, "python"))
    {
      int len = strlen(default_output_file);
      for (int i = 0; i < len; ++i)
      {
        if (default_output_file[i] == '-')
          default_output_file[i] = '_';
      }
    }

    strncpy(output_file, default_output_file, kAlquimiaMaxStringLength);
  }
  else
    strncpy(output_file, output.file, kAlquimiaMaxStringLength);
}

bool Input_IsStateSection(const char* section, char* state_name)
{
  bool result = (strstr(section, "state") != NULL);
  if (result && strstr(section, "state:") != NULL)
    strcpy(state_name, &section[6]);
  else
    state_name[0] = '\0';
  return result;
}

void Input_ParseState(const char* section,
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

bool Input_IsPropertiesSection(const char* section, char* properties_name)
{
  bool result = (strstr(section, "properties") != NULL);
  if (result && strstr(section, "properties:") != NULL)
    strcpy(properties_name, &section[11]);
  else
    properties_name[0] = '\0';
  return result;
}

void Input_ParseProperty(const char* section,
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

bool Input_IsGeochemicalConditionSection(const char* section, char* condition_name)
{
  bool result = (strstr(section, "condition:") != NULL);
  if (result)
    strcpy(condition_name, &section[10]);
  return result;
}

void Input_ParseGeochemicalCondition(const char* name,
                                     const char* value,
                                     AlquimiaGeochemicalCondition* condition)
{
  if (strstr(value, "aqueous_constraint(") != NULL)
  {
    ResizeAlquimiaAqueousConstraintVector(&condition->aqueous_constraints, condition->aqueous_constraints.size+1);
    Input_ParseAqueousConstraint(name, value, &condition->aqueous_constraints.data[condition->aqueous_constraints.size-1]);
  }
  else if (strstr(value, "mineral_constraint(") != NULL)
  {
    ResizeAlquimiaMineralConstraintVector(&condition->mineral_constraints, condition->mineral_constraints.size+1);
    Input_ParseMineralConstraint(name, value, &condition->mineral_constraints.data[condition->mineral_constraints.size-1]);
  }
}

void Input_ParseAqueousConstraint(const char* primary_species,
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

void Input_ParseMineralConstraint(const char* mineral_name,
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

