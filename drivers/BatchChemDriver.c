/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

//
// Alquimia Copyright (c) 2013-2015, The Regents of the University of California, 
// through Lawrence Berkeley National Laboratory (subject to receipt of any 
// required approvals from the U.S. Dept. of Energy).  All rights reserved.
// 
// Alquimia is available under a BSD license. See LICENSE.txt for more
// information.
//
// If you have questions about your rights to use or distribute this software, 
// please contact Berkeley Lab's Technology Transfer and Intellectual Property 
// Management at TTD@lbl.gov referring to Alquimia (LBNL Ref. 2013-119).
// 
// NOTICE.  This software was developed under funding from the U.S. Department 
// of Energy.  As such, the U.S. Government has been granted for itself and 
// others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide 
// license in the Software to reproduce, prepare derivative works, and perform 
// publicly and display publicly.  Beginning five (5) years after the date 
// permission to assert copyright is obtained from the U.S. Department of Energy, 
// and subject to any subsequent five (5) year renewals, the U.S. Government is 
// granted for itself and others acting on its behalf a paid-up, nonexclusive, 
// irrevocable, worldwide license in the Software to reproduce, prepare derivative
// works, distribute copies to the public, perform publicly and display publicly, 
// and to permit others to do so.
//

#include <limits.h>
#include <float.h>
#include "alquimia/alquimia_interface.h"
#include "alquimia/alquimia_memory.h"
#include "alquimia/alquimia_util.h"
#include "ini.h"
#include "BatchChemDriver.h"

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

static void GetIndexLabel(const char* name, char* label)
{
  size_t br1 = strstr(name, "[") - name;
  size_t br2 = strstr(name, "]") - name;
  memcpy(label, &name[br1], sizeof(char) * (br2-br1));
}

// The following macros are for parsing values that are indexed by names.
// Avert your eyes.
#define DECLARE_INDEXED_VALUE_PARSER(IndexName, name_array, num_index_var) \
static int IndexName##Index(BatchChemDriverInput* input, const char* item) \
{ \
  int i = 0; \
  for (; i < input->num_index_var; ++i) \
  { \
    if (strncmp(input->name_array[i], item, 127) == 0) \
      break; \
  } \
  if (i == input->num_index_var) \
  { \
    input->name_array[i] = AlquimiaStringDup(item); \
    ++input->num_index_var; \
  } \
  return i; \
} \
\
static void Get##IndexName##Value(BatchChemDriverInput* input,\
                                  const char* name, \
                                  const char* value, \
                                  double* array) \
{ \
  char item[128]; \
  GetIndexLabel(name, item); \
  int index = IndexName##Index(input, item); \
  array[index] = atof(value); \
}

DECLARE_INDEXED_VALUE_PARSER(Isotherm, isotherm_species, num_isotherm_species)
DECLARE_INDEXED_VALUE_PARSER(IonExchange, ion_exchange_sites, num_ion_exchange_sites)
DECLARE_INDEXED_VALUE_PARSER(SurfaceSite, surface_sites, num_surface_sites)

// Input parsing stuff. See https://github.com/benhoyt/inih for details.
static int ParseInput(void* user,
                      const char* section,
                      const char* name,
                      const char* value)
{
  BatchChemDriverInput* input = user;
#define MATCH(s, x) (strcmp(section, s) == 0) && (strcmp(name, x) == 0)
#define PARTIAL_MATCH(s, x, n) (strstr(section, s) != NULL) && (strstr(name, x, n) != NULL)
#define MATCH_N(s, x, n) (strcmp(section, s) == 0) && (strncmp(name, x, n) == 0)

  // Simulation section
  if (MATCH("simulation","description"))
    input->description = AlquimiaStringDup(value);
  else if (MATCH("simulation", "hands_off"))
    input->hands_off = ValueIsTrue(value);
  else if (MATCH("simulation","t_min"))
    input->t_min = atof(value);
  else if (MATCH("simulation","t_max"))
    input->t_max = atof(value);
  else if (MATCH("simulation","max_steps"))
    input->max_steps = atoi(value);
  else if (MATCH("simulation","timestep"))
    input->dt = atof(value);

  // Material section.
  else if (MATCH("material", "volume"))
    input->volume = atof(value);
  else if (MATCH("material", "saturation"))
    input->saturation = atof(value);
  else if (MATCH_N("material", "isotherm_kd[", 12))
    GetIsothermValue(input, name, value, input->isotherm_kd);
  else if (MATCH_N("material", "langmuir_b[", 11))
    GetIsothermValue(input, name, value, input->langmuir_b);
  else if (MATCH_N("material", "freundlich_n[", 13))
    GetIsothermValue(input, name, value, input->freundlich_n);

  // State section.
  else if (MATCH("state", "density"))
    input->water_density = atof(value);
  else if (MATCH("state", "porosity"))
    input->porosity = atof(value);
  else if (MATCH("state", "temperature"))
    input->temperature = atof(value);
  else if (MATCH("state", "pressure"))
    input->aqueous_pressure = atof(value);
  else if (MATCH_N("state", "surface_site_density[", 21))
    GetSurfaceSiteValue(input, name, value, input->surface_site_density);
  else if (MATCH_N("state", "cation_exchange_capacity[", 25))
    GetIonExchangeValue(input, name, value, input->cation_exchange_capacity);

  // Chemistry section.
  else if (MATCH("chemistry", "engine"))
    input->chemistry_engine = AlquimiaStringDup(value);
  else if (MATCH("chemistry", "input_file"))
    input->chemistry_input_file = AlquimiaStringDup(value);
  else if (MATCH("chemistry", "initial_condition"))
    input->cond_name = AlquimiaStringDup(value);

  // Output section.
  else if (MATCH("output", "verbose"))
    input->verbose = ValueIsTrue(value);
  else if (MATCH("output", "type"))
    input->output_type = AlquimiaStringDup(value);
  else if (MATCH("output", "filename"))
    input->output_file = AlquimiaStringDup(value);

  return 1;
}

BatchChemDriverInput* BatchChemDriverInput_New(const char* input_file)
{
  BatchChemDriverInput* input = malloc(sizeof(BatchChemDriverInput));
  memset(input, 0, sizeof(BatchChemDriverInput));

  // Make sure we have some meaningful defaults.
  input->hands_off = true; // Hands-off by default.
  input->t_min = 0.0;
  input->t_max = FLT_MAX;
  input->max_steps = INT_MAX;
  input->dt = FLT_MAX;
  input->num_isotherm_species = 0;
  memset(input->isotherm_kd, 0, sizeof(double) * BATCH_CHEM_INPUT_MAX);
  memset(input->langmuir_b, 0, sizeof(double) * BATCH_CHEM_INPUT_MAX);
  memset(input->freundlich_n, 0, sizeof(double) * BATCH_CHEM_INPUT_MAX);
  input->num_ion_exchange_sites = 0;
  memset(input->cation_exchange_capacity, 0, sizeof(double) * BATCH_CHEM_INPUT_MAX);
  input->num_surface_sites = 0;
  memset(input->surface_site_density, 0, sizeof(double) * BATCH_CHEM_INPUT_MAX);
  input->volume = 1.0;
  input->saturation = 1.0;
  input->water_density = 999.9720;    // density of water in kg/m**3
  input->temperature = 25.0;
  input->porosity = 1.0;
  input->aqueous_pressure = 201325.0; // pressure in Pa.
  input->output_type = NULL;
  input->output_file = NULL;

  // Fill in fields by parsing the input file.
  int error = ini_parse(input_file, ParseInput, input);
  if (error != 0)
  {
    BatchChemDriverInput_Free(input);
    alquimia_error("BatchChemDriver: Error parsing input: %s", input_file);
  }

  // Verify that our required fields are filled properly.
  if (!input->hands_off)
    alquimia_error("BatchChemDriver: simulation->hands_off must be set to true at the moment.");
  if (input->t_max <= input->t_min)
    alquimia_error("BatchChemDriver: simulation->t_max must be greater than simulation->t_min.");
  if (input->max_steps < 0)
    alquimia_error("BatchChemDriver: simulation->max_steps must be non-negative.");
  if ((input->porosity <= 0.0) || (input->porosity > 1.0))
    alquimia_error("BatchChemDriver: material->porosity must be within (0, 1].");
  if ((input->saturation <= 0.0) || (input->saturation > 1.0))
    alquimia_error("BatchChemDriver: material->saturation must be within (0, 1].");
  if (input->water_density <= 0.0)
    alquimia_error("BatchChemDriver: flow->water_density must be positive.");
  if (input->temperature <= 0.0)
    alquimia_error("BatchChemDriver: flow->temperature must be positive.");
  if (input->chemistry_engine == NULL)
    alquimia_error("BatchChemDriver: chemistry->engine not specified.");
  if (input->chemistry_input_file == NULL)
    alquimia_error("BatchChemDriver: chemistry->input_file not specified.");
  if (input->cond_name == NULL)
    alquimia_error("BatchChemDriver: chemistry->initial_condition not specified.");

  // Default description -> input file name.
  if (input->description == NULL)
    input->description = AlquimiaStringDup(input_file);

  // Default output.
  if (input->output_type == NULL)
    input->output_type = AlquimiaStringDup("gnuplot");
  char default_output_file[FILENAME_MAX];
  if (input->output_file == NULL)
  {
    // Find the last '.' in the input filename.
    int dot = strlen(input_file)-1;
    while ((dot > 0) && (input_file[dot] != '.')) --dot;
    char suffix[16];

    // Determine the suffix from the output type.
    if (AlquimiaCaseInsensitiveStringCompare(input->output_type, "gnuplot"))
      sprintf(suffix, ".gnuplot");
    else if (AlquimiaCaseInsensitiveStringCompare(input->output_type, "python"))
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
    if (AlquimiaCaseInsensitiveStringCompare(input->output_type, "python"))
    {
      int len = strlen(default_output_file);
      for (int i = 0; i < len; ++i)
      {
        if (default_output_file[i] == '-')
          default_output_file[i] = '_';
      }
    }

    input->output_file = AlquimiaStringDup(default_output_file);
  }

  return input;
}

void BatchChemDriverInput_Free(BatchChemDriverInput* input)
{
  free(input->description);
  for (int i = 0; i < input->num_isotherm_species; ++i)
    free(input->isotherm_species[i]);
  for (int i = 0; i < input->num_ion_exchange_sites; ++i)
    free(input->ion_exchange_sites[i]);
  for (int i = 0; i < input->num_surface_sites; ++i)
    free(input->surface_sites[i]);
  if (input->cond_name != NULL)
    free(input->cond_name);
  if (input->chemistry_engine != NULL)
    free(input->chemistry_engine);
  if (input->chemistry_input_file != NULL)
    free(input->chemistry_input_file);
  if (input->output_file != NULL)
    free(input->output_file);
  if (input->output_type != NULL)
    free(input->output_type);
  free(input);
}

struct BatchChemDriver 
{
  // Simulation parameters.
  char* description;
  double t_min, t_max, dt;
  int max_steps;
  bool verbose;

  // State variables.
  double water_density, porosity, temperature, aqueous_pressure;

  // Material properties.
  double volume, saturation;

  // Current sim state.
  double time;
  int step;

  // Chemistry data.
  AlquimiaProperties chem_properties;
  AlquimiaState chem_state;
  AlquimiaAuxiliaryData chem_aux_data;
  AlquimiaAuxiliaryOutputData chem_aux_output;

  // Chemistry engine -- one of each of these per thread in general.
  AlquimiaInterface chem;
  void* chem_engine;
  AlquimiaEngineStatus chem_status;

  // Chemistry metadata.
  AlquimiaSizes chem_sizes;
  AlquimiaProblemMetaData chem_metadata;

  // Initial condition.
  AlquimiaGeochemicalCondition chem_cond;

  // Stored history of data.
  int history_size, history_cap;
  double* history_times;
  AlquimiaState* chem_state_history;
  AlquimiaAuxiliaryOutputData* chem_aux_output_history;
};

BatchChemDriver* BatchChemDriver_New(BatchChemDriverInput* input)
{
  BatchChemDriver* driver = malloc(sizeof(BatchChemDriver));

  // Get basic simulation parameters.
  driver->description = AlquimiaStringDup(input->description);
  driver->t_min = input->t_min; 
  driver->t_max = input->t_max;
  driver->max_steps = input->max_steps;
  driver->dt = input->dt;
  driver->verbose = input->verbose;

  // Get chemistry state.
  driver->water_density = input->water_density;
  driver->temperature = input->temperature;
  driver->porosity = input->porosity;
  driver->aqueous_pressure = input->aqueous_pressure;

  // Get material properties.
  driver->volume = input->volume;
  driver->saturation = input->saturation;

  // Simulation state.
  driver->time = driver->t_min;
  driver->step = 0;

  // History.
  driver->history_size = 0;
  driver->history_cap = 32;
  driver->history_times = malloc(sizeof(double) * driver->history_cap);
  driver->chem_state_history = malloc(sizeof(AlquimiaState) * driver->history_cap);
  driver->chem_aux_output_history = malloc(sizeof(AlquimiaAuxiliaryOutputData) * driver->history_cap);

  // Set up the chemistry engine.
  AllocateAlquimiaEngineStatus(&driver->chem_status);
  CreateAlquimiaInterface(input->chemistry_engine, &driver->chem, &driver->chem_status);
  if (driver->chem_status.error != 0) 
  {
    alquimia_error("BatchChemDriver_New: %s", driver->chem_status.message);
    return NULL;
  }

  // Set up the engine and get storage requirements.
  AlquimiaEngineFunctionality chem_engine_functionality;
  driver->chem.Setup(input->chemistry_input_file,
                     input->hands_off,
                     &driver->chem_engine,
                     &driver->chem_sizes,
                     &chem_engine_functionality,
                     &driver->chem_status);
  if (driver->chem_status.error != 0) 
  {
    alquimia_error("BatchChemDriver_New: %s", driver->chem_status.message);
    return NULL;
  }

  // Allocate memory for the chemistry data.
  AllocateAlquimiaProblemMetaData(&driver->chem_sizes, &driver->chem_metadata);
  AllocateAlquimiaState(&driver->chem_sizes, &driver->chem_state);
  AllocateAlquimiaProperties(&driver->chem_sizes, &driver->chem_properties);
  AllocateAlquimiaAuxiliaryData(&driver->chem_sizes, &driver->chem_aux_data);
  AllocateAlquimiaAuxiliaryOutputData(&driver->chem_sizes, &driver->chem_aux_output);

  // Metadata.
  driver->chem.GetProblemMetaData(&driver->chem_engine, 
                                  &driver->chem_metadata, 
                                  &driver->chem_status);
  if (driver->chem_status.error != 0)
  {
    alquimia_error("BatchChemDriver_New: %s", driver->chem_status.message);
    return NULL;
  }

  // Initial condition.
  AllocateAlquimiaGeochemicalCondition(strlen(input->cond_name), 0, 0, &driver->chem_cond);
  strcpy(driver->chem_cond.name, input->cond_name);

  // Copy the chemistry state information in.
  // NOTE: For now, we only allow one of each of these reactions.
  if (input->num_ion_exchange_sites > 0) 
  {
    for (int i = 0; i < driver->chem_state.cation_exchange_capacity.size; ++i)
    {
      int j = IonExchangeIndex(input, driver->chem_metadata.ion_exchange_names.data[i]);
      driver->chem_state.cation_exchange_capacity.data[i] = input->cation_exchange_capacity[j];
    }
    for (int i = 0; i < driver->chem_state.surface_site_density.size; ++i)
    {
      int j = IonExchangeIndex(input, driver->chem_metadata.surface_site_names.data[i]);
      driver->chem_state.surface_site_density.data[i] = input->surface_site_density[j];
    }
  }
  else
  {
    for (int i = 0; i < driver->chem_state.cation_exchange_capacity.size; ++i)
      driver->chem_state.cation_exchange_capacity.data[i] = 0.0;
    for (int i = 0; i < driver->chem_state.surface_site_density.size; ++i)
      driver->chem_state.surface_site_density.data[i] = 0.0;
  }

  if (input->num_isotherm_species > 0) 
  {
    for (int i = 0; i < driver->chem_properties.isotherm_kd.size; ++i)
    {
      int j = IonExchangeIndex(input, driver->chem_metadata.isotherm_species_names.data[i]);
      driver->chem_properties.isotherm_kd.data[i] = input->isotherm_kd[j];
    }
    for (int i = 0; i < driver->chem_properties.langmuir_b.size; ++i)
    {
      int j = IonExchangeIndex(input, driver->chem_metadata.isotherm_species_names.data[i]);
      driver->chem_properties.langmuir_b.data[i] = input->langmuir_b[j];
    }
    for (int i = 0; i < driver->chem_properties.freundlich_n.size; ++i)
    {
      int j = IonExchangeIndex(input, driver->chem_metadata.isotherm_species_names.data[i]);
      driver->chem_properties.freundlich_n.data[i] = input->freundlich_n[j];
    }
  }

  return driver;
}

void BatchChemDriver_Free(BatchChemDriver* driver)
{
  free(driver->description);

  // Destroy initial condition.
  FreeAlquimiaGeochemicalCondition(&driver->chem_cond);

  // Destroy chemistry data.
  FreeAlquimiaState(&driver->chem_state);
  FreeAlquimiaProperties(&driver->chem_properties);
  FreeAlquimiaAuxiliaryData(&driver->chem_aux_data);
  FreeAlquimiaAuxiliaryOutputData(&driver->chem_aux_output);
  FreeAlquimiaProblemMetaData(&driver->chem_metadata);

  // Destroy history data.
  for (int i = 0; i < driver->history_size; ++i)
  {
    FreeAlquimiaState(&driver->chem_state_history[i]);
    FreeAlquimiaAuxiliaryOutputData(&driver->chem_aux_output_history[i]);
  }
  free(driver->history_times);
  free(driver->chem_state_history);
  free(driver->chem_aux_output_history);

  // Destroy chemistry engine.
  driver->chem.Shutdown(&driver->chem_engine, 
                        &driver->chem_status);
  FreeAlquimiaEngineStatus(&driver->chem_status);

  free(driver);
}

static void BatchChemDriver_RecordHistory(BatchChemDriver* driver)
{
  if (driver->history_size == driver->history_cap)
  {
    driver->history_cap *= 2;
    driver->history_times = realloc(driver->history_times, sizeof(double) * driver->history_cap);
    driver->chem_state_history = realloc(driver->chem_state_history, 
                                         sizeof(AlquimiaState) * driver->history_cap);
    driver->chem_aux_output_history = realloc(driver->chem_aux_output_history, 
                                              sizeof(AlquimiaAuxiliaryOutputData) * driver->history_cap);
  }
  driver->history_times[driver->history_size] = driver->time;
  AllocateAlquimiaState(&driver->chem_sizes, &driver->chem_state_history[driver->history_size]);
  CopyAlquimiaState(&driver->chem_state, &driver->chem_state_history[driver->history_size]);
  AllocateAlquimiaAuxiliaryOutputData(&driver->chem_sizes, &driver->chem_aux_output_history[driver->history_size]);
  CopyAlquimiaAuxiliaryOutputData(&driver->chem_aux_output, &driver->chem_aux_output_history[driver->history_size]);
  ++driver->history_size;
}

static int BatchChemDriver_Initialize(BatchChemDriver* driver)
{
  if (driver->verbose)
    printf("BatchChemDriver: initializing at time %g...\n", driver->t_min);

  // Initialize the batch.

  // Set the material properties.
  driver->chem_properties.volume = driver->volume;
  driver->chem_properties.saturation = driver->saturation;

  // Set the thermodynamic state.
  driver->chem_state.water_density = driver->water_density;
  driver->chem_state.temperature = driver->temperature;
  driver->chem_state.porosity = driver->porosity;
  driver->chem_state.aqueous_pressure = driver->aqueous_pressure;

  // Invoke the chemical initial condition.
  driver->chem.ProcessCondition(&driver->chem_engine,
                                &driver->chem_cond, 
                                &driver->chem_properties,
                                &driver->chem_state,
                                &driver->chem_aux_data,
                                &driver->chem_status);
  if (driver->chem_status.error != 0)
  {
    alquimia_error("BatchChemDriver: initialization error: %s\n", driver->chem_status.message);
  }

  // Overwrite the stuff that might have changed.
  driver->chem_state.porosity = driver->porosity;

  // Append our state and auxiliary output to our history.
  BatchChemDriver_RecordHistory(driver);

  if (driver->verbose)
    printf("BatchChemDriver: Finished initializing.\n");

  return driver->chem_status.error;
}

int BatchChemDriver_Run(BatchChemDriver* driver)
{
  if (driver->verbose)
    printf("BatchChemDriver: running %s\n", driver->description);

  double dt = driver->dt;

  // Initialize the chemistry state and set up the solution vector.
  int status = BatchChemDriver_Initialize(driver);
  if (status != 0)
    return status;

  driver->time = driver->t_min;
  driver->step = 0;
  while ((driver->time < driver->t_max) && (driver->step < driver->max_steps))
  {
    if (driver->verbose)
      printf("BatchChemDriver: step %d (t = %g, dt = %g)\n", driver->step, driver->time, dt);

    // Do the chemistry step.
    driver->chem.ReactionStepOperatorSplit(&driver->chem_engine,
                                           dt, &driver->chem_properties,
                                           &driver->chem_state,
                                           &driver->chem_aux_data,
                                           &driver->chem_status);
    if (driver->chem_status.error != 0)
    {
      status = driver->chem_status.error;
      printf("BatchChemDriver: reaction failed: %s\n", driver->chem_status.message);
      break;
    }

    // Fetch auxiliary output.
    driver->chem.GetAuxiliaryOutput(&driver->chem_engine, 
                                    &driver->chem_properties,
                                    &driver->chem_state,
                                    &driver->chem_aux_data,
                                    &driver->chem_aux_output,
                                    &driver->chem_status);
    if (driver->chem_status.error != 0)
    {
      status = driver->chem_status.error;
      printf("BatchChemDriver: auxiliary output fetch failed: %s\n", driver->chem_status.message);
      break;
    }

    // FIXME: Why is the porosity being set to zero here??
    driver->chem_state.porosity = driver->porosity;

    driver->time += dt;
    driver->step += 1;

    // Append our state and auxiliary output to our history.
    BatchChemDriver_RecordHistory(driver);
  }

  return status;
}

void BatchChemDriver_GetSoluteAndAuxData(BatchChemDriver* driver,
                                         double* time,
                                         AlquimiaVectorString* var_names,
                                         AlquimiaVectorDouble* var_data)
{
  // Destroy the contents of the vectors we're given.
  if (var_names->size > 0)
    FreeAlquimiaVectorString(var_names);
  if (var_data->size > 0)
    FreeAlquimiaVectorDouble(var_data);

  // Construct a list of all variables, which are those in the state and 
  // the auxiliary output data, and fill their data.
  int num_times = driver->history_size;
  int num_primary = driver->chem_sizes.num_primary;
  int num_sorbed = driver->chem_sizes.num_sorbed;
  int num_minerals = driver->chem_sizes.num_minerals;
  int num_surface_sites = driver->chem_sizes.num_surface_sites;
  int num_ion_exchange_sites = driver->chem_sizes.num_ion_exchange_sites;
  int num_aqueous_complexes = driver->chem_sizes.num_aqueous_complexes;
  int num_aqueous_kinetics = driver->chem_sizes.num_aqueous_kinetics;
  int num_vars = 1 +                        // time
                 num_primary +              // total mobile
                 num_sorbed +               // total immobile
                 2 * num_minerals +         // mineral volume fractions, specific surface area
                 num_surface_sites +        // surface site density
                 num_ion_exchange_sites +   // cation exchange capacity
                 1 +                        // pH
                 num_aqueous_kinetics +     // aqueous kinetic rate
                 2 * num_minerals +         // mineral saturation index, reaction rate
                 2 * num_primary +          // primary free ion concentration, activity coeff
                 2 * num_aqueous_complexes; // secondary free ion concentration, activity coeff
  int counter = 0;
  AllocateAlquimiaVectorString(num_vars, var_names);
  AllocateAlquimiaVectorDouble(num_vars * num_times, var_data);
  {
    var_names->data[counter] = AlquimiaStringDup("time");
    for (int j = 0; j < driver->history_size; ++j)
      var_data->data[num_vars*j + counter] = driver->history_times[j];
    ++counter;
  }
  for (int i = 0; i < num_primary; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "total_mobile[%s]", driver->chem_metadata.primary_names.data[i]);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < driver->history_size; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_state_history[j].total_mobile.data[i];
  }
  for (int i = 0; i < num_sorbed; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "total_immobile[%d]", i);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < driver->history_size; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_state_history[j].total_immobile.data[i];
  }
  for (int i = 0; i < num_minerals; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "mineral_volume_fractions[%s]", driver->chem_metadata.mineral_names.data[i]);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < driver->history_size; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_state_history[j].mineral_volume_fraction.data[i];
  }
  for (int i = 0; i < num_minerals; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "mineral_specific_surface_area[%s]", driver->chem_metadata.mineral_names.data[i]);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < driver->history_size; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_state_history[j].mineral_specific_surface_area.data[i];
  }
  for (int i = 0; i < num_surface_sites; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "surface_site_density[%s]", driver->chem_metadata.surface_site_names.data[i]);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < driver->history_size; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_state_history[j].surface_site_density.data[i];
  }
  for (int i = 0; i < num_ion_exchange_sites; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "cation_exchange_capacity[%s]", driver->chem_metadata.ion_exchange_names.data[i]);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < driver->history_size; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_state_history[j].cation_exchange_capacity.data[i];
  }
  {
    var_names->data[counter] = AlquimiaStringDup("pH");
    for (int j = 0; j < driver->history_size; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_aux_output.pH;
    ++counter;
  }
  for (int i = 0; i < num_aqueous_kinetics; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "aqueous_kinetic_rate[%s]", driver->chem_metadata.aqueous_kinetic_names.data[i]);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < driver->history_size; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_aux_output_history[j].aqueous_kinetic_rate.data[i];
  }
  for (int i = 0; i < num_minerals; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "mineral_saturation_index[%s]", driver->chem_metadata.mineral_names.data[i]);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < driver->history_size; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_aux_output_history[j].mineral_saturation_index.data[i];
  }
  for (int i = 0; i < num_minerals; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "mineral_reaction_rate[%s]", driver->chem_metadata.mineral_names.data[i]);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < driver->history_size; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_aux_output_history[j].mineral_reaction_rate.data[i];
  }
  for (int i = 0; i < num_primary; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "primary_free_ion_concentration[%s]", driver->chem_metadata.primary_names.data[i]);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < driver->history_size; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_aux_output_history[j].primary_free_ion_concentration.data[i];
  }
  for (int i = 0; i < num_primary; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "primary_activity_coeff[%s]", driver->chem_metadata.primary_names.data[i]);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < driver->history_size; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_aux_output_history[j].primary_activity_coeff.data[i];
  }
  for (int i = 0; i < num_aqueous_complexes; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "secondary_free_ion_concentration[%d]", i);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < driver->history_size; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_aux_output_history[j].secondary_free_ion_concentration.data[i];
  }
  for (int i = 0; i < num_aqueous_complexes; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "secondary_activity_coeff[%d]", i);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < driver->history_size; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_aux_output_history[j].secondary_activity_coeff.data[i];
  }
}

