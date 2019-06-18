/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

//
// Alquimia Copyright (c) 2013-2016, The Regents of the University of California, 
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
#include "TransportDriver.h"

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

// Input parsing stuff. See https://github.com/benhoyt/inih for details.
static int ParseInput(void* user,
                      const char* section,
                      const char* name,
                      const char* value)
{
  TransportDriverInput* input = user;
#define MATCH(s, n) (strcmp(section, s) == 0) && (strcmp(name, n) == 0)

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
  else if (MATCH("simulation", "cfl_factor"))
    input->cfl_factor = atof(value);

  // Domain section.
  else if (MATCH("domain", "x_min"))
    input->x_min = atof(value);
  else if (MATCH("domain", "x_max"))
    input->x_max = atof(value);
  else if (MATCH("domain", "num_cells"))
    input->num_cells = atoi(value);

  // Material section.
  else if (MATCH("material", "porosity"))
    input->porosity = atof(value);
  else if (MATCH("material", "saturation"))
    input->saturation = atof(value);

  // Flow section.
  else if (MATCH("flow", "temperature"))
    input->temperature = atof(value);
  else if (MATCH("flow", "velocity"))
    input->velocity = atof(value);

  // Chemistry section.
  else if (MATCH("chemistry", "engine"))
    input->chemistry_engine = AlquimiaStringDup(value);
  else if (MATCH("chemistry", "input_file"))
    input->chemistry_input_file = AlquimiaStringDup(value);
  else if (MATCH("chemistry", "initial_condition"))
    input->ic_name = AlquimiaStringDup(value);
  else if (MATCH("chemistry", "left_boundary_condition"))
    input->left_bc_name = AlquimiaStringDup(value);
  else if (MATCH("chemistry", "right_boundary_condition"))
    input->right_bc_name = AlquimiaStringDup(value);
  else if (MATCH("chemistry", "cation_exchange_capacity"))
    input->cation_exchange_capacity = atof(value);
  else if (MATCH("chemistry", "surface_site_density"))
    input->surface_site_density = atof(value);

  // Output section.
  else if (MATCH("output", "verbose"))
    input->verbose = ValueIsTrue(value);
  else if (MATCH("output", "type"))
    input->output_type = AlquimiaStringDup(value);
  else if (MATCH("output", "filename"))
    input->output_file = AlquimiaStringDup(value);

  return 1;
}

TransportDriverInput* TransportDriverInput_New(const char* input_file)
{
  TransportDriverInput* input = malloc(sizeof(TransportDriverInput));
  memset(input, 0, sizeof(TransportDriverInput));

  // Make sure we have some meaningful defaults.
  input->hands_off = true; // Hands-off by default.
  input->t_min = 0.0;
  input->max_steps = INT_MAX;
  input->dt = FLT_MAX;
  input->cfl_factor = 1.0;
  input->porosity = 1.0;
  input->saturation = 1.0;
  input->temperature = 25.0;
  input->velocity = 0.0;
  input->left_bc_name = NULL;
  input->right_bc_name = NULL;
  input->cation_exchange_capacity = 0.0;
  input->surface_site_density = 0.0;
  input->output_type = NULL;
  input->output_file = NULL;

  // Fill in fields by parsing the input file.
  int error = ini_parse(input_file, ParseInput, input);
  if (error != 0)
  {
    TransportDriverInput_Free(input);
    alquimia_error("TransportDriver: Error parsing input: %s", input_file);
  }

  // Verify that our required fields are filled properly.
  if (!input->hands_off)
    alquimia_error("TransportDriver: simulation->hands_off must be set to true at the moment.");
  if (input->t_max <= input->t_min)
    alquimia_error("TransportDriver: simulation->t_max must be greater than simulation->t_min.");
  if (input->max_steps < 0)
    alquimia_error("TransportDriver: simulation->max_steps must be non-negative.");
  if ((input->cfl_factor <= 0.0) || (input->cfl_factor > 1.0))
    alquimia_error("TransportDriver: simulation->cfl_factor must be within (0, 1].");
  if (input->x_max <= input->x_min)
    alquimia_error("TransportDriver: domain->x_max must be greater than domain->x_min.");
  if (input->num_cells <= 0)
    alquimia_error("TransportDriver: domain->num_cells must be positive.");
  if ((input->porosity <= 0.0) || (input->porosity > 1.0))
    alquimia_error("TransportDriver: material->porosity must be within (0, 1].");
  if ((input->saturation <= 0.0) || (input->saturation > 1.0))
    alquimia_error("TransportDriver: material->saturation must be within (0, 1].");
  if (input->temperature <= 0.0)
    alquimia_error("TransportDriver: flow->temperature must be positive.");
  if ((input->left_bc_name == NULL) && (input->right_bc_name == NULL) && (input->velocity != 0.0))
    alquimia_error("TransportDriver: When velocity != 0, left or right boundary condition must be given.");
  if (input->chemistry_engine == NULL)
    alquimia_error("TransportDriver: chemistry->engine not specified.");
  if (input->chemistry_input_file == NULL)
    alquimia_error("TransportDriver: chemistry->input_file not specified.");

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
    suffix[0] = '\0';

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

void TransportDriverInput_Free(TransportDriverInput* input)
{
  free(input->description);
  if (input->ic_name != NULL)
    free(input->ic_name);
  if (input->left_bc_name != NULL)
    free(input->left_bc_name);
  if (input->right_bc_name != NULL)
    free(input->right_bc_name);
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

struct TransportDriver 
{
  // Simulation parameters.
  char* description;
  TransportCoupling coupling;
  double t_min, t_max, dt, cfl;
  int max_steps;
  bool verbose;

  // 1D grid information.
  int num_cells;
  double x_min, x_max;

  // Flow velocity, temperature.
  double vx, temperature;

  // Material properties.
  double porosity, saturation;

  // Current sim state.
  double time;
  int step;

  // Per-cell chemistry data.
  AlquimiaProperties* chem_properties;
  AlquimiaState* chem_state;
  AlquimiaAuxiliaryData* chem_aux_data;
  AlquimiaAuxiliaryOutputData* chem_aux_output;

  // Chemistry engine -- one of each of these per thread in general.
  AlquimiaInterface chem;
  void* chem_engine;
  AlquimiaEngineStatus chem_status;

  // Chemistry metadata.
  AlquimiaSizes chem_sizes;
  AlquimiaProblemMetaData chem_metadata;

  // Initial and boundary conditions.
  AlquimiaGeochemicalCondition chem_ic;
  AlquimiaGeochemicalCondition chem_left_bc, chem_right_bc;
  AlquimiaState chem_left_state, chem_right_state;
  AlquimiaAuxiliaryData chem_left_aux_data, chem_right_aux_data;

  // Bookkeeping.
  AlquimiaState advected_chem_state;
  AlquimiaAuxiliaryData advected_chem_aux_data;
  double* advective_fluxes;
};

TransportDriver* TransportDriver_New(TransportDriverInput* input)
{
  TransportDriver* driver = malloc(sizeof(TransportDriver));

  // Get basic simulation parameters.
  driver->description = AlquimiaStringDup(input->description);
  driver->coupling = input->coupling;
  driver->t_min = input->t_min; 
  driver->t_max = input->t_max;
  driver->max_steps = input->max_steps;
  driver->dt = input->dt;
  driver->cfl = input->cfl_factor;
  driver->verbose = input->verbose;

  // Get grid information.
  driver->x_min = input->x_min; 
  driver->x_max = input->x_max;
  driver->num_cells = input->num_cells;

  // Get material properties.
  driver->porosity = input->porosity;
  driver->saturation = input->saturation;

  // Get flow variables.
  driver->vx = input->velocity;
  driver->temperature = input->temperature;

  // Simulation state.
  driver->time = driver->t_min;
  driver->step = 0;

  // Set up the chemistry engine.
  AllocateAlquimiaEngineStatus(&driver->chem_status);
  CreateAlquimiaInterface(input->chemistry_engine, &driver->chem, &driver->chem_status);
  if (driver->chem_status.error != 0) 
  {
    alquimia_error("TransportDriver_New: %s", driver->chem_status.message);
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
    alquimia_error("TransportDriver_New: %s", driver->chem_status.message);
    return NULL;
  }

  // If you want multiple copies of the chemistry engine with
  // OpenMP, verify: chem_data.functionality.thread_safe == true,
  // then create the appropriate number of chem status and data
  // objects.

  // Allocate memory for the chemistry data.
  AllocateAlquimiaProblemMetaData(&driver->chem_sizes, &driver->chem_metadata);
  driver->chem_properties = malloc(sizeof(AlquimiaProperties) * driver->num_cells);
  driver->chem_state = malloc(sizeof(AlquimiaState) * driver->num_cells);
  driver->chem_aux_data = malloc(sizeof(AlquimiaAuxiliaryData) * driver->num_cells);
  driver->chem_aux_output = malloc(sizeof(AlquimiaAuxiliaryOutputData) * driver->num_cells);
  for (int i = 0; i < driver->num_cells; ++i)
  {
    AllocateAlquimiaState(&driver->chem_sizes, &driver->chem_state[i]);
    AllocateAlquimiaProperties(&driver->chem_sizes, &driver->chem_properties[i]);
    AllocateAlquimiaAuxiliaryData(&driver->chem_sizes, &driver->chem_aux_data[i]);
    AllocateAlquimiaAuxiliaryOutputData(&driver->chem_sizes, &driver->chem_aux_output[i]);
  }

  // Metadata.
  driver->chem.GetProblemMetaData(&driver->chem_engine, 
                                  &driver->chem_metadata, 
                                  &driver->chem_status);
  if (driver->chem_status.error != 0)
  {
    alquimia_error("TransportDriver_New: %s", driver->chem_status.message);
    return NULL;
  }

  // Initial condition.
  AllocateAlquimiaGeochemicalCondition(strlen(input->ic_name), 0, 0, &driver->chem_ic);
  strcpy(driver->chem_ic.name, input->ic_name);

  // Boundary conditions.
  if (input->left_bc_name != NULL)
  {
    AllocateAlquimiaGeochemicalCondition(strlen(input->left_bc_name), 0, 0, &driver->chem_left_bc);
    strcpy(driver->chem_left_bc.name, input->left_bc_name);
  }
  else
    driver->chem_left_bc.name = NULL;
  AllocateAlquimiaState(&driver->chem_sizes, &driver->chem_left_state);
  AllocateAlquimiaAuxiliaryData(&driver->chem_sizes, &driver->chem_left_aux_data);
  if (input->right_bc_name != NULL)
  {
    AllocateAlquimiaGeochemicalCondition(strlen(input->right_bc_name), 0, 0, &driver->chem_right_bc);
    strcpy(driver->chem_right_bc.name, input->right_bc_name);
  }
  else
    driver->chem_right_bc.name = NULL;
  AllocateAlquimiaState(&driver->chem_sizes, &driver->chem_right_state);
  AllocateAlquimiaAuxiliaryData(&driver->chem_sizes, &driver->chem_right_aux_data);

  // Copy the miscellaneous chemistry state information in.
  // NOTE: For now, we only allow one of each of these reactions.
  for (int i = 0; i < driver->num_cells; ++i)
  {
    for (int j = 0; j < driver->chem_state[i].cation_exchange_capacity.size; ++j)
      driver->chem_state[i].cation_exchange_capacity.data[j] = input->cation_exchange_capacity;
    for (int j = 0; j < driver->chem_state[i].surface_site_density.size; ++j)
      driver->chem_state[i].surface_site_density.data[j] = input->surface_site_density;
  }
  for (int j = 0; j < driver->chem_left_state.cation_exchange_capacity.size; ++j)
    driver->chem_left_state.cation_exchange_capacity.data[j] = input->cation_exchange_capacity;
  for (int j = 0; j < driver->chem_left_state.surface_site_density.size; ++j)
    driver->chem_left_state.surface_site_density.data[j] = input->surface_site_density;
  for (int j = 0; j < driver->chem_right_state.cation_exchange_capacity.size; ++j)
    driver->chem_right_state.cation_exchange_capacity.data[j] = input->cation_exchange_capacity;
  for (int j = 0; j < driver->chem_right_state.surface_site_density.size; ++j)
    driver->chem_right_state.surface_site_density.data[j] = input->surface_site_density;

  // Bookkeeping.
  AllocateAlquimiaState(&driver->chem_sizes, &driver->advected_chem_state);
  AllocateAlquimiaAuxiliaryData(&driver->chem_sizes, &driver->advected_chem_aux_data);
  driver->advective_fluxes = malloc(sizeof(double) * driver->chem_sizes.num_primary * (driver->num_cells + 1));

  return driver;
}

void TransportDriver_Free(TransportDriver* driver)
{
  free(driver->description);

  // Destroy advection bookkeeping.
  free(driver->advective_fluxes);
  FreeAlquimiaState(&driver->advected_chem_state);
  FreeAlquimiaAuxiliaryData(&driver->advected_chem_aux_data);

  // Destroy boundary and initial conditions.
  if (driver->chem_left_bc.name != NULL)
    FreeAlquimiaGeochemicalCondition(&driver->chem_left_bc);
  FreeAlquimiaState(&driver->chem_left_state);
  FreeAlquimiaAuxiliaryData(&driver->chem_left_aux_data);
  if (driver->chem_right_bc.name != NULL)
    FreeAlquimiaGeochemicalCondition(&driver->chem_right_bc);
  FreeAlquimiaState(&driver->chem_right_state);
  FreeAlquimiaAuxiliaryData(&driver->chem_right_aux_data);
  FreeAlquimiaGeochemicalCondition(&driver->chem_ic);

  // Destroy chemistry data.
  for (int i = 0; i < driver->num_cells; ++i)
  {
    FreeAlquimiaState(&driver->chem_state[i]);
    FreeAlquimiaProperties(&driver->chem_properties[i]);
    FreeAlquimiaAuxiliaryData(&driver->chem_aux_data[i]);
    FreeAlquimiaAuxiliaryOutputData(&driver->chem_aux_output[i]);
  }
  free(driver->chem_state);
  free(driver->chem_properties);
  free(driver->chem_aux_data);
  free(driver->chem_aux_output);
  FreeAlquimiaProblemMetaData(&driver->chem_metadata);

  // Destroy chemistry engine.
  driver->chem.Shutdown(&driver->chem_engine, 
                        &driver->chem_status);
  FreeAlquimiaEngineStatus(&driver->chem_status);
  free(driver);
}

static int TransportDriver_Initialize(TransportDriver* driver)
{
  if (driver->verbose)
    printf("TransportDriver: initializing at time %g...\n", driver->t_min);

  static const double water_density = 999.9720;    // density of water in kg/m**3
  static const double aqueous_pressure = 201325.0; // pressure in Pa.

  // Determine the volume of a cell from the domain.
  double volume = (driver->x_max - driver->x_min) / driver->num_cells;

  // Initialize each cell.
  for (int i = 0; i < driver->num_cells; ++i)
  {
    if (i == 0)
    {
      // Set the material properties.
      driver->chem_properties[i].volume = volume;
      driver->chem_properties[i].saturation = driver->saturation;

      // Set the thermodynamic state.
      driver->chem_state[i].water_density = water_density;
      driver->chem_state[i].temperature = driver->temperature;
      driver->chem_state[i].porosity = driver->porosity;
      driver->chem_state[i].aqueous_pressure = aqueous_pressure;

      // Invoke the chemical initial condition.
      driver->chem.ProcessCondition(&driver->chem_engine,
                                    &driver->chem_ic, 
                                    &driver->chem_properties[i],
                                    &driver->chem_state[i],
                                    &driver->chem_aux_data[i],
                                    &driver->chem_status);
      if (driver->chem_status.error != 0)
      {
        printf("TransportDriver: initialization error: %s\n", 
               driver->chem_status.message);
        break;
      }

      // Overwrite the stuff that might have changed.
      driver->chem_state[i].porosity = driver->porosity;
    }
    else
    {
      // Just copy the contents of cell 0 over.
      CopyAlquimiaProperties(&driver->chem_properties[0], &driver->chem_properties[i]);
      CopyAlquimiaState(&driver->chem_state[0], &driver->chem_state[i]);
      CopyAlquimiaAuxiliaryData(&driver->chem_aux_data[0], &driver->chem_aux_data[i]);
    }
  }

  // Initialize the left and right boundary states.
  int num_primary = driver->chem_sizes.num_primary;

  // Left.
  driver->chem_left_state.water_density = water_density;
  driver->chem_left_state.temperature = driver->temperature;
  driver->chem_left_state.porosity = driver->porosity;
  driver->chem_left_state.aqueous_pressure = aqueous_pressure;
  if (driver->chem_left_bc.name != NULL)
  {
    driver->chem.ProcessCondition(&driver->chem_engine,
                                  &driver->chem_left_bc, 
                                  &driver->chem_properties[0],
                                  &driver->chem_left_state,
                                  &driver->chem_left_aux_data,
                                  &driver->chem_status);
    if (driver->chem_status.error != 0)
    {
      printf("TransportDriver: boundary condition error at leftmost interface: %s\n",
             driver->chem_status.message);
      return driver->chem_status.error;
    }

    // Overwrite things that might have changed.
    // FIXME: Why??
    driver->chem_left_state.porosity = driver->porosity;
  }
  else
  {
    for (int c = 0; c < num_primary; ++c)
      driver->chem_left_state.total_mobile.data[c] = driver->chem_state[0].total_mobile.data[c];
  }

  // Right.
  driver->chem_right_state.water_density = water_density;
  driver->chem_right_state.temperature = driver->temperature;
  driver->chem_right_state.porosity = driver->porosity;
  driver->chem_right_state.aqueous_pressure = aqueous_pressure;
  if (driver->chem_right_bc.name != NULL)
  {
    driver->chem.ProcessCondition(&driver->chem_engine,
                                  &driver->chem_right_bc, 
                                  &driver->chem_properties[driver->num_cells-1],
                                  &driver->chem_right_state,
                                  &driver->chem_right_aux_data,
                                  &driver->chem_status);
    if (driver->chem_status.error != 0)
    {
      printf("TransportDriver: boundary condition error at rightmost interface: %s\n",
             driver->chem_status.message);
      return driver->chem_status.error;
    }

    // Overwrite things that might have changed.
    // FIXME: Why??
    driver->chem_right_state.porosity = driver->porosity;
  }
  else
  {
    for (int c = 0; c < num_primary; ++c)
      driver->chem_right_state.total_mobile.data[c] = driver->chem_state[driver->num_cells-1].total_mobile.data[c];
  }

  if (driver->verbose)
    printf("TransportDriver: Finished initializing.\n");

  return driver->chem_status.error;
}

// Here's a finite volume method for advecting concentrations using the 
// 1D flow velocity ux in the given cell.
static int ComputeAdvectiveFluxes(TransportDriver* driver, 
                                  double dx,
                                  double* advective_fluxes)
{
  int num_primary = driver->chem_sizes.num_primary;

  // Estimate the concentrations on each interface, beginning with those 
  // at the boundaries. Note that properties are taken from the leftmost
  // and rightmost cells. Then compute fluxes.

  // Left boundary.
  double phi_l = driver->chem_left_state.porosity;
  for (int c = 0; c < num_primary; ++c)
    advective_fluxes[c] = phi_l * driver->vx * driver->chem_left_state.total_mobile.data[c];

  // Right boundary.
  double phi_r = driver->chem_right_state.porosity;
  for (int c = 0; c < num_primary; ++c)
    advective_fluxes[num_primary*driver->num_cells+c] = phi_r * driver->vx * driver->chem_right_state.total_mobile.data[c];

  // Interior interfaces.
  for (int i = 0; i < driver->num_cells-1; ++i)
  {
    for (int c = 0; c < num_primary; ++c)
    {
      // Figure out the upstream and downstream concentrations.
      double up_conc, down_conc;
      if (driver->vx >= 0.0)
      {
        up_conc = driver->chem_state[i].total_mobile.data[c];
        down_conc = driver->chem_state[i+1].total_mobile.data[c];
      }
      else
      {
        up_conc = driver->chem_state[i+1].total_mobile.data[c];
        down_conc = driver->chem_state[i].total_mobile.data[c];
      }

      // We use a simple upwinding method to determine the interface
      // concentration.
      double interface_conc = up_conc;
#if 0
      // We use the method of Gupta et al, 1991 to construct an explicit, 
      // third-order TVD estimate of the concentration at the interface.

      // Construct the limiter.
      double conc = driver->chem_state[i].total_mobile.data[c];
      double r_num = (conc - up_conc) / (2.0 * dx);
      double r_denom = (down_conc - conc) / (2.0 * dx);
      double r = (r_denom != 0.0) ? (r_num / r_denom) : 0.0;
      double beta = Max(0.0, Min(2.0, Min(2.0*r, (2.0 + r) / 3.0)));

      // Compute the interface concentration.
      double interface_conc = conc + 0.5 * beta * (down_conc - conc);
#endif

      // Compute the interface porosity.
      double phi = 0.5 * (driver->chem_state[i].porosity + driver->chem_state[i+1].porosity);

      // Compute the flux.
      advective_fluxes[num_primary * (i+1) + c] = phi * driver->vx * interface_conc;
    }
  }

  return 0;
}

static int Run_OperatorSplit(TransportDriver* driver)
{
  double dx = (driver->x_max - driver->x_min) / driver->num_cells;
  double dt = Min(driver->dt, driver->cfl * dx / driver->vx);
  int num_primary = driver->chem_sizes.num_primary;

  // Initialize the chemistry state in each cell, and set up the solution vector.
  int status = TransportDriver_Initialize(driver);
  if (status != 0)
    return status;

  driver->time = driver->t_min;
  driver->step = 0;
  while ((driver->time < driver->t_max) && (driver->step < driver->max_steps))
  {
    if (driver->verbose)
      printf("TransportDriver: step %d (t = %g, dt = %g)\n", driver->step, driver->time, dt);

    // Compute the advective fluxes.
    if (driver->vx != 0.0)
    {
      status = ComputeAdvectiveFluxes(driver, dx, driver->advective_fluxes);
      if (status != 0) break;
    }

    for (int i = 0; i < driver->num_cells; ++i)
    {
      double phi = driver->chem_state[i].porosity;

      // Advect the state in this cell using a finite volume method with 
      // the advective fluxes we've computed.
      CopyAlquimiaState(&driver->chem_state[i], &driver->advected_chem_state);
      CopyAlquimiaAuxiliaryData(&driver->chem_aux_data[i], &driver->advected_chem_aux_data);
      for (int c = 0; c < num_primary; ++c)
      {
        double F_left = driver->advective_fluxes[num_primary*i+c];
        double F_right = driver->advective_fluxes[num_primary*(i+1)+c];
//printf("F_left = %g, F_right = %g\n", F_left, F_right);
        double phiC = phi * driver->advected_chem_state.total_mobile.data[c] -
                      dt * (F_right - F_left) / dx;
        driver->advected_chem_state.total_mobile.data[c] = phiC / phi;
        driver->advected_chem_state.total_mobile.data[c] = Max(driver->advected_chem_state.total_mobile.data[c], 1e-20); // Floor it at a tiny positive value.
//printf("phi = %g, C = %g, phiC = %g, phi* = %g\n", phi, driver->advected_chem_state.total_mobile.data[c], phiC, driver->advected_chem_state.total_mobile.data[c]);
      }

      // Do the chemistry step, using the advected state as input.
      driver->chem.ReactionStepOperatorSplit(&driver->chem_engine,
                                             dt, &driver->chem_properties[i],
                                             &driver->advected_chem_state,
                                             &driver->advected_chem_aux_data,
                                             &driver->chem_status);
      if (driver->chem_status.error != 0)
      {
        status = driver->chem_status.error;
        printf("TransportDriver: operator-split reaction in cell %d: %s\n", 
               i, driver->chem_status.message);
        break;
      }

      // FIXME: Why is the porosity being set to zero here??
      driver->advected_chem_state.porosity = driver->porosity;

      // Copy the advected/reacted state back into place.
      CopyAlquimiaState(&driver->advected_chem_state, &driver->chem_state[i]);
      CopyAlquimiaAuxiliaryData(&driver->advected_chem_aux_data, &driver->chem_aux_data[i]);

      // Fetch auxiliary output.
      driver->chem.GetAuxiliaryOutput(&driver->chem_engine, 
                                      &driver->chem_properties[i],
                                      &driver->chem_state[i],
                                      &driver->chem_aux_data[i],
                                      &driver->chem_aux_output[i],
                                      &driver->chem_status);
      if (driver->chem_status.error != 0)
      {
        status = driver->chem_status.error;
        printf("TransportDriver: auxiliary output fetch failed: %s\n", driver->chem_status.message);
        break;
      }
    }

    if (status != 0) break;

    driver->time += dt;
    driver->step += 1;
  }

  return status;
}

static int Run_GlobalImplicit(TransportDriver* driver)
{
  alquimia_error("Globally implicit transport is not yet supported.");
  return -1;
}

int TransportDriver_Run(TransportDriver* driver)
{
  if (driver->verbose)
    printf("TransportDriver: running %s\n", driver->description);
  if (driver->coupling == TRANSPORT_OPERATOR_SPLIT)
    return Run_OperatorSplit(driver);
  else
    return Run_GlobalImplicit(driver);
}

void TransportDriver_GetSoluteAndAuxData(TransportDriver* driver,
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
  int num_cells = driver->num_cells;
  int num_primary = driver->chem_sizes.num_primary;
  int num_sorbed = driver->chem_sizes.num_sorbed;
  int num_minerals = driver->chem_sizes.num_minerals;
  int num_surface_sites = driver->chem_sizes.num_surface_sites;
  int num_ion_exchange_sites = driver->chem_sizes.num_ion_exchange_sites;
  int num_aqueous_complexes = driver->chem_sizes.num_aqueous_complexes;
  int num_aqueous_kinetics = driver->chem_sizes.num_aqueous_kinetics;
  int num_vars = 1 +                        // grid cell locations 
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
  AllocateAlquimiaVectorDouble(num_vars * driver->num_cells, var_data);
  {
    var_names->data[counter] = AlquimiaStringDup("x");
    for (int j = 0; j < num_cells; ++j)
      var_data->data[num_vars*j + counter] = driver->x_min + (j+0.5) * (driver->x_max - driver->x_min) / driver->num_cells;
    ++counter;
  }
  for (int i = 0; i < num_primary; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "total_mobile[%s]", driver->chem_metadata.primary_names.data[i]);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < num_cells; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_state[j].total_mobile.data[i];
  }
  for (int i = 0; i < num_sorbed; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "total_immobile[%d]", i);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < num_cells; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_state[j].total_immobile.data[i];
  }
  for (int i = 0; i < num_minerals; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "mineral_volume_fractions[%s]", driver->chem_metadata.mineral_names.data[i]);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < num_cells; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_state[j].mineral_volume_fraction.data[i];
  }
  for (int i = 0; i < num_minerals; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "mineral_specific_surface_area[%s]", driver->chem_metadata.mineral_names.data[i]);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < num_cells; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_state[j].mineral_specific_surface_area.data[i];
  }
  for (int i = 0; i < num_surface_sites; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "surface_site_density[%s]", driver->chem_metadata.surface_site_names.data[i]);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < num_cells; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_state[j].surface_site_density.data[i];
  }
  for (int i = 0; i < num_ion_exchange_sites; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "cation_exchange_capacity[%s]", driver->chem_metadata.ion_exchange_names.data[i]);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < num_cells; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_state[j].cation_exchange_capacity.data[i];
  }
  {
    var_names->data[counter] = AlquimiaStringDup("pH");
    for (int j = 0; j < num_cells; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_aux_output[j].pH;
    ++counter;
  }
  for (int i = 0; i < num_aqueous_kinetics; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "aqueous_kinetic_rate[%s]", driver->chem_metadata.aqueous_kinetic_names.data[i]);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < num_cells; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_aux_output[j].aqueous_kinetic_rate.data[i];
  }
  for (int i = 0; i < num_minerals; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "mineral_saturation_index[%s]", driver->chem_metadata.mineral_names.data[i]);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < num_cells; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_aux_output[j].mineral_saturation_index.data[i];
  }
  for (int i = 0; i < num_minerals; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "mineral_reaction_rate[%s]", driver->chem_metadata.mineral_names.data[i]);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < num_cells; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_aux_output[j].mineral_reaction_rate.data[i];
  }
  for (int i = 0; i < num_primary; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "primary_free_ion_concentration[%s]", driver->chem_metadata.primary_names.data[i]);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < num_cells; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_aux_output[j].primary_free_ion_concentration.data[i];
  }
  for (int i = 0; i < num_primary; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "primary_activity_coeff[%s]", driver->chem_metadata.primary_names.data[i]);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < num_cells; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_aux_output[j].primary_activity_coeff.data[i];
  }
  for (int i = 0; i < num_aqueous_complexes; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "secondary_free_ion_concentration[%d]", i);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < num_cells; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_aux_output[j].secondary_free_ion_concentration.data[i];
  }
  for (int i = 0; i < num_aqueous_complexes; ++i, ++counter)
  {
    char var_name[1024];
    snprintf(var_name, 1023, "secondary_activity_coeff[%d]", i);
    var_names->data[counter] = AlquimiaStringDup(var_name);
    for (int j = 0; j < num_cells; ++j)
      var_data->data[num_vars*j + counter] = driver->chem_aux_output[j].secondary_activity_coeff.data[i];
  }
}

