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

#ifndef ALQUIMIA_BATCH_CHEM_DRIVER_H_
#define ALQUIMIA_BATCH_CHEM_DRIVER_H_

#include "alquimia/alquimia_containers.h"

// This type stores the metadata for a reactive batch chemistgry simulation.
typedef struct BatchChemDriver BatchChemDriver;

// This is a ceiling on the number of primary/secondary species, complexes, etc.
#define BATCH_CHEM_INPUT_MAX 128

// This type holds simulation input information for the batch chemistry driver.
typedef struct
{
  // Problem description.
  char* description;

  // ---------------------
  // Simulation parameters
  // ---------------------

  // Do we use hands-off mode to read property data directly from the engine's
  // input file?
  bool hands_off;

  // Start/stop time, max number of steps.
  double t_min, t_max;
  int max_steps;

  // Specified timestep.
  double dt;

  // --------
  // Metadata
  // --------
  int num_isotherm_species, num_ion_exchange_sites, num_surface_sites;
  char* isotherm_species[BATCH_CHEM_INPUT_MAX];
  char* ion_exchange_sites[BATCH_CHEM_INPUT_MAX];
  char* surface_sites[BATCH_CHEM_INPUT_MAX];

  // ---------------
  // State variables
  // ---------------
  double water_density, temperature, porosity, aqueous_pressure;
  double surface_site_density[BATCH_CHEM_INPUT_MAX];
  double cation_exchange_capacity[BATCH_CHEM_INPUT_MAX];

  // -------------------
  // Material properties
  // -------------------
  double volume, saturation;
  double isotherm_kd[BATCH_CHEM_INPUT_MAX];
  double langmuir_b[BATCH_CHEM_INPUT_MAX];
  double freundlich_n[BATCH_CHEM_INPUT_MAX];

  // ----------------------
  // Geochemical condition.
  // ----------------------
  char* cond_name;

  // ---------------------
  // Chemistry engine info
  // ---------------------
  char* chemistry_engine;
  char* chemistry_input_file;

  // ------------------
  // Output information
  // ------------------
  bool verbose;
  char* output_file;
  char* output_type;

} BatchChemDriverInput;

// Parses an input file and produces a BatchChemDriverInput.
BatchChemDriverInput* BatchChemDriverInput_New(const char* input_file);

// Frees the given BatchChemDriverInput.
void BatchChemDriverInput_Free(BatchChemDriverInput* input);

// Creates a new BatchChemDriver object given input parameters. 
BatchChemDriver* BatchChemDriver_New(BatchChemDriverInput* input);

// Destroys the given BatchChemDriver object, freeing its resources.
void BatchChemDriver_Free(BatchChemDriver* driver);

// Takes a single step in the given batch simulation.
int BatchChemDriver_Run(BatchChemDriver* driver);

// Runs the batcbatch simulation as defined by input.
int BatchChemDriver_Run(BatchChemDriver* driver);

// Retrieves all solute data and auxiliary data corresponding to the 
// internal state of the batch chemistry model, placing the names of the variables
// into var_names, and placing the various data into components of the 
// multi-component vector var_data. These vectors must be destroyed by the 
// caller.
void BatchChemDriver_GetSoluteAndAuxData(BatchChemDriver* driver,
                                         double* time,
                                         AlquimiaVectorString* var_names,
                                         AlquimiaVectorDouble* var_data);

#endif
