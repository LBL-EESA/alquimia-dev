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

#ifndef ALQUIMIA_TRANSPORT_DRIVER_H_
#define ALQUIMIA_TRANSPORT_DRIVER_H_

#include "alquimia/alquimia_containers.h"

// This type stores the metadata for a reactive transport simulation.
typedef struct TransportDriver TransportDriver;

// Types of coupling for time integration.
typedef enum
{
  TRANSPORT_OPERATOR_SPLIT,
  TRANSPORT_GLOBAL_IMPLICIT
} TransportCoupling;
  
// This type holds simulation input information for the transport driver.
typedef struct
{
  // Problem description.
  char* description;

  // ---------------------
  // Simulation parameters
  // ---------------------
  TransportCoupling coupling;

  // Do we use hands-off mode to read property data directly from the engine's
  // input file?
  bool hands_off;

  // Start/stop time, max number of steps.
  double t_min, t_max;
  int max_steps;

  // Specified timestep.
  double dt;

  // CFL factor.
  double cfl_factor;

  // --------------------
  // Computational domain 
  // --------------------
  double x_min, x_max;
  int num_cells;

  // --------------
  // Flow variables
  // --------------
  double velocity, temperature;

  // -------------------
  // Material properties
  // -------------------
  double porosity, volume, saturation;

  // -------------------------------
  // Initial and boundary conditions
  // -------------------------------
  char* ic_name;
  char* left_bc_name;
  char* right_bc_name;

  // ---------------------
  // Chemistry engine info
  // ---------------------
  char* chemistry_engine;
  char* chemistry_input_file;

  // -----------------
  // Geochemistry data
  // -----------------
  double cation_exchange_capacity;
  double surface_site_density;

  // ------------------
  // Output information
  // ------------------
  bool verbose;
  char* output_file;
  char* output_type;

} TransportDriverInput;

// Parses an input file and produces a TransportDriverInput.
TransportDriverInput* TransportDriverInput_New(const char* input_file);

// Frees the given TransportDriverInput.
void TransportDriverInput_Free(TransportDriverInput* input);

// Creates a new TransportDriver object given input parameters. 
TransportDriver* TransportDriver_New(TransportDriverInput* input);

// Destroys the given TransportDriver object, freeing its resources.
void TransportDriver_Free(TransportDriver* driver);

// Runs the transport simulation as defined by input.
int TransportDriver_Run(TransportDriver* driver);

// Retrieves all solute data and auxiliary data corresponding to the 
// internal state of the transport model, placing the names of the variables
// into var_names, and placing the various data into components of the 
// multi-component vector var_data. These vectors must be destroyed by the 
// caller.
void TransportDriver_GetSoluteAndAuxData(TransportDriver* driver,
                                         double* time,
                                         AlquimiaVectorString* var_names,
                                         AlquimiaVectorDouble* var_data);

#endif
