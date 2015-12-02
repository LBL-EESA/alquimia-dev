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

#include "alquimia/alquimia.h"
#include "ini.h"
#include "TransportInput.h"

struct TransportInput
{
  TransportInputCoupling coupling;
  double t_min, t_max;
  int max_steps;
  double cfl;
  int order;

  double x_min, x_max;
  int num_cells;

  double T, ux;

  char chem_engine[FILENAME_MAX]; 
  char chem_input[FILENAME_MAX]; 
  char chem_ic[FILENAME_MAX]; 
  char chem_left_bc[FILENAME_MAX]; 
  char chem_right_bc[FILENAME_MAX]; 

  char output_file[FILENAME_MAX]; 
  char output_type[FILENAME_MAX]; 
};

TransportInput* TransportInput_New(const char* input_file)
{
  TransportInput* input = malloc(sizeof(TransportInput));
  return input;
}

void TransportInput_Free(TransportInput* input)
{
  free(input);
}

void TransportInput_GetSimParameters(TransportInput* input,
                                     TransportInputCoupling* coupling,
                                     double* t_min, 
                                     double* t_max,
                                     int* max_steps,
                                     double* cfl_factor,
                                     int* order_of_accuracy)
{
  *coupling = input->coupling;
  *t_min = input->t_min;
  *t_max = input->t_max;
  *max_steps = input->max_steps;
  *cfl_factor = input->cfl;
  *order_of_accuracy = input->order;
}

void TransportInput_GetDomain(TransportInput* input,
                              double* x_min, 
                              double* x_max, 
                              int* num_cells)
{
  *x_min = input->x_min;
  *x_max = input->x_max;
  *num_cells = input->num_cells;
}

void TransportInput_GetFlow(TransportInput* input,
                            double* temperature,
                            double* x_velocity)
{
  *temperature = input->T;
  *x_velocity = input->ux;
}

void TransportInput_GetChemistry(TransportInput* input,
                                 char** chemistry_engine,
                                 char** chemistry_input_file,
                                 char** chemical_ic_name,
                                 char** left_chemical_bc_name,
                                 char** right_chemical_bc_name)
{
  *chemistry_engine = input->chem_engine;
  *chemistry_input_file = input->chem_input;
  *chemical_ic_name = input->chem_ic;
  *left_chemical_bc_name = input->chem_left_bc;
  *right_chemical_bc_name = input->chem_right_bc;
}

void TransportInput_GetOutput(TransportInput* input,
                              char** filename,
                              char** output_type)
{
  *filename = input->output_file;
  *output_type = input->output_type;
}
