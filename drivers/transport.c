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

#include "petsc.h"
#include "alquimia/alquimia_memory.h"
#include "alquimia/alquimia_util.h"
#include "TransportDriver.h"
#include "DriverOutput.h"

void Usage()
{
  printf("transport: usage:\n");
  printf("transport <input_file>\n\n");
  exit(0);
}

int main(int argc, char* argv[]) 
{
  if (argc == 1)
    Usage();

  // Initialize PETSc/MPI for command line options and engines that
  // require it.
  char help[] = "Alquimia advective, nondispersive reactive transport driver";
  PetscInitialize(&argc, &argv, (char*)0, help);
  PetscInitializeFortran();

  char input_file[FILENAME_MAX];
  strncpy(input_file, argv[1], FILENAME_MAX-1);

  // Parse the input file.
  TransportDriverInput* input = TransportDriverInput_New(input_file);
  if (input == NULL)
    alquimia_error("transport: error encountered reading input file '%s'.", input_file);

  // Set up output.
  DriverOutput* output = NULL;
  if (AlquimiaCaseInsensitiveStringCompare(input->output_type, "python"))
    output = PythonDriverOutput_New();
  else if (AlquimiaCaseInsensitiveStringCompare(input->output_type, "gnuplot"))
    output = GnuplotDriverOutput_New();

  // Create a TransportDriver from the parsed input.
  TransportDriver* transport = TransportDriver_New(input);

  // Run the simulation.
  int status = TransportDriver_Run(transport);

  // Get the solution out of the driver and write it out.
  if (output != NULL)
  {
    double final_time;
    AlquimiaVectorString var_names = {.size = 0};
    AlquimiaVectorDouble var_data = {.size = 0};
    TransportDriver_GetSoluteAndAuxData(transport, &final_time, &var_names, &var_data);
    DriverOutput_WriteMulticompVector(output, input->output_file, var_names, var_data);

    FreeAlquimiaVectorString(&var_names);
    FreeAlquimiaVectorDouble(&var_data);
  }

  // Clean up.
  TransportDriverInput_Free(input);
  TransportDriver_Free(transport);
  PetscInt petsc_error = PetscFinalize();
  if (status == EXIT_SUCCESS && petsc_error == 0) 
    printf("Success!\n");
  else 
    printf("Failed!\n");

  return status;
}  

