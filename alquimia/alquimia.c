/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
** Alquimia Copyright (c) 2013-2016, The Regents of the University of California, 
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

#include <stdarg.h>
#include "alquimia/alquimia.h"

// Error handler.
static alquimia_error_handler_function error_handler = NULL;

// Default error handler.
static void default_error_handler(const char* message)
{
  printf("alquimia: Fatal error: %s\n", message);
  exit(-1);
}

void alquimia_error(const char* message, ...)
{
  // Set the default error handler if no handler is set.
  if (error_handler == NULL)
    error_handler = default_error_handler;

  // Extract the variadic arguments and splat them into a string.
  char err[1024];
  va_list argp;
  va_start(argp, message);
  vsnprintf(err, 1023, message, argp);
  va_end(argp);

  // Call the handler.
  error_handler(err);
}

void alquimia_set_error_handler(alquimia_error_handler_function handler)
{
  error_handler = handler;
}

// Abort function.
static alquimia_abort_function abort_function = NULL;

// Default abort function.
static void default_abort_function(const char* message)
{
  printf("alquimia: aborting: %s\n", message);
  abort();
}

void alquimia_abort(const char* message, ...)
{
  // Set the default abort function if none is set.
  if (abort_function == NULL)
    abort_function = default_abort_function;

  // Extract the variadic arguments and splat them into a string.
  char err[1024];
  va_list argp;
  va_start(argp, message);
  vsnprintf(err, 1023, message, argp);
  va_end(argp);

  // Call the function.
  abort_function(err);
}

void alquimia_set_abort_function(alquimia_abort_function function)
{
  abort_function = function;
}

void alquimia_version_fprintf(const char* exe_name, FILE* stream)
{
  if (stream == NULL) return;
    fprintf(stream, "%s v%s\n", exe_name, ALQUIMIA_VERSION);
}


