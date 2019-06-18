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

#ifndef _GNU_SOURCE
#define _GNU_SOURCE // for older GNU compilers
#endif

#ifndef __USE_GNU
#define __USE_GNU   // for newer GNU compilers
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include <fenv.h>
#include "fpe.h"

#ifdef __APPLE__
#include <xmmintrin.h>
#else

#endif

void EnableFPE()
{
#ifdef __APPLE__
  // Catch all the interesting ones.
  _MM_SET_EXCEPTION_MASK(_MM_MASK_INEXACT | _MM_MASK_UNDERFLOW);
#else

  // Clear existing exceptions.
  feclearexcept(FE_ALL_EXCEPT);

  int flags = FE_DIVBYZERO |
              FE_INVALID   |
//              FE_UNDERFLOW |
              FE_OVERFLOW  ;

  // Enable only those we've selected.
  feenableexcept(flags);
#endif
}

void DisableFPE()
{
  fesetenv(FE_DFL_ENV);
}

// The following suspend/restore mechanism uses standard C99 machinery.
static fenv_t fpe_env;

void SuspendFPE()
{
  // Hold exceptions till further notice.
  feholdexcept(&fpe_env);
}

void ResumeFPE()
{
  // Clear all exception flags.
  feclearexcept(FE_ALL_EXCEPT);

  // Now restore the previous floating point environment.
  fesetenv(&fpe_env);
}

#ifdef __cplusplus
}
#endif
