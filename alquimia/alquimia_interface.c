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
** 
** Authors: Benjamin Andre <bandre@lbl.gov>
*/

#include "alquimia/alquimia_interface.h"

#include "alquimia/pflotran_alquimia_interface.h"
#include "alquimia/crunch_alquimia_interface.h"

#include "alquimia/alquimia_util.h"
#include "alquimia/alquimia_constants.h"

void CreateAlquimiaInterface(const char* const engine_name,
                             AlquimiaInterface* interface,
                             AlquimiaEngineStatus* status) {
  interface->Setup = NULL;
  interface->Shutdown = NULL;
  interface->ProcessCondition = NULL;
  interface->ReactionStepOperatorSplit = NULL;
  interface->GetAuxiliaryOutput = NULL;
  interface->GetProblemMetaData = NULL;

  if (AlquimiaCaseInsensitiveStringCompare(engine_name,
                                           kAlquimiaStringPFloTran)) {
#if ALQUIMIA_HAVE_PFLOTRAN
    interface->Setup = &pflotran_alquimia_setup;
    interface->Shutdown = &pflotran_alquimia_shutdown;
    interface->ProcessCondition = &pflotran_alquimia_processcondition;
    interface->ReactionStepOperatorSplit = &pflotran_alquimia_reactionstepoperatorsplit;
    interface->GetAuxiliaryOutput = &pflotran_alquimia_getauxiliaryoutput;
    interface->GetProblemMetaData = &pflotran_alquimia_getproblemmetadata;
    status->error = kAlquimiaNoError;
    snprintf(status->message, kAlquimiaMaxStringLength,
             "CreateAlquimiaInterface() : successfully created PFloTran interface.\n");
#else
    status->error = kAlquimiaErrorInvalidEngine;
    snprintf(status->message, kAlquimiaMaxStringLength,
             "\nERROR : CreateAlquimiaInterface() : PFloTran interface requested, but alquimia was not compiled with PFloTran!\n");
#endif

  } else if (AlquimiaCaseInsensitiveStringCompare(engine_name,
                                                  kAlquimiaStringCrunchFlow)) {
#if ALQUIMIA_HAVE_CRUNCHFLOW
    /* interface->Setup = ...; */
    interface->Setup = &crunch_alquimia_setup;
    interface->Shutdown = &crunch_alquimia_shutdown;
    interface->ProcessCondition = &crunch_alquimia_processcondition;
    interface->ReactionStepOperatorSplit = &crunch_alquimia_reactionstepoperatorsplit;
    interface->GetAuxiliaryOutput = &crunch_alquimia_getauxiliaryoutput;
    interface->GetProblemMetaData = &crunch_alquimia_getproblemmetadata;
    status->error = kAlquimiaNoError;
    snprintf(status->message, kAlquimiaMaxStringLength,
             "CreateAlquimiaInterface() : successfully created CrunchFlow interface.\n");
#else
    status->error = kAlquimiaErrorInvalidEngine;
    snprintf(status->message, kAlquimiaMaxStringLength,
             "\nERROR : CreateAlquimiaInterface() : CrunchFlow interface requested, but alquimia was not compiled with CrunchFlow!\n");
#endif

  } else {
    status->error = kAlquimiaErrorInvalidEngine;
    snprintf(status->message, kAlquimiaMaxStringLength,
             "\nERROR : CreateAlquimiaInterface() : Invalid interface name '%s'.\n  Valid names are:\n    '%s'\n    '%s'\n",
             engine_name, kAlquimiaStringPFloTran, kAlquimiaStringCrunchFlow);
  }

}  /* end CreateAlquimiaInterface() */
