/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "alquimia_interface.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "pflotran_alquimia_interface.h"

#include "alquimia_util.h"
#include "alquimia_containers.h"
#include "alquimia_constants.h"

void CreateAlquimiaInterface(const char* const engine_name,
                             struct AlquimiaInterface* interface,
                             struct AlquimiaEngineStatus* status) {
  interface->Setup = NULL;
  interface->Shutdown = NULL;
  interface->ProcessCondition = NULL;
  interface->ReactionStepOperatorSplit = NULL;
  interface->GetAuxiliaryOutput = NULL;
  interface->GetProblemMetaData = NULL;

  if (AlquimiaCaseInsensitiveStringCompare(engine_name,
                                           kAlquimiaStringPFloTran)) {
#ifdef HAVE_PFLOTRAN
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
#ifdef HAVE_CRUNCH
    //interface->Setup = ...;
    status->error = kAlquimiaErrorInvalidEngine;
    snprintf(status->message, kAlquimiaMaxStringLength,
             "\nERROR : CreateAlquimiaInterface() : CrunchFlow interface requested, but has not been implemented!\n");
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

}  // end CreateAlquimiaInterface()
