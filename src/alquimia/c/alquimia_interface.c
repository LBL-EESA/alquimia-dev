/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "alquimia_interface.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "pflotran_alquimia_interface.h"

#include "alquimia_containers.h"
#include "alquimia_constants.h"


void CreateAlquimiaInterface(const char* engine_name,
                             struct AlquimiaInterface* interface,
                             struct AlquimiaEngineStatus_C* status) {

  if (strncmp(engine_name, kPFloTran, strlen(kPFloTran)) == 0) {
#ifdef HAVE_PFLOTRAN
    interface->Setup = &pflotran_alquimia_setup;
    interface->Shutdown = &pflotran_alquimia_shutdown;
    interface->ProcessCondition = &pflotran_alquimia_processcondition;
    interface->ReactionStepOperatorSplit = &pflotran_alquimia_reactionstepoperatorsplit;
    interface->GetAuxiliaryOutput = &pflotran_alquimia_getauxiliaryoutput;
    interface->GetEngineMetaData = &pflotran_alquimia_getenginemetadata;
    interface->GetPrimaryNameFromIndex = &pflotran_alquimia_getprimarynamefromindex;
    status->error = ALQUIMIA_NO_ERROR;
    snprintf(status->message, ALQUIMIA_MAX_STRING_LENGTH,
             "CreateAlquimiaInterface() : successfully created PFloTran interface.\n");
#else
    status->error = ALQUIMIA_ERROR_INVALID_ENGINE;
    snprintf(status->message, ALQUIMIA_MAX_STRING_LENGTH,
             "\nERROR : CreateAlquimiaInterface() : PFloTran interface requested, but alquimia was not compiled with PFloTran!\n");
#endif

  } else if (strncmp(engine_name, kCrunchFlow, strlen(kCrunchFlow)) == 0) {
#ifdef HAVE_CRUNCH
    //interface->Setup = ...;
#else
      status->error = ALQUIMIA_ERROR_INVALID_ENGINE;
      snprintf(status->message, ALQUIMIA_MAX_STRING_LENGTH,
            "\nERROR : CreateAlquimiaInterface() : CrunchFlow interface requested, but alquimia was not compiled with CrunchFlow!\n");
#endif

  } else {
        snprintf(status->message, ALQUIMIA_MAX_STRING_LENGTH,
            "\nERROR : CreateAlquimiaInterface() : Invalid interface name '%s'.\n  Valid names are:\n    '%s'\n    '%s'\n",
            engine_name, kPFloTran, kCrunchFlow);
  }
  
}  // end CreateAlquimiaInterface()
