/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
** Alquimia Copyright (c) 2013, The Regents of the University of California, 
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

#ifndef ALQUIMIA_CRUNCH_INTERFACE_H_
#define ALQUIMIA_CRUNCH_INTERFACE_H_

/*******************************************************************************
 **
 ** C declarations of the crunchflow alquimia interface
 **
 ******************************************************************************/

#include "alquimia_interface.h"

#include "alquimia_containers.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef HAVE_CRUNCH
  void crunch_alquimia_setup(
      const char* input_filename,
      bool* hands_off,
      void* cf_engine_state,
      struct AlquimiaSizes* sizes,
      struct AlquimiaEngineFunctionality* functionality,
      struct AlquimiaEngineStatus* status);
  void crunch_alquimia_shutdown(
      void* cf_engine_state,
      struct AlquimiaEngineStatus* status);
  void crunch_alquimia_processcondition(
      void* cf_engine_state,
      struct AlquimiaGeochemicalCondition* condition,
      struct AlquimiaProperties* props,
      struct AlquimiaState* state,
      struct AlquimiaAuxiliaryData* aux_data,
      struct AlquimiaEngineStatus* status);
  void crunch_alquimia_reactionstepoperatorsplit(
      void* cf_engine_state,
      double* delta_t,
      struct AlquimiaProperties* properties,
      struct AlquimiaState* state,
      struct AlquimiaAuxiliaryData* aux_data,
      struct AlquimiaEngineStatus* status);
  void crunch_alquimia_getauxiliaryoutput(
      void* cf_engine_state,
      struct AlquimiaProperties* properties,
      struct AlquimiaState* state,
      struct AlquimiaAuxiliaryData* aux_data,
      struct AlquimiaAuxiliaryOutputData* aux_out,
      struct AlquimiaEngineStatus* status);
  void crunch_alquimia_getproblemmetadata(
      void* cf_engine_state,
      struct AlquimiaProblemMetaData* meta_data,
      struct AlquimiaEngineStatus* status);
#endif

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif  /* ALQUIMIA_CRUNCH_INTERFACE_H_ */
