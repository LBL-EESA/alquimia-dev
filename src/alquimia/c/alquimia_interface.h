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

#ifndef ALQUIMIA_INTERFACE_H_
#define ALQUIMIA_INTERFACE_H_

/*******************************************************************************
 **
 ** C implementation of the alquimia interface.
 **
 ******************************************************************************/

#include "alquimia_containers.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  /* NOTE(bja): AlquimiaData is a convenience container, and not a data
     structure that should be passed between the driver and
     engine. Conditions are not included here because they are short
     lived data structures. Need one data object per thread. */
  struct AlquimiaData {
    void* engine_state;
    struct AlquimiaSizes sizes;
    struct AlquimiaEngineFunctionality functionality;
    struct AlquimiaState state;
    struct AlquimiaProperties properties;
    struct AlquimiaAuxiliaryData aux_data;
    struct AlquimiaProblemMetaData meta_data;
    struct AlquimiaAuxiliaryOutputData aux_output;
  };

  /* NOTE(bja): The alquimia interface should contain nothing but
     function pointers. Inorder to thread chemistry, we need just one
     interface object. */ 
  struct AlquimiaInterface {
    /* read data files/structures, initialize memory, basis management
       (includes reading database, swapping basis, etc.) */
    void (*Setup)(
        const char* input_filename,
        bool* hands_off,
        void* pft_engine_state,
        struct AlquimiaSizes* sizes,
        struct AlquimiaEngineFunctionality* functionality,
        struct AlquimiaEngineStatus* status);

    /* gracefully shutdown the engine, cleanup memory */
    void (*Shutdown)(
      void* pft_engine_state,
      struct AlquimiaEngineStatus* status);

    /* constrain processing for boundary/initial constraints. Called
       once for each IC/BC. */
    void (*ProcessCondition)(
        void* pft_engine_state,
        struct AlquimiaGeochemicalCondition* condition,
        struct AlquimiaProperties* props,
        struct AlquimiaState* state,
        struct AlquimiaAuxiliaryData* aux_data,
        struct AlquimiaEngineStatus* status);

    /* take one (or more?) reaction steps in operator split mode */
    void (*ReactionStepOperatorSplit)(
        void* pft_engine_state,
        double* delta_t,
        struct AlquimiaProperties* props,
        struct AlquimiaState* state,
        struct AlquimiaAuxiliaryData* aux_data,
        struct AlquimiaEngineStatus* status);
    
    /* Access to user selected geochemical data for output, i.e. pH, 
       mineral SI, reaction rates */
    void (*GetAuxiliaryOutput)(
        void* pft_engine_state,
        struct AlquimiaProperties* props,
        struct AlquimiaState* state,
        struct AlquimiaAuxiliaryData* aux_data,
        struct AlquimiaAuxiliaryOutputData* aux_out,
        struct AlquimiaEngineStatus* status);
    
    void (*GetProblemMetaData)(
        void* pft_engine_state,
        struct AlquimiaProblemMetaData* meta_data,
        struct AlquimiaEngineStatus* status);
  };


  void CreateAlquimiaInterface(const char* const engine_name,
                               struct AlquimiaInterface* interface,
                               struct AlquimiaEngineStatus* status);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* ALQUIMIA_INTERFACE_H_ */
