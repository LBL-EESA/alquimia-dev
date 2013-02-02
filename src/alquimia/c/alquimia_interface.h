/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
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
    struct AlquimiaMaterialProperties material_properties;
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
        struct AlquimiaMaterialProperties* material_props,
        struct AlquimiaState* state,
        struct AlquimiaAuxiliaryData* aux_data,
        struct AlquimiaEngineStatus* status);

    /* take one (or more?) reaction steps in operator split mode */
    void (*ReactionStepOperatorSplit)(
        void* pft_engine_state,
        double* delta_t,
        struct AlquimiaMaterialProperties* material_props,
        struct AlquimiaState* state,
        struct AlquimiaAuxiliaryData* aux_data,
        struct AlquimiaEngineStatus* status);
    
    /* Access to user selected geochemical data for output, i.e. pH, 
       mineral SI, reaction rates */
    void (*GetAuxiliaryOutput)(
        void* pft_engine_state,
        struct AlquimiaMaterialProperties* material_props,
        struct AlquimiaState* state,
        struct AlquimiaAuxiliaryData* aux_data,
        struct AlquimiaAuxiliaryOutputData* aux_out,
        struct AlquimiaEngineStatus* status);
    
    void (*GetProblemMetaData)(
        void* pft_engine_state,
        struct AlquimiaProblemMetaData* meta_data,
        struct AlquimiaEngineStatus* status);
  };


  void CreateAlquimiaInterface(const char* engine_name,
                               struct AlquimiaInterface* interface,
                               struct AlquimiaEngineStatus* status);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  // ALQUIMIA_INTERFACE_H_
