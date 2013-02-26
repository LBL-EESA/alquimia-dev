/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_PFLOTRAN_INTERFACE_H_
#define ALQUIMIA_PFLOTRAN_INTERFACE_H_

/*******************************************************************************
 **
 ** C declarations of the pflotran alquimia interface
 **
 ******************************************************************************/

#include "alquimia_interface.h"

#include "alquimia_containers.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef HAVE_PFLOTRAN
  void pflotran_alquimia_setup(
      const char* input_filename,
      void* pft_engine_state,
      struct AlquimiaSizes* sizes,
      struct AlquimiaEngineFunctionality* functionality,
      struct AlquimiaEngineStatus* status);
  void pflotran_alquimia_shutdown(
      void* pft_engine_state,
      struct AlquimiaEngineStatus* status);
  void pflotran_alquimia_processcondition(
      void* pft_engine_state,
      struct AlquimiaGeochemicalCondition* condition,
      struct AlquimiaMaterialProperties* material_props,
      struct AlquimiaState* state,
      struct AlquimiaAuxiliaryData* aux_data,
      struct AlquimiaEngineStatus* status);
  void pflotran_alquimia_reactionstepoperatorsplit(
      void* pft_engine_state,
      double* delta_t,
      struct AlquimiaMaterialProperties* material_properties,
      struct AlquimiaState* state,
      struct AlquimiaAuxiliaryData* aux_data,
      struct AlquimiaEngineStatus* status);
  void pflotran_alquimia_getauxiliaryoutput(
      void* pft_engine_state,
      struct AlquimiaMaterialProperties* material_properties,
      struct AlquimiaState* state,
      struct AlquimiaAuxiliaryData* aux_data,
      struct AlquimiaAuxiliaryOutputData* aux_out,
      struct AlquimiaEngineStatus* status);
  void pflotran_alquimia_getproblemmetadata(
      void* pft_engine_state,
      struct AlquimiaProblemMetaData* meta_data,
      struct AlquimiaEngineStatus* status);
#endif

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif  // ALQUIMIA_PFLOTRAN_INTERFACE_H_
