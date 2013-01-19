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

  void pflotran_alquimia_setup(
      const char* input_filename,
      void* pft_engine_state,
      struct AlquimiaSizes_C* sizes,
      struct AlquimiaEngineStatus_C* status) ;
  void pflotran_alquimia_shutdown(
      void* pft_engine_state,
      struct AlquimiaEngineStatus_C* status);
  void pflotran_alquimia_processcondition(
      void* pft_engine_state,
      struct AlquimiaGeochemicalCondition_C* condition,
      struct AlquimiaMaterialProperties_C* material_props,
      struct AlquimiaState_C* state,
      struct AlquimiaAuxiliaryData_C* aux_data,
      struct AlquimiaEngineStatus_C* status);
  void pflotran_alquimia_reactionstepoperatorsplit(
      void* pft_engine_state,
      double* delta_t,
      struct AlquimiaMaterialProperties_C* material_properties,
      struct AlquimiaState_C* state,
      struct AlquimiaAuxiliaryData_C* aux_data,
      struct AlquimiaEngineStatus_C* status);
  void pflotran_alquimia_getauxiliaryoutput(
      void* pft_engine_state,
      struct AlquimiaEngineStatus_C* status);
  void pflotran_alquimia_getenginemetadata(
      void* pft_engine_state,
      struct AlquimiaSizes_C* sizes,
      struct AlquimiaMetaData_C* metadata,
      struct AlquimiaEngineStatus_C* status);
  void pflotran_alquimia_getprimarynamefromindex(
      void* pft_engine_state,
      int* primary_index,
      char* primary_name,
      struct AlquimiaEngineStatus_C* status);
#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif  // ALQUIMIA_PFLOTRAN_INTERFACE_H_
