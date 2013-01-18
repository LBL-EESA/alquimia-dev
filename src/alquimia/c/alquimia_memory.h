/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_C_MEMORY_H_
#define ALQUIMIA_C_MEMORY_H_

#include "alquimia_containers.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
  
  /* State */
  void AllocateAlquimiaState(const struct AlquimiaSizes_C* sizes,
                             struct AlquimiaState_C* state);

  void FreeAlquimiaState(struct AlquimiaState_C* state);

  /* Auxiliary Data */ 
  void AllocateAlquimiaAuxiliaryData(const struct AlquimiaSizes_C* sizes,
                                     struct AlquimiaAuxiliaryData_C* aux_data);
  void FreeAlquimiaAuxiliaryData(struct AlquimiaAuxiliaryData_C* aux_data);

  /* material properties */
  void AllocateAlquimiaMaterialProperties(
      const struct AlquimiaSizes_C* sizes,
      struct AlquimiaMaterialProperties_C* material_props);
  void FreeAlquimiaMaterialProperties(
      struct AlquimiaMaterialProperties_C* material_props);

  /* MetaData */
  void AllocateAlquimiaMetaData(const struct AlquimiaSizes_C* sizes,
                                struct AlquimiaMetaData_C* meta_data);

  void FreeAlquimiaMetaData(const struct AlquimiaSizes_C* sizes,
                            struct AlquimiaMetaData_C* metda_data);

  /* Status */
  void AllocateAlquimiaEngineStatus(struct AlquimiaEngineStatus_C* status);

  void FreeAlquimiaEngineStatus(struct AlquimiaEngineStatus_C* status);

  /* Geochemical conditions/constraints */
  void AllocateAlquimiaGeochemicalConditionList(
      const int num_conditions,
      struct AlquimiaGeochemicalConditionList_C* condition_list);

  void AllocateAlquimiaGeochemicalCondition(
      const char* name, const int num_constraints,
      struct AlquimiaGeochemicalCondition_C* condition);

  void AllocateAlquimiaGeochemicalConstraint(
      struct AlquimiaGeochemicalConstraint_C* constraint);

  void FreeAlquimiaGeochemicalConditionList(
      struct AlquimiaGeochemicalConditionList_C* conditions);

  void FreeAlquimiaGeochemicalCondition(
      struct AlquimiaGeochemicalCondition_C* condition);

  void FreeAlquimiaGeochemicalConstraint(
      struct AlquimiaGeochemicalConstraint_C* constraint);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* ALQUIMIA_C_MEMORY_H_ */
