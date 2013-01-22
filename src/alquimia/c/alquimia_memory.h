/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_C_MEMORY_H_
#define ALQUIMIA_C_MEMORY_H_

#include "alquimia_containers.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
  
  void AllocateDoubleArray(const int size, int *array_size, double** array);
  void AllocateIntArray(const int size, int *array_size, int** array);

  /* Alquimia Interface */
  void AllocateAlquimiaInterface(struct AlquimiaInterface* interface);

  void FreeAlquimiaInterface(struct AlquimiaInterface* interface);

  /* State */
  void AllocateAlquimiaState(const struct AlquimiaSizes* sizes,
                             struct AlquimiaState* state);

  void FreeAlquimiaState(struct AlquimiaState* state);

  /* Auxiliary Data */ 
  void AllocateAlquimiaAuxiliaryData(const struct AlquimiaSizes* sizes,
                                     struct AlquimiaAuxiliaryData* aux_data);
  void FreeAlquimiaAuxiliaryData(struct AlquimiaAuxiliaryData* aux_data);

  /* material properties */
  void AllocateAlquimiaMaterialProperties(
      const struct AlquimiaSizes* sizes,
      struct AlquimiaMaterialProperties* material_props);
  void FreeAlquimiaMaterialProperties(
      struct AlquimiaMaterialProperties* material_props);

  /* MetaData */
  void AllocateAlquimiaMetaData(const struct AlquimiaSizes* sizes,
                                struct AlquimiaMetaData* meta_data);

  void FreeAlquimiaMetaData(struct AlquimiaMetaData* metda_data);

  /* Status */
  void AllocateAlquimiaEngineStatus(struct AlquimiaEngineStatus* status);

  void FreeAlquimiaEngineStatus(struct AlquimiaEngineStatus* status);

  /* Geochemical conditions/constraints */
  void AllocateAlquimiaGeochemicalConditionList(
      const int num_conditions,
      struct AlquimiaGeochemicalConditionList* condition_list);

  void AllocateAlquimiaGeochemicalCondition(
      const char* name, const int num_constraints,
      struct AlquimiaGeochemicalCondition* condition);

  void AllocateAlquimiaGeochemicalConstraint(
      struct AlquimiaGeochemicalConstraint* constraint);

  void FreeAlquimiaGeochemicalConditionList(
      struct AlquimiaGeochemicalConditionList* conditions);

  void FreeAlquimiaGeochemicalCondition(
      struct AlquimiaGeochemicalCondition* condition);

  void FreeAlquimiaGeochemicalConstraint(
      struct AlquimiaGeochemicalConstraint* constraint);


  /* Data */
  void AllocateAlquimiaData(struct AlquimiaData* data);
  void FreeAlquimiaData(struct AlquimiaData* data);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* ALQUIMIA_C_MEMORY_H_ */
