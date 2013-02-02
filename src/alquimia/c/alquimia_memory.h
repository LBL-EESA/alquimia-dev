/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_C_MEMORY_H_
#define ALQUIMIA_C_MEMORY_H_

#include "alquimia_interface.h"
#include "alquimia_containers.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
  
  /* Alquimia Vectors */
  void AllocateAlquimiaVectorDouble(const int size, struct AlquimiaVectorDouble* vector);
  void FreeAlquimiaVectorDouble(struct AlquimiaVectorDouble* vector);

  void AllocateAlquimiaVectorInt(const int size, struct AlquimiaVectorInt* vector);
  void FreeAlquimiaVectorInt(struct AlquimiaVectorInt* vector);

  void AllocateAlquimiaVectorString(const int size, struct AlquimiaVectorString* vector);
  void FreeAlquimiaVectorString(struct AlquimiaVectorString* vector);

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

  /* Problem Meta Data */
  void AllocateAlquimiaProblemMetaData(const struct AlquimiaSizes* sizes,
                                       struct AlquimiaProblemMetaData* meta_data);

  void FreeAlquimiaProblemMetaData(struct AlquimiaProblemMetaData* metda_data);

  /* Status */
  void AllocateAlquimiaEngineStatus(struct AlquimiaEngineStatus* status);

  void FreeAlquimiaEngineStatus(struct AlquimiaEngineStatus* status);

  /* Auxiliary Output Data */ 
  void AllocateAlquimiaAuxiliaryOutputData(
      const struct AlquimiaSizes* sizes,
      struct AlquimiaAuxiliaryOutputData* aux_output);
  void FreeAlquimiaAuxiliaryOutputData(
      struct AlquimiaAuxiliaryOutputData* aux_output);

  /* Geochemical conditions/constraints */
  void AllocateAlquimiaGeochemicalConditionVector(
      const int num_conditions,
      struct AlquimiaGeochemicalConditionVector* condition_list);
  void AllocateAlquimiaGeochemicalCondition(
      const int size_name,
      const int num_aqueous_constraints, const int num_mineral_constraints,
      struct AlquimiaGeochemicalCondition* condition);
  void AllocateAlquimiaAqueousConstraint(
      struct AlquimiaAqueousConstraint* constraint);
  void AllocateAlquimiaMineralConstraint(
      struct AlquimiaMineralConstraint* constraint);

  void FreeAlquimiaGeochemicalConditionVector(
      struct AlquimiaGeochemicalConditionVector* condition_list);
  void FreeAlquimiaGeochemicalCondition(
      struct AlquimiaGeochemicalCondition* condition);
  void FreeAlquimiaAqueousConstraintVector(
      struct AlquimiaAqueousConstraintVector* vector);
  void FreeAlquimiaAqueousConstraint(
      struct AlquimiaAqueousConstraint* constraint);
  void FreeAlquimiaMineralConstraintVector(
      struct AlquimiaMineralConstraintVector* vector);
  void FreeAlquimiaMineralConstraint(
      struct AlquimiaMineralConstraint* constraint);


  /* Data */
  void AllocateAlquimiaData(struct AlquimiaData* data);
  void FreeAlquimiaData(struct AlquimiaData* data);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* ALQUIMIA_C_MEMORY_H_ */
