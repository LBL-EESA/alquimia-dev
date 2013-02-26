/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_C_UTIL_H_
#define ALQUIMIA_C_UTIL_H_

#include "alquimia_containers.h"
#include "alquimia_interface.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  bool AlquimiaCaseInsensitiveStringCompare(const char* const str1,
                                            const char* const str2);

  void AlquimiaEngineIndexFromName(const char* const name,
                                   const struct AlquimiaVectorString* const names,
                                   const struct AlquimiaVectorInt* const indices,
                                   int* index);

  void AlquimiaNameFromEngineIndex(const int index,
                                   const struct AlquimiaVectorString* const names,
                                   const struct AlquimiaVectorInt* const indices,
                                   char* name);

  void PrintAlquimiaVectorDouble(const char* const name,
                                 const struct AlquimiaVectorDouble* const vector);
  void PrintAlquimiaVectorInt(const char* const name,
                              const struct AlquimiaVectorInt* const vector);
  void PrintAlquimiaVectorString(const char* const name,
                                 const struct AlquimiaVectorString* const vector);

  void PrintAlquimiaData(const struct AlquimiaData* const data);
  void PrintAlquimiaSizes(const struct AlquimiaSizes* const sizes);
  void PrintAlquimiaProblemMetaData(
      const struct AlquimiaProblemMetaData* const meta_data);
  void PrintAlquimiaMaterialProperties(
      const struct AlquimiaMaterialProperties* const mat_prop);
  void PrintAlquimiaEngineFunctionality(
      const struct AlquimiaEngineFunctionality* const functionality);
  void PrintAlquimiaState(const struct AlquimiaState* const state);
  void PrintAlquimiaAuxiliaryData(const struct AlquimiaAuxiliaryData* const aux_data);
  void PrintAlquimiaAuxiliaryOutputData(
      const struct AlquimiaAuxiliaryOutputData* const aux_output);
  void PrintAlquimiaGeochemicalConditionVector(
      const struct AlquimiaGeochemicalConditionVector* condition_list);
  void PrintAlquimiaGeochemicalCondition(
      const struct AlquimiaGeochemicalCondition* const condition);
  void PrintAlquimiaAqueousConstraint(
      const struct AlquimiaAqueousConstraint* const constraint);
  void PrintAlquimiaMineralConstraint(
      const struct AlquimiaMineralConstraint* const constraint);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* ALQUIMIA_C_UTIL_H_ */
