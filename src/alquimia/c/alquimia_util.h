/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_C_UTIL_H_
#define ALQUIMIA_C_UTIL_H_

#include "alquimia_containers.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  bool AlquimiaCaseInsensitiveStringCompare(const char* const str1,
                                            const char* const str2);

  void PrintAlquimiaVectorDouble(const char* const name,
                                 const struct AlquimiaVectorDouble* const vector);
  void PrintAlquimiaVectorInt(const char* const name,
                              const struct AlquimiaVectorInt* const vector);
  void PrintAlquimiaVectorString(const char* const name,
                                 const struct AlquimiaVectorString* const vector);

  void PrintAlquimiaSizes(const struct AlquimiaSizes* const sizes);
  void PrintAlquimiaMetaData(const struct AlquimiaMetaData* const meta_data);
  void PrintAlquimiaState(const struct AlquimiaState* const state);
  void PrintAlquimiaAuxiliaryData(const struct AlquimiaAuxiliaryData* const aux_data);
  void PrintAlquimiaAuxiliaryOutputData(
      const struct AlquimiaAuxiliaryOutputData* const aux_output);
  void PrintAlquimiaGeochemicalConditionList(
      const struct AlquimiaGeochemicalConditionList* condition_list);
  void PrintAlquimiaGeochemicalCondition(
      const struct AlquimiaGeochemicalCondition* const condition);
  void PrintAlquimiaGeochemicalConstraint(
      const struct AlquimiaGeochemicalConstraint* const constraint);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* ALQUIMIA_C_UTIL_H_ */
