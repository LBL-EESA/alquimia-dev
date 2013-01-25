/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_C_UTIL_H_
#define ALQUIMIA_C_UTIL_H_

#include "alquimia_containers.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  bool AlquimiaCaseInsensitiveStringCompare(const char* const str1,
                                            const char* const str2);

  void PrintAlquimiaDoubleArray(const char* name, const int size, const double* array);
  void PrintAlquimiaIntArray(const char* name, const int size, const int* array);
  void PrintAlquimiaStringArray(const char* name, const int size, char** array);

  void PrintAlquimiaSizes(const struct AlquimiaSizes* sizes);
  void PrintAlquimiaMetaData(const struct AlquimiaMetaData* meta_data);
  void PrintAlquimiaState(const struct AlquimiaState* state);
  void PrintAlquimiaAuxiliaryData(const struct AlquimiaAuxiliaryData* aux_data);
  void PrintAlquimiaAuxiliaryOutputData(
      const struct AlquimiaAuxiliaryOutputData* aux_output);
  void PrintAlquimiaGeochemicalConditionList(
      const struct AlquimiaGeochemicalConditionList* condition_list);
  void PrintAlquimiaGeochemicalCondition(
      const struct AlquimiaGeochemicalCondition* condition);
  void PrintAlquimiaGeochemicalConstraint(
      const struct AlquimiaGeochemicalConstraint* constraint);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* ALQUIMIA_C_UTIL_H_ */
