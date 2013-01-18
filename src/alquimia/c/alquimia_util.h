/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_C_UTIL_H_
#define ALQUIMIA_C_UTIL_H_

#include "alquimia_containers.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  void PrintAlquimiaSizes(const struct AlquimiaSizes_C* sizes);
  void PrintAlquimiaMetaData(const struct AlquimiaSizes_C* sizes,
                             const struct AlquimiaMetaData_C* meta_data);
  void PrintAlquimiaState(const struct AlquimiaSizes_C* sizes,
                          const struct AlquimiaState_C* state);
  void PrintAlquimiaAuxiliaryData(const struct AlquimiaSizes_C* sizes,
                                  const struct AlquimiaAuxiliaryData_C* aux_data);
  void PrintAlquimiaGeochemicalConditionList(
      const struct AlquimiaGeochemicalConditionList_C* condition_list);
  void PrintAlquimiaGeochemicalCondition(
      const struct AlquimiaGeochemicalCondition_C* condition);
  void PrintAlquimiaGeochemicalConstraint(
      const struct AlquimiaGeochemicalConstraint_C* constraint);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* ALQUIMIA_C_UTIL_H_ */
