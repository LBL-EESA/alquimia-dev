/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_C_UTIL_H_
#define ALQUIMIA_C_UTIL_H_

#include "alquimia_containers.h"

void PrintAlquimiaSizes(const struct AlquimiaSizes_C* sizes);
void PrintAlquimiaMetaData(const struct AlquimiaSizes_C* sizes,
                           const struct AlquimiaMetaData_C* meta_data);
void PrintAlquimiaState(const struct AlquimiaSizes_C* sizes,
                        const struct AlquimiaState_C* state);


void PrintAlquimiaGeochemicalConditionList(
    const struct AlquimiaGeochemicalConditionList_C* condition_list);
void PrintAlquimiaGeochemicalCondition(
    const struct AlquimiaGeochemicalCondition_C* condition);
void PrintAlquimiaGeochemicalConstraint(
    const struct AlquimiaGeochemicalConstraint_C* constraint);

#endif  /* ALQUIMIA_C_UTIL_H_ */
