/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_C_MEMORY_H_
#define ALQUIMIA_C_MEMORY_H_

#include "alquimia_containers.h"

/* State */
void AllocateAlquimiaState(const struct AlquimiaSizes_C* sizes,
                           struct AlquimiaState_C* state);

void FreeAlquimiaState(struct AlquimiaState_C* state);

/* MetaData */
void AllocateAlquimiaMetaData(const struct AlquimiaSizes_C* sizes,
                              struct AlquimiaMetaData_C* meta_data);

void FreeAlquimiaMetaData(const struct AlquimiaSizes_C* sizes,
                          struct AlquimiaMetaData_C* metda_data);

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


#endif  /* ALQUIMIA_C_MEMORY_H_ */
