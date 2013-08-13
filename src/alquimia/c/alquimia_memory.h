/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
** Alquimia Copyright (c) 2013, The Regents of the University of California, 
** through Lawrence Berkeley National Laboratory (subject to receipt of any 
** required approvals from the U.S. Dept. of Energy).  All rights reserved.
** 
** If you have questions about your rights to use or distribute this software, 
** please contact Berkeley Lab's Technology Transfer and Intellectual Property 
** Management at TTD@lbl.gov referring to Alquimia (LBNL Ref. 2013-119).
** 
** NOTICE.  This software was developed under funding from the U.S. Department 
** of Energy.  As such, the U.S. Government has been granted for itself and 
** others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide 
** license in the Software to reproduce, prepare derivative works, and perform 
** publicly and display publicly.  Beginning five (5) years after the date 
** permission to assert copyright is obtained from the U.S. Department of Energy, 
** and subject to any subsequent five (5) year renewals, the U.S. Government is 
** granted for itself and others acting on its behalf a paid-up, nonexclusive, 
** irrevocable, worldwide license in the Software to reproduce, prepare derivative
** works, distribute copies to the public, perform publicly and display publicly, 
** and to permit others to do so.
** 
** Authors: Benjamin Andre <bandre@lbl.gov>
*/

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
  void AllocateAlquimiaState(const struct AlquimiaSizes* const sizes,
                             struct AlquimiaState* state);

  void FreeAlquimiaState(struct AlquimiaState* state);

  /* Auxiliary Data */ 
  void AllocateAlquimiaAuxiliaryData(const struct AlquimiaSizes* const sizes,
                                     struct AlquimiaAuxiliaryData* aux_data);
  void FreeAlquimiaAuxiliaryData(struct AlquimiaAuxiliaryData* aux_data);

  /* material properties */
  void AllocateAlquimiaMaterialProperties(
      const struct AlquimiaSizes* const sizes,
      struct AlquimiaMaterialProperties* material_props);
  void FreeAlquimiaMaterialProperties(
      struct AlquimiaMaterialProperties* material_props);

  /* Problem Meta Data */
  void AllocateAlquimiaProblemMetaData(const struct AlquimiaSizes* const sizes,
                                       struct AlquimiaProblemMetaData* meta_data);

  void FreeAlquimiaProblemMetaData(struct AlquimiaProblemMetaData* metda_data);

  /* Status */
  void AllocateAlquimiaEngineStatus(struct AlquimiaEngineStatus* status);

  void FreeAlquimiaEngineStatus(struct AlquimiaEngineStatus* status);

  /* Auxiliary Output Data */ 
  void AllocateAlquimiaAuxiliaryOutputData(
      const struct AlquimiaSizes* const sizes,
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
