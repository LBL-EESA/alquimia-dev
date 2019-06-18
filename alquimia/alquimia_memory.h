/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
** Alquimia Copyright (c) 2013-2016, The Regents of the University of California, 
** through Lawrence Berkeley National Laboratory (subject to receipt of any 
** required approvals from the U.S. Dept. of Energy).  All rights reserved.
** 
** Alquimia is available under a BSD license. See LICENSE.txt for more
** information.
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

#include "alquimia/alquimia_interface.h"
#include "alquimia/alquimia_containers.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
  
  /* Alquimia Vectors */
  void AllocateAlquimiaVectorDouble(const int size, AlquimiaVectorDouble* vector);
  void FreeAlquimiaVectorDouble(AlquimiaVectorDouble* vector);

  void AllocateAlquimiaVectorInt(const int size, AlquimiaVectorInt* vector);
  void FreeAlquimiaVectorInt(AlquimiaVectorInt* vector);

  void AllocateAlquimiaVectorString(const int size, AlquimiaVectorString* vector);
  void FreeAlquimiaVectorString(AlquimiaVectorString* vector);

  /* State */
  void AllocateAlquimiaState(const AlquimiaSizes* const sizes,
                             AlquimiaState* state);

  void FreeAlquimiaState(AlquimiaState* state);

  /* Auxiliary Data */ 
  void AllocateAlquimiaAuxiliaryData(const AlquimiaSizes* const sizes,
                                     AlquimiaAuxiliaryData* aux_data);
  void FreeAlquimiaAuxiliaryData(AlquimiaAuxiliaryData* aux_data);

  /* Properties */
  void AllocateAlquimiaProperties(const AlquimiaSizes* const sizes,
                                  AlquimiaProperties* props);
  void FreeAlquimiaProperties(AlquimiaProperties* props);

  /* Problem Meta Data */
  void AllocateAlquimiaProblemMetaData(const AlquimiaSizes* const sizes,
                                       AlquimiaProblemMetaData* meta_data);

  void FreeAlquimiaProblemMetaData(AlquimiaProblemMetaData* metda_data);

  /* Status */
  void AllocateAlquimiaEngineStatus(AlquimiaEngineStatus* status);

  void FreeAlquimiaEngineStatus(AlquimiaEngineStatus* status);

  /* Auxiliary Output Data */ 
  void AllocateAlquimiaAuxiliaryOutputData(const AlquimiaSizes* const sizes,
                                           AlquimiaAuxiliaryOutputData* aux_output);
  void FreeAlquimiaAuxiliaryOutputData(AlquimiaAuxiliaryOutputData* aux_output);

  /* Geochemical conditions/constraints */
  void AllocateAlquimiaGeochemicalConditionVector(const int num_conditions,
                                                  AlquimiaGeochemicalConditionVector* condition_list);
  void AllocateAlquimiaGeochemicalCondition(const int size_name,
                                            const int num_aqueous_constraints, 
                                            const int num_mineral_constraints,
                                            AlquimiaGeochemicalCondition* condition);
  void AllocateAlquimiaAqueousConstraintVector(const int num_constraints,
                                               AlquimiaAqueousConstraintVector* constraint_list);
  void AllocateAlquimiaAqueousConstraint(AlquimiaAqueousConstraint* constraint);
  void AllocateAlquimiaMineralConstraintVector(const int num_constraints,
                                               AlquimiaMineralConstraintVector* constraint_list);
  void AllocateAlquimiaMineralConstraint(AlquimiaMineralConstraint* constraint);

  void FreeAlquimiaGeochemicalConditionVector(AlquimiaGeochemicalConditionVector* condition_list);
  void FreeAlquimiaGeochemicalCondition(AlquimiaGeochemicalCondition* condition);
  void FreeAlquimiaAqueousConstraintVector(AlquimiaAqueousConstraintVector* vector);
  void FreeAlquimiaAqueousConstraint(AlquimiaAqueousConstraint* constraint);
  void FreeAlquimiaMineralConstraintVector(AlquimiaMineralConstraintVector* vector);
  void FreeAlquimiaMineralConstraint(AlquimiaMineralConstraint* constraint);


  /* Data */
  void AllocateAlquimiaData(AlquimiaData* data);
  void FreeAlquimiaData(AlquimiaData* data);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* ALQUIMIA_C_MEMORY_H_ */
