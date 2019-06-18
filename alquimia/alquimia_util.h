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

#ifndef ALQUIMIA_C_UTIL_H_
#define ALQUIMIA_C_UTIL_H_

#include "alquimia/alquimia_containers.h"
#include "alquimia/alquimia_interface.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  // String manipulation functions.
  char* AlquimiaStringDup(const char* str);
  bool AlquimiaCaseInsensitiveStringCompare(const char* const str1,
                                            const char* const str2);

  void AlquimiaFindIndexFromName(const char* const name,
                                 const AlquimiaVectorString* const names,
                                 int* index);

  // Functions for copying Alquimia containers.
  void CopyAlquimiaVectorDouble(const AlquimiaVectorDouble* const source,
                                AlquimiaVectorDouble* destination);
  void CopyAlquimiaVectorInt(const AlquimiaVectorInt* const source,
                             AlquimiaVectorInt* destination);
  void CopyAlquimiaVectorString(const AlquimiaVectorString* const source,
                                AlquimiaVectorString* destination);
  void CopyAlquimiaSizes(const AlquimiaSizes* const source, 
                         AlquimiaSizes* destination);
  void CopyAlquimiaProblemMetaData(const AlquimiaProblemMetaData* const source, 
                                   AlquimiaProblemMetaData* destination);
  void CopyAlquimiaProperties(const AlquimiaProperties* const source, 
                              AlquimiaProperties* destination);
  void CopyAlquimiaEngineFunctionality(const AlquimiaEngineFunctionality* const source, 
                                       AlquimiaEngineFunctionality* destination);
  void CopyAlquimiaState(const AlquimiaState* const source, 
                         AlquimiaState* destination);
  void CopyAlquimiaAuxiliaryData(const AlquimiaAuxiliaryData* const source, 
                                 AlquimiaAuxiliaryData* destination);
  void CopyAlquimiaAuxiliaryOutputData(const AlquimiaAuxiliaryOutputData* const source, 
                                       AlquimiaAuxiliaryOutputData* destination);
  void CopyAlquimiaGeochemicalCondition(const AlquimiaGeochemicalCondition* const source, 
                                        AlquimiaGeochemicalCondition* destination);
  void CopyAlquimiaGeochemicalConditionVector(const AlquimiaGeochemicalConditionVector* source, 
                                              AlquimiaGeochemicalConditionVector* destination);
  void CopyAlquimiaAqueousConstraint(const AlquimiaAqueousConstraint* const source, 
                                     AlquimiaAqueousConstraint* destination);
  void CopyAlquimiaAqueousConstraintVector(const AlquimiaAqueousConstraintVector* const source, 
                                           AlquimiaAqueousConstraintVector* destination);
  void CopyAlquimiaMineralConstraint(const AlquimiaMineralConstraint* const source, 
                                     AlquimiaMineralConstraint* destination);
  void CopyAlquimiaMineralConstraintVector(const AlquimiaMineralConstraintVector* const source, 
                                           AlquimiaMineralConstraintVector* destination);

  // Functions for resizing Alquimia vectors.
  void ResizeAlquimiaVectorDouble(AlquimiaVectorDouble* vec, int new_size);
  void ResizeAlquimiaVectorInt(AlquimiaVectorInt* vec, int new_size);
  void ResizeAlquimiaVectorString(AlquimiaVectorString* vec, int new_size);
  void ResizeAlquimiaGeochemicalConditionVector(AlquimiaGeochemicalConditionVector* vec, int new_size);
  void ResizeAlquimiaAqueousConstraintVector(AlquimiaAqueousConstraintVector* vec, int new_size);
  void ResizeAlquimiaMineralConstraintVector(AlquimiaMineralConstraintVector* vec, int new_size);

  // The following functions write data to the given FILE.
  void PrintAlquimiaVectorDouble(const char* const name,
                                 const AlquimiaVectorDouble* const vector,
                                 FILE* file);
  void PrintAlquimiaVectorInt(const char* const name,
                              const AlquimiaVectorInt* const vector,
                              FILE* file);
  void PrintAlquimiaVectorString(const char* const name,
                                 const AlquimiaVectorString* const vector,
                                 FILE* file);

  void PrintAlquimiaData(const AlquimiaData* const data, FILE* file);
  void PrintAlquimiaSizes(const AlquimiaSizes* const sizes, FILE* file);
  void PrintAlquimiaProblemMetaData(const AlquimiaProblemMetaData* const meta_data, FILE* file);
  void PrintAlquimiaProperties(const AlquimiaProperties* const prop, FILE* file);
  void PrintAlquimiaEngineFunctionality(const AlquimiaEngineFunctionality* const functionality, FILE* file);
  void PrintAlquimiaState(const AlquimiaState* const state, FILE* file);
  void PrintAlquimiaAuxiliaryData(const AlquimiaAuxiliaryData* const aux_data, FILE* file);
  void PrintAlquimiaAuxiliaryOutputData(const AlquimiaAuxiliaryOutputData* const aux_output, FILE* file);
  void PrintAlquimiaGeochemicalConditionVector(const AlquimiaGeochemicalConditionVector* condition_list, FILE* file);
  void PrintAlquimiaGeochemicalCondition(const AlquimiaGeochemicalCondition* const condition, FILE* file);
  void PrintAlquimiaAqueousConstraint(const AlquimiaAqueousConstraint* const constraint, FILE* file);
  void PrintAlquimiaMineralConstraint(const AlquimiaMineralConstraint* const constraint, FILE* file);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* ALQUIMIA_C_UTIL_H_ */
