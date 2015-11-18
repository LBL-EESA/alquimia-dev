/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
** Alquimia Copyright (c) 2013, The Regents of the University of California, 
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

  bool AlquimiaCaseInsensitiveStringCompare(const char* const str1,
                                            const char* const str2);

  void AlquimiaFindIndexFromName(const char* const name,
                                   const struct AlquimiaVectorString* const names,
                                   int* index);

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
  void PrintAlquimiaProperties(
      const struct AlquimiaProperties* const prop);
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
