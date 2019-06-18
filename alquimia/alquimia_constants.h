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

#ifndef ALQUIMIA_CONSTANTS_H_
#define ALQUIMIA_CONSTANTS_H_

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* String lengths */
extern const int kAlquimiaMaxStringLength;
extern const int kAlquimiaMaxWordLength;

/* Geochemistry Engine Strings */
extern const char* kAlquimiaStringPFloTran;
extern const char* kAlquimiaStringCrunchFlow;
extern const char* kAlquimiaStringTotalAqueous;
extern const char* kAlquimiaStringTotalSorbed;
extern const char* kAlquimiaStringTotalAqueousPlusSorbed;
extern const char* kAlquimiaStringFree;
extern const char* kAlquimiaStringPH;
extern const char* kAlquimiaStringMineral;
extern const char* kAlquimiaStringGas;
extern const char* kAlquimiaStringCharge;

/* Error Codes */
extern const int kAlquimiaNoError;
extern const int kAlquimiaErrorInvalidEngine;
extern const int kAlquimiaErrorUnknownConstraintName;
extern const int kAlquimiaErrorUnsupportedFunctionality;
extern const int kAlquimiaErrorEngineIntegrity;

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif     /* ALQUIMIA_CONSTANTS_H_ */
