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


#include "alquimia/alquimia_constants.h"

/* String lengths */
const int kAlquimiaMaxStringLength = 512;
const int kAlquimiaMaxWordLength = 32;

/* Geochemistry Engine Strings */
const char* kAlquimiaStringPFloTran = "PFloTran";
const char* kAlquimiaStringCrunchFlow = "CrunchFlow";
const char* kAlquimiaStringTotal = "total_aqueous";
const char* kAlquimiaStringTotalSorbed = "total_sorbed";
const char* kAlquimiaStringTotalAqueousPlusSorbed = "total_aqueous_plus_sorbed";
const char* kAlquimiaStringFree = "free";
const char* kAlquimiaStringPH = "pH";
const char* kAlquimiaStringMineral = "mineral";
const char* kAlquimiaStringGas = "gas";
const char* kAlquimiaStringCharge = "charge";

/* Error Codes */
const int kAlquimiaNoError = 0;
const int kAlquimiaErrorInvalidEngine = 1;
const int kAlquimiaErrorUnknownConstraintName = 2;
const int kAlquimiaErrorUnsupportedFunctionality = 3;
const int kAlquimiaErrorEngineIntegrity = 4577;

