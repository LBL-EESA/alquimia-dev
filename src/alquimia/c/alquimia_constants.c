/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "alquimia_constants.h"

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
const int kAlquimiaErrorEngineIntegrity = 4577;

