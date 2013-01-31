/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_STRINGS_H_
#define ALQUIMIA_STRINGS_H_

#include "stddef.h"

/* String lengths */
extern const size_t kAlquimiaMaxStringLength;
extern const size_t kAlquimiaMaxWordLength;

/* Geochemistry Engine Strings */
extern const char* kAlquimiaStringPFloTran;
extern const char* kAlquimiaStringCrunchFlow;
extern const char* kAlquimiaStringTotal;
extern const char* kAlquimiaStringTotalSorbed;
extern const char* kAlquimiaStringFree;
extern const char* kAlquimiaStringPH;
extern const char* kAlquimiaStringMineral;
extern const char* kAlquimiaStringGas;
extern const char* kAlquimiaStringCharge;

/* Error Codes */
extern const int kAlquimiaNoError;
extern const int kAlquimiaErrorInvalidEngine;
extern const int kAlquimiaErrorUnknownConstraintName;
extern const int kAlquimiaErrorEngineIntegrity;


#endif     /* ALQUIMIA_STRINGS_H_ */
