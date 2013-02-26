/* -*-  mode: c; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*******************************************************************************
 **
 **  Unit tests for alquimia C utilities
 **
 *******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <stdbool.h>

#include "alquimia_util.h"
#include "alquimia_memory.h"

#include "pflotran_alquimia_interface.h"

#include "alquimia_constants.h"
#include "alquimia_containers.h"

void test_AlquimiaCaseInsensitiveStringCompare(void) {
  char* str1;
  char* str2;
  
  str1 = (char*) calloc(kAlquimiaMaxStringLength, sizeof(char));
  str2 = (char*) calloc(kAlquimiaMaxStringLength, sizeof(char));
  
  strncpy(str1, "Hello", kAlquimiaMaxStringLength);
  strncpy(str2, "World", kAlquimiaMaxStringLength);

  assert(AlquimiaCaseInsensitiveStringCompare(str1, str2) == false);

  strncpy(str2, "Goodbye", kAlquimiaMaxStringLength);
  assert(AlquimiaCaseInsensitiveStringCompare(str1, str2) == false);

  strncpy(str2, "Hello", kAlquimiaMaxStringLength);
  assert(AlquimiaCaseInsensitiveStringCompare(str1, str2) == true);

  strncpy(str2, "hELLO", kAlquimiaMaxStringLength);
  assert(AlquimiaCaseInsensitiveStringCompare(str1, str2) == true);
}  // end test_AlquimiaCaseInsensitiveStringCompare

void test_AlquimiaVectors(void) {
  int size;
  struct AlquimiaVectorDouble dvector;
  struct AlquimiaVectorInt ivector;
  struct AlquimiaVectorString svector;

  size = -1;
  AllocateAlquimiaVectorDouble(size, &dvector);
  assert(dvector.size == 0);
  assert(dvector.data == NULL);

  size = 5;
  AllocateAlquimiaVectorDouble(size, &dvector);
  assert(dvector.size == size);
  assert(dvector.data != NULL);
  FreeAlquimiaVectorDouble(&dvector);

  size = -1;
  AllocateAlquimiaVectorInt(size, &ivector);
  assert(ivector.size == 0);
  assert(ivector.data == NULL);

  size = 5;
  AllocateAlquimiaVectorInt(size, &ivector);
  assert(ivector.size == size);
  assert(ivector.data != NULL);
  FreeAlquimiaVectorInt(&ivector);

  size = -1;
  AllocateAlquimiaVectorString(size, &svector);
  assert(svector.size == 0);
  assert(svector.data == NULL);

  size = 5;
  AllocateAlquimiaVectorString(size, &svector);
  assert(svector.size == size);
  assert(svector.data != NULL);
  FreeAlquimiaVectorString(&svector);
}  // end test_AlquimiaVectors()

void test_AlquimiaNameIndexMapping(void) {
  int i, id, size;
  char* name;
  struct AlquimiaVectorString names;

  name = (char*) calloc(kAlquimiaMaxStringLength, sizeof(char));
  size = 5;
  AllocateAlquimiaVectorString(size, &names);
  for (i = 0; i < size; ++i) {
    snprintf(names.data[i], kAlquimiaMaxStringLength, "name_%d", i+10);
  }
  //PrintAlquimiaVectorString("names", &names);

  strncpy(name, "foo", kAlquimiaMaxStringLength);
  AlquimiaFindIndexFromName(name, &names, &id);
  assert(id == -1);

  strncpy(name, "name_13", kAlquimiaMaxStringLength);
  AlquimiaFindIndexFromName(name, &names, &id);
  assert(id == 3);

  FreeAlquimiaVectorString(&names);
}  // end test_AlquimiaNameIndexMapping()

void test_CreateAlquimiaInterface(void) {
  struct AlquimiaEngineStatus status;
  struct AlquimiaInterface interface;
  char* name;

  AllocateAlquimiaEngineStatus(&status);
  name = (char*) calloc(kAlquimiaMaxStringLength, sizeof(char));

  strncpy(name, "junk", kAlquimiaMaxStringLength);
  CreateAlquimiaInterface(name, &interface, &status);
  assert(status.error == kAlquimiaErrorInvalidEngine);
  assert(interface.Setup == NULL);
  assert(interface.Shutdown == NULL);
  assert(interface.ProcessCondition == NULL);
  assert(interface.ReactionStepOperatorSplit == NULL);
  assert(interface.GetAuxiliaryOutput == NULL);
  assert(interface.GetProblemMetaData == NULL);

  strncpy(name, "pflotran", kAlquimiaMaxStringLength);
  CreateAlquimiaInterface(name, &interface, &status);
#ifdef HAVE_PFLOTRAN
  assert(status.error == kAlquimiaNoError);
  assert(interface.Setup == &pflotran_alquimia_setup);
  assert(interface.Shutdown == &pflotran_alquimia_shutdown);
  assert(interface.ProcessCondition == &pflotran_alquimia_processcondition);
  assert(interface.ReactionStepOperatorSplit == &pflotran_alquimia_reactionstepoperatorsplit);
  assert(interface.GetAuxiliaryOutput == &pflotran_alquimia_getauxiliaryoutput);
  assert(interface.GetProblemMetaData == &pflotran_alquimia_getproblemmetadata);
#else
  assert(status.error == kAlquimiaErrorInvalidEngine);
  assert(interface.Setup == NULL);
  assert(interface.Shutdown == NULL);
  assert(interface.ProcessCondition == NULL);
  assert(interface.ReactionStepOperatorSplit == NULL);
  assert(interface.GetAuxiliaryOutput == NULL);
  assert(interface.GetProblemMetaData == NULL);
#endif

  strncpy(name, "crunchflow", kAlquimiaMaxStringLength);
  CreateAlquimiaInterface(name, &interface, &status);
  assert(status.error == kAlquimiaErrorInvalidEngine);
  assert(interface.Setup == NULL);
  assert(interface.Shutdown == NULL);
  assert(interface.ProcessCondition == NULL);
  assert(interface.ReactionStepOperatorSplit == NULL);
  assert(interface.GetAuxiliaryOutput == NULL);
  assert(interface.GetProblemMetaData == NULL);


}  // end test_CreateAlquimiaInterface()

int main(int argc, char** argv) {
  (void) argc;
  (void) argv;

  printf("Testing alquimia c utilites.\n");

  test_AlquimiaCaseInsensitiveStringCompare();
  test_AlquimiaVectors();
  test_AlquimiaNameIndexMapping();
  test_CreateAlquimiaInterface();

  printf("All tests passed.\n");
  return EXIT_SUCCESS;
}  // end main()
