/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_INTERFACE_H_
#define ALQUIMIA_INTERFACE_H_

/*******************************************************************************
 **
 ** C implementation of the alquimia interface.
 **
 ******************************************************************************/

#include "alquimia_containers.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  void CreateAlquimiaInterface(const char* engine_name,
                               struct AlquimiaInterface* interface,
                               struct AlquimiaEngineStatus_C* status);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  // ALQUIMIA_INTERFACE_H_
