/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_CXX_CRUNCH_INTERFACE_H_
#define ALQUIMIA_CXX_CRUNCH_INTERFACE_H_

/*******************************************************************************
 **
 ** C++ implementation of the crunch alquimia interface
 **
 ******************************************************************************/

#include "alquimia_interface.h"

extern "C" {
  void crunchalquimia_setup_ (char* input_filename,
                              AlquimiaMetaData_C* meta_data,
                              AlquimiaSizes_C* sizes) ;
  void crunchalquimia_processconstraint_ () ;
  void crunchalquimia_reactionstepoperatorsplit_ () ;
  void crunchalquimia_getauxiliaryoutput_ () ;
  void crunchalquimia_getenginefunctionality_ (AlquimiaMetaData_C* metadata) ;
}

namespace alquimia {

class CrunchAlquimiaInterface : public AlquimiaInterface {
 public:
  CrunchAlquimiaInterface();
  virtual ~CrunchAlquimiaInterface();

  void Setup(const std::string& input_file,
             AlquimiaSizes_C* sizes);

  void ProcessCondition(AlquimiaGeochemicalCondition_C* condition,
                        AlquimiaMaterialProperties_C* material_props,
                        AlquimiaState_C* state,
                        AlquimiaAuxiliaryData_C* aux_data,
                        AlquimiaEngineStatus_C* status);

  void ReactionStepOperatorSplit(
      double delta_t,
      AlquimiaMaterialProperties_C* material_props,
      AlquimiaState_C* state,
      AlquimiaAuxiliaryData_C* aux_data,
      AlquimiaEngineStatus_C* status);

  void GetAuxiliaryOutput(AlquimiaAuxiliaryData_C* aux_data);

  void GetEngineMetaData(AlquimiaSizes_C* sizes,
                         AlquimiaMetaData_C* meta_data);

 protected:

 private:
};

}  // namespace alquimia
#endif  // ALQUIMIA_CXX_CRUNCH_INTERFACE_H_
