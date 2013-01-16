/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_CXX_PFLOTRAN_INTERFACE_H_
#define ALQUIMIA_CXX_PFLOTRAN_INTERFACE_H_

/*******************************************************************************
 **
 ** C++ implementation of the pflotran alquimia interface
 **
 ******************************************************************************/

#include "alquimia_interface.h"

#include "alquimia_containers.h"

extern "C" {
  void pflotran_alquimia_setup(char* input_filename,
                               void* pft_engine_state,
                               AlquimiaSizes_C* sizes) ;
  void pflotran_alquimia_shutdown(void* pft_engine_state);
  void pflotran_alquimia_processcondition(void* pft_engine_state,
                                          AlquimiaGeochemicalCondition_C* condition,
                                          AlquimiaSizes_C* sizes,
                                          AlquimiaState_C* state,
                                          AlquimiaEngineStatus_C* status);
  void pflotran_alquimia_reactionstepoperatorsplit(
      void* pft_engine_state,
      const double delta_t,
      AlquimiaMaterialProperties_C* material_properties,
      AlquimiaState_C* state,
      AlquimiaAuxiliaryData_C* aux_data,
      AlquimiaEngineStatus_C* status);
  void pflotran_alquimia_getauxiliaryoutput(void* pft_engine_state) ;
  void pflotran_alquimia_getenginemetadata(void* pft_engine_state,
                                           AlquimiaSizes_C* sizes,
                                           AlquimiaMetaData_C* metadata) ;
  void pflotran_alquimia_getprimarynamefromindex(void* pft_engine_state,
                                                 int* primary_index,
                                                 char* primary_name);
}

namespace alquimia {


class PFloTranAlquimiaInterface : public AlquimiaInterface {
 public:
  PFloTranAlquimiaInterface();
  virtual ~PFloTranAlquimiaInterface();

  void Setup(const std::string& input_file,
             AlquimiaSizes_C* sizes);

  void ProcessCondition(AlquimiaGeochemicalCondition_C* condition,
                        AlquimiaSizes_C* sizes,
                        AlquimiaState_C* state,
                        AlquimiaEngineStatus_C* status);

  void ReactionStepOperatorSplit(
      const double delta_t,
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
#endif  // ALQUIMIA_CXX_PFLOTRAN_INTERFACE_H_
