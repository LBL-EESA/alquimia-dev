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
  void pflotranalquimia_setup_ (char* input_filename,
                                AlquimiaSizes_C* sizes) ;
  void pflotranalquimia_processconstraint_ () ;
  void pflotranalquimia_reactionstepoperatorsplit_ () ;
  void pflotranalquimia_getauxiliaryoutput_ () ;
  void pflotranalquimia_getenginemetadata_ (AlquimiaSizes_C* sizes,
                                            AlquimiaMetaData_C* metadata) ;
  void pflotranalquimia_getprimarynamefromindex_(
      int* primary_index, char* primary_name);
}

namespace alquimia {


class PFloTranAlquimiaInterface : public AlquimiaInterface {
 public:
  PFloTranAlquimiaInterface();
  virtual ~PFloTranAlquimiaInterface();

  void Setup(const std::string& input_file,
             AlquimiaSizes_C* sizes);

  void ProcessCondition(const AlquimiaGeochemicalCondition_C& condition,
                         AlquimiaState_C* state);

  void ReactionStepOperatorSplit(
      const double delta_t,
      const AlquimiaAuxiliaryData_C& aux_data,
      const AlquimiaMaterialProperties_C& material_props,
      AlquimiaState_C* state,
      AlquimiaEngineStatus_C* status);

  void GetAuxiliaryOutput(AlquimiaAuxiliaryData_C* aux_data);


  void GetEngineMetaData(AlquimiaSizes_C* sizes,
                         AlquimiaMetaData_C* meta_data);

 protected:

 private:

};

}  // namespace alquimia
#endif  // ALQUIMIA_CXX_PFLOTRAN_INTERFACE_H_
