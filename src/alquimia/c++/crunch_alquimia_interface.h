/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_CXX_CRUNCH_INTERFACE_H_
#define ALQUIMIA_CXX_CRUNCH_INTERFACE_H_

/*******************************************************************************
 **
 ** C++ implementation of the crunch alquimia interface
 **
 ******************************************************************************/

#include "alquimia_interface.h"

namespace alquimia {

class CrunchAlquimiaInterface : public AlquimiaInterface {
 public:
  CrunchAlquimiaInterface();
  virtual ~CrunchAlquimiaInterface();

  void Setup(const std::string& input_file,
             AlquimiaMetaData* meta_data,
             AlquimiaEngineStatus* status);

  void ProcessConstraint(const AlquimiaGeochemicalCondition& condition,
                         AlquimiaState* state);

  void ReactionStepOperatorSplit(const double delta_t,
                                 const AlquimiaAuxiliaryData& aux_data,
                                 const AlquimiaMaterialProperties& material_props,
                                 AlquimiaState* state,
                                 AlquimiaEngineStatus* status);

  void GetAuxiliaryOutput(AlquimiaAuxiliaryData* aux_data);


 protected:

 private:
};

}  // namespace alquimia
#endif  // ALQUIMIA_CXX_CRUNCH_INTERFACE_H_
