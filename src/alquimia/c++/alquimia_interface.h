/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_CXX_ALQUIMIA_INTERFACE_H_
#define ALQUIMIA_CXX_ALQUIMIA_INTERFACE_H_

/*******************************************************************************
 **
 ** C++ definition of the alquimia interface as an abstract base class
 **
 ******************************************************************************/

#include <string>

#include "alquimia_containers.h"

// Base class defining the alquimia interface
namespace alquimia {

class AlquimiaInterface {
 public:
  AlquimiaInterface() {};
  virtual ~AlquimiaInterface() { };  // delete engine_state_?

  // read data files/structures, initialize memory, basis management
  // (includes reading database, swapping basis, etc.)
  virtual void Setup(const std::string& input_file,
                     AlquimiaSizes_C* sizes,
                     AlquimiaEngineStatus_C* status) = 0;

  virtual void Shutdown(AlquimiaEngineStatus_C* status) = 0;

  // constrain processing for boundary/initial constraints. Called
  // once for each IC/BC.
  virtual void ProcessCondition(AlquimiaGeochemicalCondition_C* condition,
                                AlquimiaMaterialProperties_C* material_props,
                                AlquimiaState_C* state,
                                AlquimiaAuxiliaryData_C* aux_data,
                                AlquimiaEngineStatus_C* status) = 0;

  // take one (or more?) reaction steps in operator split mode
  virtual void ReactionStepOperatorSplit(
      double delta_t,
      AlquimiaMaterialProperties_C* material_props,
      AlquimiaState_C* state,
      AlquimiaAuxiliaryData_C* aux_data,
      AlquimiaEngineStatus_C* status) = 0;

  // Access to user selected geochemical data for output, i.e. pH,
  // mineral SI, reaction rates
  virtual void GetAuxiliaryOutput(AlquimiaAuxiliaryData_C* aux_data,
                                  AlquimiaEngineStatus_C* status) = 0;

  virtual void GetEngineMetaData(AlquimiaSizes_C* sizes,
                                 AlquimiaMetaData_C* meta_data,
                                 AlquimiaEngineStatus_C* status) = 0;

 protected:
  void* engine_state(void) {
    return engine_state_;
  };

  void set_engine_state(void* engine_state) {
    engine_state_ = engine_state;
  }

 private:
  void* engine_state_;

};

}  // namespace alquimia
#endif  // ALQUIMIA_CXX_ALQUIMIA_INTERFACE_H_
