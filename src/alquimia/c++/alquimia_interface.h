/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_CXX_ALQUIMIA_INTERFACE_H_
#define ALQUIMIA_CXX_ALQUIMIA_INTERFACE_H_

/*******************************************************************************
 **
 ** C++ definition of the alquimia interface as an abstract base class
 **
 ******************************************************************************/

#include "alquimia_containers.h"

// Base class defining the alquimia interface
namespace alquimia {

class AlquimiaInterface {
 public:
  AlquimiaInterface() {};
  virtual ~AlquimiaInterface() {};

  // read data files/structures, initialize memory, basis management
  // (includes reading database, swapping basis, etc.)
  virtual void Setup(const std::string& input_file,
                     AlquimiaMetaData* meta_data,
                     AlquimiaEngineStatus* status) = 0;

  // constrain processing for boundary/initial constraints. Called
  // once for each IC/BC.
  virtual void ProcessConstraint(const AlquimiaGeochemicalCondition& condition,
                                 AlquimiaState* state) = 0;

  // take one (or more?) reaction steps in operator split mode
  virtual void ReactionStepOperatorSplit(const double delta_t,
                                         const AlquimiaAuxiliaryData& aux_data,
                                         const AlquimiaMaterialProperties& material_props,
                                         AlquimiaState* state,
                                         AlquimiaEngineStatus* status) = 0;

  // Access to user selected geochemical data for output, i.e. pH,
  // mineral SI, reaction rates
  virtual void GetAuxiliaryOutput(AlquimiaAuxiliaryData* aux_data) = 0;

 protected:

 private:

};

}  // namespace alquimia
#endif  // ALQUIMIA_CXX_ALQUIMIA_INTERFACE_H_
