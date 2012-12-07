/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_CXX_PFLOTRAN_INTERFACE_H_
#define ALQUIMIA_CXX_PFLOTRAN_INTERFACE_H_

/*******************************************************************************
 **
 ** C++ implementation of the pflotran alquimia interface
 **
 ******************************************************************************/

#include "alquimia_interface.h"

namespace alquimia {

class PFloTranAlquimiaInterface : public AlquimiaInterface {
 public:
  PFloTranAlquimiaInterface();
  virtual ~PFloTranAlquimiaInterface();

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
  void SetupDataTransfer(const int data_size);
  // NOTE(bja): probably need multiple munge data functions in each
  // direction, i.e. MungeOperatorSplittingDataAlquimiaToPflotran,
  // MungeGlobalImplicitDataAlquimiaToPflotran,
  // MungeConstraintDataAlquimiaToPflotran, etc
  void MungeDataAlquimiaToPflotran(
      const AlquimiaAuxiliaryData& aux_data,
      const AlquimiaMaterialProperties& material_props,
      const AlquimiaState& aux_data);
  void MungeDataPflotranToAlquimia(
      AlquimiaAuxiliaryData* aux_data,
      AlquimiaState* aux_data);

  // NOTE(bja): will probably require multiple data transfer vectors...
  int data_transfer_size(void) {
    return data_transfer_size_;
  }
  void data_transfer_size(const int data_size) {
    data_transfer_size_ = data_size;
  }

  int data_transfer_size_;
  double* pflotran_data_;
};

}  // namespace alquimia
#endif  // ALQUIMIA_CXX_PFLOTRAN_INTERFACE_H_
