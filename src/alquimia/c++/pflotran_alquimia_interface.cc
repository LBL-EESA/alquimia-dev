/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "pflotran_alquimia_interface.h"

#include "pflotran_alquimia_interface.h"
#include "alquimia_containers.h"

namespace alquimia {

PFloTranAlquimiaInterface::PFloTranAlquimiaInterface() 
    : AlquimiaInterface() {

}  // end PFloTranAlquimiaInterface()

PFloTranAlquimiaInterface::~PFloTranAlquimiaInterface() {
  delete pflotran_data_;
}  // end ~PFloTranAlquimiaInterface()

void PFloTranAlquimiaInterface::Setup(
    const std::string& input_file,
    AlquimiaMetaData* meta_data,
    AlquimiaEngineStatus* status) {
  // call pflotran's init function
  
  // get the size of data that pflotran expects to be passed
  int data_size = 1;
  SetupDataTransfer(data_size);
}  // end Setup()

void PFloTranAlquimiaInterface::ProcessConstraint(
    const AlquimiaGeochemicalCondition& condition,
    AlquimiaState* state) {

}  // end ProcessConstraint()

void PFloTranAlquimiaInterface::ReactionStepOperatorSplit(
    const double delta_t,
    const AlquimiaAuxiliaryData& aux_data,
    const AlquimiaMaterialProperties& material_props,
    AlquimiaState* state,
    AlquimiaEngineStatus* status) {

}  // end ReactionStepOperatorSplit()

void PFloTranAlquimiaInterface::GetAuxiliaryOutput(
    AlquimiaAuxiliaryData* aux_data) {

}  // end GetAuxiliaryOutput()

void PFloTranAlquimiaInterface::SetupDataTransfer(const int data_size) {
  data_transfer_size(data_size);
  pflotran_data_ = new double[data_transfer_size()];
}  // end SetupDataTransfer()

void PFloTranAlquimiaInterface::MungeDataAlquimiaToPflotran(
    const AlquimiaAuxiliaryData& aux_data,
    const AlquimiaMaterialProperties& material_props,
    const AlquimiaState& state) {

}  // end MungeDataAlquimiaToPflotran()

void PFloTranAlquimiaInterface::MungeDataPflotranToAlquimia(
    AlquimiaAuxiliaryData* aux_data,
    AlquimiaState* state) {

}  // end MungeDataPflotranToAlquimia()

}  //  namespace alquimia

