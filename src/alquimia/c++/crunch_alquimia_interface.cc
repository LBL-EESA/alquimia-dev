/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "crunch_alquimia_interface.h"

#include <iostream>

#include "crunch_alquimia_interface.h"
#include "alquimia_containers.h"

namespace alquimia {

CrunchAlquimiaInterface::CrunchAlquimiaInterface() 
    : AlquimiaInterface() {
  std::cout << "CrunchAlquimiaInterface has not been implemented!\n";
}  // end CrunchAlquimiaInterface()

CrunchAlquimiaInterface::~CrunchAlquimiaInterface() {

}  // end ~CrunchAlquimiaInterface()

void CrunchAlquimiaInterface::Setup(
    const std::string& input_file,
    AlquimiaMetaData_C* meta_data,
    AlquimiaSizes_C* sizes) {

}  // end Setup()

void CrunchAlquimiaInterface::ProcessCondition(
    AlquimiaGeochemicalCondition_C* condition,
    AlquimiaSizes_C* sizes,
    AlquimiaState_C* state) {

}  // end ProcessCondition()

void CrunchAlquimiaInterface::ReactionStepOperatorSplit(
    const double delta_t,
    const AlquimiaAuxiliaryData_C& aux_data,
    const AlquimiaMaterialProperties_C& material_props,
    AlquimiaState_C* state,
    AlquimiaEngineStatus_C* status) {

}  // end ReactionStepOperatorSplit()

void CrunchAlquimiaInterface::GetAuxiliaryOutput(
    AlquimiaAuxiliaryData_C* aux_data) {

}  // end GetAuxiliaryOutput()


}  //  namespace alquimia

