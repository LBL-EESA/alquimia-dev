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
    AlquimiaMetaData* meta_data,
    AlquimiaEngineStatus* status) {

}  // end Setup()

void CrunchAlquimiaInterface::ProcessConstraint(
    const AlquimiaGeochemicalCondition& condition,
    AlquimiaState* state) {

}  // end ProcessConstraint()

void CrunchAlquimiaInterface::ReactionStepOperatorSplit(
    const double delta_t,
    const AlquimiaAuxiliaryData& aux_data,
    const AlquimiaMaterialProperties& material_props,
    AlquimiaState* state,
    AlquimiaEngineStatus* status) {

}  // end ReactionStepOperatorSplit()

void CrunchAlquimiaInterface::GetAuxiliaryOutput(
    AlquimiaAuxiliaryData* aux_data) {

}  // end GetAuxiliaryOutput()


}  //  namespace alquimia

