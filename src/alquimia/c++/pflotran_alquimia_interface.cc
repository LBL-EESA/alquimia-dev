/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "pflotran_alquimia_interface.h"

#include <iostream>
#include <cstring>
#include <string>

#include "alquimia_containers.h"
#include "alquimia_util.h"

#include "pflotran_alquimia_interface.h"

namespace alquimia {

/******************************************************************************/
PFloTranAlquimiaInterface::PFloTranAlquimiaInterface() 
    : AlquimiaInterface() {

}  // end PFloTranAlquimiaInterface()

/******************************************************************************/
PFloTranAlquimiaInterface::~PFloTranAlquimiaInterface() {

}  // end ~PFloTranAlquimiaInterface()

/******************************************************************************/
void PFloTranAlquimiaInterface::Setup(
    const std::string& input_file,
    AlquimiaSizes_C* sizes,
    AlquimiaEngineStatus_C* status) {
  std::cout << "PFloTranAlquimiaInterface::Setup() :  '"
            << input_file << "'\n";
  // copy the c++ string into a c-style char* that can be passed to fortran
  char* inputfile = new char [input_file.size() + 1];
  strcpy(inputfile, input_file.c_str());

  // NOTE(bja): this is required to prevent crashing, but I'm not
  // entirely sure why it's necessary, and is creating a memory leak?
  void* pft_engine_state = new void*;

  // initialize pflotran
  pflotran_alquimia_setup(inputfile, pft_engine_state, sizes, status);
  set_engine_state(pft_engine_state);
  delete inputfile;
  sizes_ = sizes;
}  // end Setup()

/******************************************************************************/
void PFloTranAlquimiaInterface::Shutdown(AlquimiaEngineStatus_C* status) {
  pflotran_alquimia_shutdown(engine_state(), status);
}  // end Shutdown()

/******************************************************************************/
void PFloTranAlquimiaInterface::ProcessCondition(
    AlquimiaGeochemicalCondition_C* condition,
    AlquimiaMaterialProperties_C* material_props,
    AlquimiaState_C* state,
    AlquimiaAuxiliaryData_C* aux_data,
    AlquimiaEngineStatus_C* status) {
  //std::cout << "PFloTranAlquimiaInterface::ProcessCondition() : " << std::endl;
  //std::cout << "  Processing '" << condition->name << "'" << std::endl;
  pflotran_alquimia_processcondition(engine_state(), condition, material_props,
                                     state, aux_data, status);
}  // end ProcessCondition()

/******************************************************************************/
void PFloTranAlquimiaInterface::ReactionStepOperatorSplit(
    double delta_t,
    AlquimiaMaterialProperties_C* material_props,
    AlquimiaState_C* state,
    AlquimiaAuxiliaryData_C* aux_data,
    AlquimiaEngineStatus_C* status) {

  pflotran_alquimia_reactionstepoperatorsplit(engine_state(),
                                              &delta_t,
                                              material_props,
                                              state,
                                              aux_data,
                                              status);
}  // end ReactionStepOperatorSplit()

/******************************************************************************/
void PFloTranAlquimiaInterface::GetAuxiliaryOutput(
    AlquimiaAuxiliaryData_C* aux_data,
    AlquimiaEngineStatus_C* status) {
  static_cast<void>(status);
  static_cast<void>(aux_data);
}  // end GetAuxiliaryOutput()

/******************************************************************************/
void PFloTranAlquimiaInterface::GetEngineMetaData(
    AlquimiaSizes_C* sizes,
    AlquimiaMetaData_C* meta_data,
    AlquimiaEngineStatus_C* status) {
  //std::cout << "PFloTranAlquimiaInterface::GetEngineMetaData() :\n";

  pflotran_alquimia_getenginemetadata(engine_state(), sizes, meta_data, status);

  // NOTE(bja) : can't get arrays of strings to pass gracefully
  // between c and fortran, so for now we loop through and request the
  // names one at a time

  // NOTE(bja): the indices in meta_data.primary_indices already have
  // the engine base, so we don't need to do any conversions!
  for (int i = 0; i < sizes->num_primary; ++i) {
    pflotran_alquimia_getprimarynamefromindex(
        engine_state(),
        &(meta_data->primary_indices[i]), 
        meta_data->primary_names[i],
        status);
  }
}  // end GetEngineMetaData()

}  //  namespace alquimia

