/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "pflotran_alquimia_interface.h"

#include <iostream>
#include <cstring>
#include <string>

#include "alquimia_containers.h"

#include "pflotran_alquimia_interface.h"

namespace alquimia {

/******************************************************************************/
PFloTranAlquimiaInterface::PFloTranAlquimiaInterface() 
    : AlquimiaInterface() {

}  // end PFloTranAlquimiaInterface()

/******************************************************************************/
PFloTranAlquimiaInterface::~PFloTranAlquimiaInterface() {
  pflotran_alquimia_shutdown(engine_state());
}  // end ~PFloTranAlquimiaInterface()

/******************************************************************************/
void PFloTranAlquimiaInterface::Setup(
    const std::string& input_file,
    AlquimiaSizes_C* sizes) {
  std::cout << "PFloTranAlquimiaInterface::Setup() :  '"
            << input_file << "'\n";
  // copy the c++ string into a c-style char* that can be passed to fortran
  char* inputfile = new char [input_file.size() + 1];
  strcpy(inputfile, input_file.c_str());

  // NOTE(bja): this is required to prevent crashing, but I'm not
  // entirely sure why it's necessary, and is creating a memory leak?
  void* pft_engine_state = new void*;

  // initialize pflotran
  pflotran_alquimia_setup(inputfile, pft_engine_state, sizes);
  set_engine_state(pft_engine_state);
  delete inputfile;
}  // end Setup()

/******************************************************************************/
void PFloTranAlquimiaInterface::ProcessCondition(
    AlquimiaGeochemicalCondition_C* condition,
    AlquimiaSizes_C* sizes,
    AlquimiaState_C* state,
    AlquimiaEngineStatus_C* status) {
  std::cout << "PFloTranAlquimiaInterface::ProcessCondition() : " << std::endl;
  std::cout << "  Processing '" << condition->name << "'" << std::endl;
  pflotran_alquimia_processcondition(engine_state(), condition, sizes, state, status);
}  // end ProcessCondition()

/******************************************************************************/
void PFloTranAlquimiaInterface::ReactionStepOperatorSplit(
    const double delta_t,
    AlquimiaMaterialProperties_C* material_props,
    AlquimiaState_C* state,
    AlquimiaAuxiliaryData_C* aux_data,
    AlquimiaEngineStatus_C* status) {
  pflotran_alquimia_reactionstepoperatorsplit(engine_state(),
                                              delta_t,
                                              material_props,
                                              state,
                                              aux_data,
                                              status);
}  // end ReactionStepOperatorSplit()

/******************************************************************************/
void PFloTranAlquimiaInterface::GetAuxiliaryOutput(
    AlquimiaAuxiliaryData_C* aux_data) {
  static_cast<void>(aux_data);
}  // end GetAuxiliaryOutput()

/******************************************************************************/
void PFloTranAlquimiaInterface::GetEngineMetaData(
    AlquimiaSizes_C* sizes,
    AlquimiaMetaData_C* meta_data) {
  std::cout << "PFloTranAlquimiaInterface::GetEngineMetaData() :\n";
  pflotran_alquimia_getenginemetadata(engine_state(), sizes, meta_data);

  // NOTE(bja) : can't get arrays of strings to pass gracefully
  // between c and fortran, so for now we loop through and request the
  // names one at a time

  // NOTE(bja): the indices in meta_data.primary_indices already have
  // the engine base, so we don't need to do any conversions!
  for (int i = 0; i < sizes->num_primary; ++i) {
    pflotran_alquimia_getprimarynamefromindex(engine_state(),
        &(meta_data->primary_indices[i]), meta_data->primary_names[i]);
  }
}  // end GetEngineMetaData()

}  //  namespace alquimia

