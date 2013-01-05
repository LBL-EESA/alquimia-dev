/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "pflotran_alquimia_interface.h"

#include <iostream>
#include <cstring>
#include <string>

#include "alquimia_containers.h"

#include "pflotran_alquimia_interface.h"

namespace alquimia {

PFloTranAlquimiaInterface::PFloTranAlquimiaInterface() 
    : AlquimiaInterface() {

}  // end PFloTranAlquimiaInterface()

PFloTranAlquimiaInterface::~PFloTranAlquimiaInterface() {

}  // end ~PFloTranAlquimiaInterface()

void PFloTranAlquimiaInterface::Setup(
    const std::string& input_file,
    AlquimiaSizes_C* sizes) {
  std::cout << "PFloTranAlquimiaInterface::Setup() :  '"
            << input_file << "'\n";
  // copy the c++ string into a c-style char* that can be passed to fortran
  char* inputfile = new char [input_file.size() + 1];
  strcpy(inputfile, input_file.c_str());

  // call pflotran's init function
  pflotranalquimia_setup(inputfile, sizes);

  delete inputfile;
}  // end Setup()

void PFloTranAlquimiaInterface::ProcessCondition(
    AlquimiaGeochemicalCondition_C* condition,
    AlquimiaSizes_C* sizes,
    AlquimiaState_C* state) {
  std::cout << "PFloTranAlquimiaInterface::ProcessCondition() : " << std::endl;
  std::cout << "  Processing '" << condition->name << "'" << std::endl;
  pflotranalquimia_processcondition(condition, sizes, state);
}  // end ProcessCondition()

void PFloTranAlquimiaInterface::ReactionStepOperatorSplit(
    const double delta_t,
    const AlquimiaAuxiliaryData_C& aux_data,
    const AlquimiaMaterialProperties_C& material_props,
    AlquimiaState_C* state,
    AlquimiaEngineStatus_C* status) {
  static_cast<void>(delta_t);
  static_cast<void>(aux_data);
  static_cast<void>(material_props);
  static_cast<void>(state);
  static_cast<void>(status);
}  // end ReactionStepOperatorSplit()

void PFloTranAlquimiaInterface::GetAuxiliaryOutput(
    AlquimiaAuxiliaryData_C* aux_data) {
  static_cast<void>(aux_data);
}  // end GetAuxiliaryOutput()

void PFloTranAlquimiaInterface::GetEngineMetaData(
    AlquimiaSizes_C* sizes,
    AlquimiaMetaData_C* meta_data) {
  std::cout << "PFloTranAlquimiaInterface::GetEngineMetaData() :\n";
  pflotranalquimia_getenginemetadata(sizes, meta_data);
  // NOTE(bja): the indices in meta_data.primary_indices already have
  // the engine base, so we don't need to do any conversions!
  for (int i = 0; i < sizes->num_primary; ++i) {
    pflotranalquimia_getprimarynamefromindex(
        &(meta_data->primary_indices[i]), meta_data->primary_names[i]);
  }
}  // end GetEngineMetaData()

}  //  namespace alquimia

