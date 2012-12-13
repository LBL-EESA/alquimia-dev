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
    AlquimiaMetaData_C* meta_data,
    AlquimiaSizes_C* sizes) {
  // call pflotran 's init function
  std::cout << "PFloTranAlquimiaInterface::Setup() :  '"
            << input_file << "'\n";
  char* inputfile = new char [input_file.size() + 1];
  strcpy(inputfile, input_file.c_str());
  pflotranalquimia_setup_(inputfile, meta_data, sizes);
  delete inputfile;
  std::cout << "PFloTranAlquimiaInterface::Setup() :\n";
  std::cout << "  sizes :\n"
            << "    num primary species : " << sizes->num_primary << "\n"
            << "    num kinetic minerals : " << sizes->num_kinetic_minerals << "\n"
            << "    num aqueous complexes : " << sizes->num_aqueous_complexes << "\n"
            << "    num surface sites : " << sizes->num_surface_sites << "\n"
            << "    num ion exchange sites : " << sizes->num_ion_exchange_sites << "\n";

  std::cout << "  meta data :\n"
            << "    thread_safe : " << meta_data->thread_safe << "\n"
            << "    temperature_dependent : " << meta_data->temperature_dependent << "\n"
            << "    pressure_dependent : " << meta_data->pressure_dependent << "\n"
            << "    porosity_update  : " << meta_data->porosity_update << "\n";

}  // end Setup()

void PFloTranAlquimiaInterface::ProcessCondition(
    const AlquimiaGeochemicalCondition_C& condition,
    AlquimiaState_C* state) {

}  // end ProcessCondition()

void PFloTranAlquimiaInterface::ReactionStepOperatorSplit(
    const double delta_t,
    const AlquimiaAuxiliaryData_C& aux_data,
    const AlquimiaMaterialProperties_C& material_props,
    AlquimiaState_C* state,
    AlquimiaEngineStatus_C* status) {

}  // end ReactionStepOperatorSplit()

void PFloTranAlquimiaInterface::GetAuxiliaryOutput(
    AlquimiaAuxiliaryData_C* aux_data) {

}  // end GetAuxiliaryOutput()

}  //  namespace alquimia

