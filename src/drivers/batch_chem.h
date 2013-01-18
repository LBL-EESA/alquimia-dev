/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_DRIVERS_CC_BATCH_CHEM_H_
#define ALQUIMIA_DRIVERS_CC_BATCH_CHEM_H_

#include <string>
#include <vector>

#include "alquimia_containers.h"

#include "demo_containers.h"

int CommandLineOptions(int argc, char** argv,
                       std::string* verbosity_name,
                       std::string* input_file_name,
                       std::string* template_file_name,
                       bool* debug_batch_driver);

void SetupChemistryOutput(void);

void SetTimeUnits(const std::string& output_time_units,
                  char* time_units,
                  double* time_units_conversion);

void WriteOutputHeader(std::fstream* text_output,
                       const char time_units,
                       const AlquimiaSizes_C& sizes,
                       const AlquimiaMetaData_C& meta_data);

void WriteOutput(std::fstream* text_output,
                 const double time,
                 const AlquimiaSizes_C& sizes,
                 const AlquimiaState_C& state);

void CopyDemoStateToAlquimiaState(
    const alquimia::drivers::utilities::DemoState& demo_state,
    const AlquimiaSizes_C& alquimia_sizes,
    AlquimiaState_C* alquimia_state);

void CopyDemoMaterialPropertiesToAlquimiaMaterials(
    const alquimia::drivers::utilities::DemoMaterialProperties& material_props,
    const AlquimiaSizes_C& alquimia_sizes,
    AlquimiaMaterialProperties_C* alquimia_material_props);

void CopyDemoConditionsToAlquimiaConditions(
    const alquimia::drivers::utilities::DemoConditions& demo_conditions,
    AlquimiaGeochemicalConditionList_C* alquimia_conditions);


#endif  /* ALQUIMIA_DRIVERS_CC_BATCH_CHEM_H_ */
