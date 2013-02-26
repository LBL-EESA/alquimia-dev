/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_DRIVERS_CC_BATCH_CHEM_H_
#define ALQUIMIA_DRIVERS_CC_BATCH_CHEM_H_

#include <string>
#include <vector>

#include "alquimia_constants.h"
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
                       const AlquimiaProblemMetaData& meta_data,
                       const AlquimiaSizes& sizes,
                       bool* write_pH);

void WriteOutput(std::fstream* text_output,
                 const double time,
                 const AlquimiaState& state,
                 const AlquimiaAuxiliaryOutputData& aux_output,
                 const bool write_pH);

void CopyDemoStateToAlquimiaState(
    const alquimia::drivers::utilities::DemoState& demo_state,
    const AlquimiaProblemMetaData* const alquimia_meta_data,
    AlquimiaState* alquimia_state);

void CopyDemoMaterialPropertiesToAlquimiaMaterials(
    const alquimia::drivers::utilities::DemoMaterialProperties& material_props,
    const AlquimiaProblemMetaData& alquimia_meta_data,
    AlquimiaMaterialProperties* alquimia_material_props);

void CopyDemoConditionsToAlquimiaConditions(
    const alquimia::drivers::utilities::DemoConditions& demo_conditions,
    AlquimiaGeochemicalConditionVector* alquimia_conditions);

void CopyDemoAqueousConstraintsToAlquimia(
    const std::vector<alquimia::drivers::utilities::DemoAqueousConstraint>& demo_aqueous_constraints,
    AlquimiaAqueousConstraintVector* alquimia_aqueous_constraints);

void CopyDemoMineralConstraintsToAlquimia(
    const std::vector<alquimia::drivers::utilities::DemoMineralConstraint>& demo_mineral_constraints,
    AlquimiaMineralConstraintVector* alquimia_mineral_constraints);

#endif  /* ALQUIMIA_DRIVERS_CC_BATCH_CHEM_H_ */
