/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_DRIVERS_CC_BATCH_CHEM_H_
#define ALQUIMIA_DRIVERS_CC_BATCH_CHEM_H_

#include <string>
#include <vector>
#include <fstream>

#include "petscsys.h"

#include "alquimia_constants.h"
#include "alquimia_containers.h"

#include "demo_containers.h"
#include "demo_output.h"

int CommandLineOptions(char const * const argv,
                       std::string* input_file_name,
                       std::string* template_file_name,
                       PetscBool* debug_batch_driver);

void SetupChemistryOutput(void);

void SetTimeUnits(const std::string& output_time_units,
                  char* time_units,
                  double* time_units_conversion);

int BatchChemWithAlquimia(
    const alquimia::drivers::utilities::DemoSimulation& demo_simulation,
    const alquimia::drivers::utilities::DemoState& demo_state,
    const alquimia::drivers::utilities::DemoMaterialProperties& demo_material_props,
    const alquimia::drivers::utilities::DemoConditions& demo_conditions,
    alquimia::drivers::utilities::DemoOutput* output);

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
