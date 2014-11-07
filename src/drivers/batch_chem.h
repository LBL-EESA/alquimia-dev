/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
** Alquimia Copyright (c) 2013, The Regents of the University of California, 
** through Lawrence Berkeley National Laboratory (subject to receipt of any 
** required approvals from the U.S. Dept. of Energy).  All rights reserved.
** 
** Alquimia is available under a BSD license. See LICENSE.txt for more
** information.
**
** If you have questions about your rights to use or distribute this software, 
** please contact Berkeley Lab's Technology Transfer and Intellectual Property 
** Management at TTD@lbl.gov referring to Alquimia (LBNL Ref. 2013-119).
** 
** NOTICE.  This software was developed under funding from the U.S. Department 
** of Energy.  As such, the U.S. Government has been granted for itself and 
** others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide 
** license in the Software to reproduce, prepare derivative works, and perform 
** publicly and display publicly.  Beginning five (5) years after the date 
** permission to assert copyright is obtained from the U.S. Department of Energy, 
** and subject to any subsequent five (5) year renewals, the U.S. Government is 
** granted for itself and others acting on its behalf a paid-up, nonexclusive, 
** irrevocable, worldwide license in the Software to reproduce, prepare derivative
** works, distribute copies to the public, perform publicly and display publicly, 
** and to permit others to do so.
** 
** Authors: Benjamin Andre <bandre@lbl.gov>
*/

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
