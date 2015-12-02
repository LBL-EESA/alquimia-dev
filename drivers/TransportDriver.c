/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

//
// Alquimia Copyright (c) 2013-2015, The Regents of the University of California, 
// through Lawrence Berkeley National Laboratory (subject to receipt of any 
// required approvals from the U.S. Dept. of Energy).  All rights reserved.
// 
// Alquimia is available under a BSD license. See LICENSE.txt for more
// information.
//
// If you have questions about your rights to use or distribute this software, 
// please contact Berkeley Lab's Technology Transfer and Intellectual Property 
// Management at TTD@lbl.gov referring to Alquimia (LBNL Ref. 2013-119).
// 
// NOTICE.  This software was developed under funding from the U.S. Department 
// of Energy.  As such, the U.S. Government has been granted for itself and 
// others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide 
// license in the Software to reproduce, prepare derivative works, and perform 
// publicly and display publicly.  Beginning five (5) years after the date 
// permission to assert copyright is obtained from the U.S. Department of Energy, 
// and subject to any subsequent five (5) year renewals, the U.S. Government is 
// granted for itself and others acting on its behalf a paid-up, nonexclusive, 
// irrevocable, worldwide license in the Software to reproduce, prepare derivative
// works, distribute copies to the public, perform publicly and display publicly, 
// and to permit others to do so.
//

#include "alquimia/alquimia_interface.h"
#include "alquimia/alquimia_memory.h"
#include "TransportDriver.h"

struct TransportDriver 
{
  TransportInput* input;

  // Simulation parameters.
  TransportInputCoupling coupling;
  double t_min, t_max, cfl;
  int max_steps, order;

  // 1D grid information.
  int num_cells;
  double x_min, x_max;

  // Flow velocity.
  double ux;

  // Current sim state.
  double time;
  int step;

  // Solution vector and auxiliary data.
  int num_primary;
  AlquimiaVectorDouble solution;
  AlquimiaVectorDouble aux_doubles;
  AlquimiaVectorInt aux_ints;

  // Chemistry engine -- one of each of these per thread in general.
  AlquimiaInterface chem;
  AlquimiaEngineStatus chem_status;

  // Chemistry data.
  AlquimiaData chem_data;
  AlquimiaGeochemicalCondition chem_ic;
  AlquimiaGeochemicalConditionVector chem_bcs;
};

TransportDriver* TransportDriver_New(TransportInput* input)
{
  TransportDriver* driver = malloc(sizeof(TransportDriver));
  driver->input = input;

  // Get basic simulation parameters.
  TransportInput_GetSimParameters(driver->input, &driver->coupling, 
                                  &driver->t_min, &driver->t_max, 
                                  &driver->max_steps, &driver->cfl, &driver->order);

  // Get grid information.
  TransportInput_GetDomain(input, &driver->x_min, &driver->x_max, &driver->num_cells);

  // Set up the chemistry engine.
  AllocateAlquimiaEngineStatus(&driver->chem_status);
  char *chem_engine, *chem_input, *chem_ic, *chem_left_bc, *chem_right_bc;
  TransportInput_GetChemistry(driver->input, &chem_engine, &chem_input,
                              &chem_ic, &chem_left_bc, &chem_right_bc);
  CreateAlquimiaInterface(chem_engine, &driver->chem, &driver->chem_status);
  if (driver->chem_status.error != 0) 
  {
    alquimia_error("TransportDriver_New: %s", driver->chem_status.message);
    return NULL;
  }

  // Set up the engine and get storage requirements.
  driver->chem.Setup(chem_input,
                     &driver->chem_data.engine_state,
                     &driver->chem_data.sizes,
                     &driver->chem_data.functionality,
                     &driver->chem_status);
  if (driver->chem_status.error != 0) 
  {
    alquimia_error("TransportDriver_New: %s", driver->chem_status.message);
    return NULL;
  }

  // If you want multiple copies of the chemistry engine with
  // OpenMP, verify: chem_data.functionality.thread_safe == true,
  // then create the appropriate number of chem status and data
  // objects.

  // chem_data.sizes was set by Setup(), so now we can allocate
  // memory for alquimia data transfer containers.
  AllocateAlquimiaData(&driver->chem_data);

  // Initialize a solution vector.
  driver->num_primary = driver->chem_data.sizes.num_primary;
  AllocateAlquimiaVectorDouble(driver->num_primary * driver->num_cells, 
                               &driver->solution);

  // Initialize auxiliary data.
  AllocateAlquimiaVectorDouble(driver->chem_data.sizes.num_aux_doubles * driver->num_cells,
                               &driver->aux_doubles);
  AllocateAlquimiaVectorInt(driver->chem_data.sizes.num_aux_integers * driver->num_cells,
                            &driver->aux_ints);

  return driver;
}

void TransportDriver_Free(TransportDriver* driver)
{
  // Destroy input data.
  TransportInput_Free(driver->input);

  // Destroy solution vector and aux data.
  FreeAlquimiaVectorDouble(&driver->solution);
  FreeAlquimiaVectorDouble(&driver->aux_doubles);
  FreeAlquimiaVectorInt(&driver->aux_ints);

  // Destroy chemistry stuff.
  driver->chem.Shutdown(&driver->chem_data.engine_state, 
                        &driver->chem_status);
  FreeAlquimiaData(&driver->chem_data);
  FreeAlquimiaEngineStatus(&driver->chem_status);
  free(driver);
}

static int TransportDriver_Initialize(TransportDriver* driver)
{
  return 0;
}

static int Run_OperatorSplit(TransportDriver* driver)
{
  double dx = (driver->x_max - driver->x_min) / driver->num_cells;
  double dt = driver->cfl * dx / driver->ux;

  // Initialize the chemistry state in each cell, and set up the solution vector.
  int status = TransportDriver_Initialize(driver);
  if (status != 0)
    return status;

  driver->time = driver->t_min;
  driver->step = 0;
  while ((driver->time < driver->t_max) && (driver->step < driver->max_steps))
  {
    // Do the advection step.
    if (driver->order > 1) // TVD scheme 
    {
    }
    else // simple upwinding
    {
    }

    // Do the chemistry step, using the advection step as input.

    if (status != 0) break;

    driver->time += dt;
    driver->step += 1;
  }

  return status;
}

static int Run_GlobalImplicit(TransportDriver* driver)
{
  alquimia_error("Globally implicit transport is not yet supported.");
  return -1;
}

int TransportDriver_Run(TransportDriver* driver)
{
  if (driver->coupling == TRANSPORT_OPERATOR_SPLIT)
    return Run_OperatorSplit(driver);
  else
    return Run_GlobalImplicit(driver);
}

void TransportDriver_GetSoluteAndAuxData(TransportDriver* driver,
                                         double* time,
                                         AlquimiaVectorString* var_names,
                                         AlquimiaVectorDouble* var_data)
{
  // FIXME: Right now we only write out the solution.
  AllocateAlquimiaVectorString(driver->num_primary, var_names);
  AllocateAlquimiaVectorDouble(driver->num_primary * driver->num_cells, 
                               var_data);
}

#if 0
// Demo-specific helper rountines for setting up alquimia structs
void CopyDemoStateToAlquimiaState(
    const alquimia::drivers::utilities::DemoState& demo_state,
    const AlquimiaProblemMetaData* const alquimia_meta_data,
    AlquimiaState* alquimia_state) {
  alquimia_state->water_density = demo_state.water_density;
  alquimia_state->porosity = demo_state.porosity;
  alquimia_state->temperature = demo_state.temperature;
  alquimia_state->aqueous_pressure = demo_state.aqueous_pressure;
  assert(static_cast<int>(demo_state.cec.size()) == 
         alquimia_state->cation_exchange_capacity.size);
  for (size_t i = 0; i < demo_state.cec.size(); ++i) {
    alquimia_state->cation_exchange_capacity.data[i] = demo_state.cec.at(i);
  }
  assert(static_cast<int>(demo_state.site_density.size()) == 
         alquimia_state->surface_site_density.size);

  char* name;
  name = (char*) calloc(kAlquimiaMaxStringLength, sizeof(char));
  // loop through the *engine's* surface site list
  for (int i = 0; i < alquimia_meta_data->surface_site_names.size; ++i) {
    strncpy(name, alquimia_meta_data->surface_site_names.data[i],
            kAlquimiaMaxStringLength);
    std::map<std::string, double>::const_iterator site;
    site = demo_state.site_density.find(name);
    if (site != demo_state.site_density.end()) {
      alquimia_state->surface_site_density.data[i] = site->second;
    } else {
      std::stringstream message;
      message << "ERROR: chemistry engine expects surface site '" 
              << alquimia_meta_data->surface_site_names.data[i]
              << "', but it was not found in the driver state surface site list."
              << std::endl;
      throw std::runtime_error(message.str());
    }
  }
  free(name);
}  // end CopyDemoStateToAlquimiaState()

void CopyDemoPropertiesToAlquimiaMaterials(
    const alquimia::drivers::utilities::DemoProperties& demo_props,
    const AlquimiaProblemMetaData& alquimia_meta_data,
    AlquimiaProperties* alquimia_props) {
  alquimia_props->volume = demo_props.volume;
  alquimia_props->saturation = demo_props.saturation;

  if (static_cast<size_t>(alquimia_meta_data.isotherm_species_names.size) != 
      demo_props.isotherm_species.size()) {
    std::stringstream message;
    message << "ERROR: chemistry engine expects " 
            << alquimia_meta_data.isotherm_species_names.size << " isotherm species. "
            << " but input file only contains " 
            << demo_props.isotherm_species.size() << std::endl;
    throw std::runtime_error(message.str());
  }

  char* name;
  name = (char*) calloc(kAlquimiaMaxStringLength, sizeof(char));
  // loop through each species in the *engine's* isotherm list
  for (int i = 0; i < alquimia_meta_data.isotherm_species_names.size; ++i) {
    // save the isotherm species id
    strncpy(name, alquimia_meta_data.isotherm_species_names.data[i],
            kAlquimiaMaxStringLength);
    for (size_t j = 0; j < demo_props.isotherm_species.size(); ++j) {
      if (demo_props.isotherm_species.at(j).compare(name) == 0) {
        
        alquimia_props->isotherm_kd.data[i] = 
            demo_props.isotherm_kd.at(j);
        alquimia_props->freundlich_n.data[i] = 
            demo_props.freundlich_n.at(j);
        alquimia_props->langmuir_b.data[i] = 
            demo_props.langmuir_b.at(j);
        break;
      }
    }
  }
  free(name);
}  // end CopyDemoPropertiesToAlquimiaMaterials()


void CopyDemoConditionsToAlquimiaConditions(
    const alquimia::drivers::utilities::DemoConditions& demo_conditions,
    AlquimiaGeochemicalConditionVector* alquimia_conditions) {
  namespace util = alquimia::drivers::utilities;

  AllocateAlquimiaGeochemicalConditionVector(demo_conditions.size(),
                                             alquimia_conditions);

  // copy the geochemical conditions
  util::DemoConditions::const_iterator demo_cond;
  unsigned int i_cond;
  for (demo_cond = demo_conditions.begin(), i_cond = 0;
       i_cond < demo_conditions.size(); ++i_cond, ++demo_cond) {
    std::cout << "    " << demo_cond->first << " : " << i_cond << std::endl;
    AlquimiaGeochemicalCondition* condition =
        &(alquimia_conditions->data[i_cond]);
    AllocateAlquimiaGeochemicalCondition(kAlquimiaMaxStringLength,
                                         demo_cond->second.aqueous_constraints.size(),
                                         demo_cond->second.mineral_constraints.size(),
                                         condition);
    unsigned int max_copy_length = std::min(kAlquimiaMaxStringLength, static_cast<int>(demo_cond->first.size()));
    strncpy(condition->name, demo_cond->first.c_str(), max_copy_length);
    CopyDemoAqueousConstraintsToAlquimia(demo_cond->second.aqueous_constraints,
                                         &condition->aqueous_constraints);
    CopyDemoMineralConstraintsToAlquimia(demo_cond->second.mineral_constraints,
                                         &condition->mineral_constraints);
  }
}  // end CopyDemoConditionsToAlquimiaConditions()

void CopyDemoAqueousConstraintsToAlquimia(
    const std::vector<alquimia::drivers::utilities::DemoAqueousConstraint>& demo_aqueous_constraints,
    AlquimiaAqueousConstraintVector* alquimia_aqueous_constraints) {
  // loop through aqueous constraints
  for (unsigned int i = 0; i < demo_aqueous_constraints.size(); ++i) {

    // easier to work with a pointer to one constraint than the entire vector
    AlquimiaAqueousConstraint* constraint =
        &(alquimia_aqueous_constraints->data[i]);
    
    // allocate memory fo the current constraint
    AllocateAlquimiaAqueousConstraint(constraint);

    //
    // copy demo constraint to alquimia constraint
    //

    // name
    int max_copy_length = std::min(kAlquimiaMaxStringLength,
                                   static_cast<int>(demo_aqueous_constraints.at(i).primary_species_name.size()));
    std::strncpy(constraint->primary_species_name,
                 demo_aqueous_constraints.at(i).primary_species_name.c_str(),
                 max_copy_length);

    // constraint type
    max_copy_length = std::min(kAlquimiaMaxStringLength,
                               static_cast<int>(demo_aqueous_constraints.at(i).constraint_type.size()));
    std::strncpy(constraint->constraint_type,
                 demo_aqueous_constraints.at(i).constraint_type.c_str(),
                 max_copy_length);

    // associated species
    max_copy_length = std::min(kAlquimiaMaxStringLength,
                               static_cast<int>(demo_aqueous_constraints.at(i).associated_species.size()));
    std::strncpy(constraint->associated_species,
                 demo_aqueous_constraints.at(i).associated_species.c_str(),
                 max_copy_length);

    // constraint value
    constraint->value = demo_aqueous_constraints.at(i).value;
  }
}  // end CopyDemoAqueousConstraintsToAlquimia()

void CopyDemoMineralConstraintsToAlquimia(
    const std::vector<alquimia::drivers::utilities::DemoMineralConstraint>& demo_mineral_constraints,
    AlquimiaMineralConstraintVector* alquimia_mineral_constraints) {

  //assert(demo_mineral_constraints.size() == alquimia_mineral_constraints->size);

  // loop through mineral constraints
  for (unsigned int i = 0; i < demo_mineral_constraints.size(); ++i) {

    // easier to work with a single pointer instead of the full array
    AlquimiaMineralConstraint* constraint =
        &(alquimia_mineral_constraints->data[i]);

    // create memory in the constraint
    AllocateAlquimiaMineralConstraint(constraint);

    //
    // copy demo constraint to alquimia constraint
    //

    // name
    int max_copy_length = std::min(kAlquimiaMaxStringLength,
                                   static_cast<int>(demo_mineral_constraints.at(i).mineral_name.size()));
    std::strncpy(constraint->mineral_name,
                 demo_mineral_constraints.at(i).mineral_name.c_str(),
                 max_copy_length);
    
    constraint->volume_fraction = demo_mineral_constraints.at(i).volume_fraction;
    constraint->specific_surface_area = demo_mineral_constraints.at(i).specific_surface_area;
  }
}  // end CopyDemoMineralConstraintsToAlquimia()

int RunSimulation(TransportSim* sim, FILE* output)
{
  // Create the alquimia structures.  NOTE: chem_data are for a single
  // grid cell. For openmp/threaded code, you'll need one chem_status
  // and chem_data per thead, but we don't know if the engine is
  // thread safe until after the call to Setup()!
  //
  AlquimiaInterface chem;
  AlquimiaEngineStatus chem_status;
  AlquimiaData chem_data;
  AlquimiaGeochemicalConditionVector alquimia_conditions;

  // All alquimia function calls require a status object.
  AllocateAlquimiaEngineStatus(&chem_status);
  // Create the chemistry engine
  CreateAlquimiaInterface(demo_simulation.engine.c_str(), &chem, &chem_status);
  if (chem_status.error != 0) {
    std::cout << chem_status.message << std::endl;
    return chem_status.error;
  }

  // setup the engine and get the memory requirements
  chem.Setup(demo_simulation.engine_inputfile.c_str(),
             &chem_data.engine_state,
             &chem_data.sizes,
             &chem_data.functionality,
             &chem_status);
  if (chem_status.error != 0) {
    std::cout << chem_status.message << std::endl;
    PrintAlquimiaSizes(&chem_data.sizes);
    return chem_status.error;
  }

  // if you want multiple copies of the chemistry engine with
  // OpenMP, verify: chem_data.functionality.thread_safe == true,
  // then create the appropriate number of chem status and data
  // objects

  // chem_data.sizes was set by Setup(), so now we can allocate
  // memory for alquimia data transfer containers.
  AllocateAlquimiaData(&chem_data);

  // allocate the remaining memory in the driver (mesh dependent)

  // get the problem meta data (species and mineral names, etc)
  chem.GetProblemMetaData(&chem_data.engine_state,
                          &chem_data.meta_data,
                          &chem_status);
  if (chem_status.error != 0) {
    std::cout << chem_status.message << std::endl;
    PrintAlquimiaProblemMetaData(&chem_data.meta_data);
    return chem_status.error;
  }
  //PrintAlquimiaProblemMetaData(&chem_data.meta_data);

  // finish initializing the driver, e.g. verify material
  // properties, species names, etc

  // initialize the alquimia state and material properties with
  // appropriate values from the driver's memory.
  CopyDemoStateToAlquimiaState(demo_state, &chem_data.meta_data,
                               &chem_data.state);
  CopyDemoPropertiesToAlquimiaMaterials(
      demo_props, chem_data.meta_data, &chem_data.properties);

  PrintAlquimiaData(&chem_data);

  //
  // prepare for constraint processing
  //

  // Read the geochemical conditions from the driver's native format
  // and store them in alquimia's format
  CopyDemoConditionsToAlquimiaConditions(demo_conditions, &alquimia_conditions);

  PrintAlquimiaGeochemicalConditionVector(&alquimia_conditions);

  for (int i = 0; i < alquimia_conditions.size; ++i) {
    // ask the engine to process the geochemical conditions
    if (demo_simulation.initial_condition.compare(
            alquimia_conditions.data[i].name) == 0) {
      // for batch, only care about the IC. If the conditions have
      // spatially dependent values, then the state and material
      // properties need to be updated here!
      chem.ProcessCondition(&chem_data.engine_state,
                            &(alquimia_conditions.data[i]),
                            &chem_data.properties,
                            &chem_data.state,
                            &chem_data.aux_data,
                            &chem_status);
      if (chem_status.error != 0) {
        PrintAlquimiaData(&chem_data);
        PrintAlquimiaGeochemicalCondition(&(alquimia_conditions.data[i]));
        std::cout << chem_status.message << std::endl;
        return chem_status.error;
      }
    }
    // store the processed geochemical conditions in driver's memory
  }
  // we are done with the conditions (data is stored in driver state
  // vectors at this point).
  FreeAlquimiaGeochemicalConditionVector(&alquimia_conditions);

  char time_units;
  double time_units_conversion;  // [time_units / sec]
  SetTimeUnits(demo_simulation.time_units,
               &time_units, &time_units_conversion);

  // set delta t: [time_units]/[time_units/sec] = [sec]
  double delta_t = demo_simulation.delta_t / time_units_conversion;
  double time = 0.0;
  // get initial pH
  chem.GetAuxiliaryOutput(&chem_data.engine_state,
                          &chem_data.properties,
                          &chem_data.state,
                          &chem_data.aux_data,
                          &chem_data.aux_output,
                          &chem_status);
  // save the IC to our output file
  output->WriteHeader(time_units, chem_data.meta_data, chem_data.sizes);
  output->Write(time, chem_data.state, chem_data.aux_output);
  std::cout << "Starting reaction stepping (OS) with dt = " << delta_t << " [s]\n";
  for (int t = 0; t < demo_simulation.num_time_steps; ++t) {
    time += delta_t;
    //std::cout << "reaction step : " << t << "  time: " << time << std::endl;
    if (false) {
      PrintAlquimiaState(&chem_data.state);
      PrintAlquimiaAuxiliaryData(&chem_data.aux_data);
    }
    // unpack from driver memory, since this is batch, no unpacking
    chem.ReactionStepOperatorSplit(&chem_data.engine_state,
                                   &delta_t,
                                   &chem_data.properties,
                                   &chem_data.state,
                                   &chem_data.aux_data,
                                   &chem_status);
    if (chem_status.error != 0) {
      std::cout << chem_status.message << std::endl;
      PrintAlquimiaState(&chem_data.state);
      PrintAlquimiaAuxiliaryData(&chem_data.aux_data);
      return chem_status.error;
    }
    chem.GetAuxiliaryOutput(&chem_data.engine_state,
                            &chem_data.properties,
                            &chem_data.state,
                            &chem_data.aux_data,
                            &chem_data.aux_output,
                            &chem_status);
    if (chem_status.error != 0) {
      std::cout << chem_status.message << std::endl;
      return chem_status.error;
    }
    double out_time = time * time_units_conversion;  // [sec]*[time_units/sec]
    output->Write(out_time, chem_data.state, chem_data.aux_output);
    std::cout << "  step = " << std::setw(6) << std::right << t 
              << "    time = " << std::setw(8) << std::right
              << time*time_units_conversion << " [" << time_units 
              << "]    newton = " << std::setw(4) << std::right 
              << chem_status.num_newton_iterations << std::endl; 
    // repack into driver memory...
  }
  std::cout << std::endl;

  // cleanup memory
  chem.Shutdown(&chem_data.engine_state, &chem_status);
  FreeAlquimiaData(&chem_data);
  FreeAlquimiaEngineStatus(&chem_status);
  return 0;
}
#endif

