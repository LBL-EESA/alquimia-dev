/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "batch_chem.h"

#ifdef WINDOWS
#include "xgetopt.h"
#else
#include <unistd.h>
#endif

#include <cstdlib>
#include <cctype>
#include <cstring>
#include <cassert>

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>

#include <mpi.h>

#include <petscsys.h>

#include "alquimia_memory.h"
#include "alquimia_util.h"
#include "alquimia_containers.h"
#include "alquimia_interface.h"
#include "alquimia_constants.h"

#include "cfg_reader.h"
#include "demo_containers.h"
#include "demo_utils.h"
#include "string_tokenizer.h"

int main(int argc, char** argv) {
  namespace util = alquimia::drivers::utilities;
  std::stringstream message;

  bool debug_batch_driver(false);
  std::string verbosity_name("");
  std::string input_file_name("");
  std::string template_file_name("");
  int error = EXIT_SUCCESS;

  error = CommandLineOptions(argc, argv,
                             &verbosity_name,
                             &input_file_name,
                             &template_file_name,
                             &debug_batch_driver);
  if (error != EXIT_SUCCESS) {
    exit(error);
  }

  // initialize petsc/mpi for engines that require it (pflotran)
  char help[] = "petsc help string";
  PetscInitialize(&argc, &argv, (char*)0, help);

  try {
    //
    // Read the demo driver input file
    //
    util::DemoConfigReader cfg_reader;

    if (!template_file_name.empty()) {
      cfg_reader.WriteTemplateFile(template_file_name);
      exit(EXIT_SUCCESS);
    }

    util::DemoSimulation demo_simulation;
    util::DemoState demo_state;
    util::DemoMaterialProperties demo_material_props;
    util::DemoConditions demo_conditions;

    if (!input_file_name.empty()) {
      cfg_reader.set_debug(false);
      cfg_reader.ReadInputFile(input_file_name,
                               &demo_simulation, &demo_state,
                               &demo_material_props, &demo_conditions);
      if (debug_batch_driver) {
        std::cout << "- Input File ---------------------------------------------------------\n";
        demo_simulation.Print();
        demo_state.Print();
        demo_material_props.Print();
        util::PrintGeochemicalConditions(demo_conditions);
        std::cout << "--------------------------------------------------------- Input File -\n";
      }
    }

    //
    // open the output file
    //
    std::fstream text_output;
    size_t position = input_file_name.find_last_of('.');
    std::string text_output_name = input_file_name.substr(0, position) + ".txt";
    text_output.open(text_output_name.c_str(), std::fstream::out);

    //
    // run the demo batch chemistry
    //
    error = BatchChemWithAlquimia(demo_simulation, demo_state,
                                  demo_material_props, demo_conditions,
                                  &text_output);

    text_output.close();

  } catch (const std::runtime_error& rt_error) {
    std::cout << rt_error.what();
    error = EXIT_FAILURE;
  } catch (const std::logic_error& lg_error) {
    std::cout << lg_error.what();
    error = EXIT_FAILURE;
  } catch (const std::exception& e) {
    std::cout << e.what();
    error = EXIT_FAILURE;
  }

  PetscInt petsc_error = PetscFinalize();

  if (error == 0 && petsc_error == 0) {
    std::cout << "Success!\n";
  } else {
    std::cout << "Failed!\n";
  }

  return error;
}  // end main()


/*******************************************************************************
 **
 **  Demonstrate setting up the alquimia interface and data objects,
 **  then running simple batch chemistry reaction steps.
 **
 *******************************************************************************/
int BatchChemWithAlquimia(
    const alquimia::drivers::utilities::DemoSimulation& demo_simulation,
    const alquimia::drivers::utilities::DemoState& demo_state,
    const alquimia::drivers::utilities::DemoMaterialProperties& demo_material_props,
    const alquimia::drivers::utilities::DemoConditions& demo_conditions,
    std::fstream* text_output) {
  //
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

  //
  // prepare for constraint processing
  //

  // initialize the alquimia state and material properties with
  // appropriate values from the driver's memory.
  CopyDemoStateToAlquimiaState(demo_state, &chem_data.meta_data,
                               &chem_data.state);
  CopyDemoMaterialPropertiesToAlquimiaMaterials(
      demo_material_props, chem_data.meta_data, &chem_data.material_properties);

  PrintAlquimiaData(&chem_data);

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
                            &chem_data.material_properties,
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
                          &chem_data.material_properties,
                          &chem_data.state,
                          &chem_data.aux_data,
                          &chem_data.aux_output,
                          &chem_status);
  // save the IC to our output file
  bool write_pH;
  WriteOutputHeader(text_output, time_units,
                    chem_data.meta_data, chem_data.sizes, &write_pH);
  WriteOutput(text_output, time,
              chem_data.state, chem_data.aux_output, write_pH);
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
                                   &chem_data.material_properties,
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
                            &chem_data.material_properties,
                            &chem_data.state,
                            &chem_data.aux_data,
                            &chem_data.aux_output,
                            &chem_status);
    if (chem_status.error != 0) {
      std::cout << chem_status.message << std::endl;
      return chem_status.error;
    }
    double out_time = time * time_units_conversion;  // [sec]*[time_units/sec]
    WriteOutput(text_output, out_time,
                chem_data.state, chem_data.aux_output, write_pH);
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

} // end BatchChemWithAlquimia()

/*******************************************************************************
 **
 **  Commandline
 **
 *******************************************************************************/
int CommandLineOptions(int argc, char** argv,
                       std::string* verbosity_name,
                       std::string* input_file_name,
                       std::string* template_file_name,
                       bool* debug_batch_driver)
{
  int error = -2;
  int option;
  extern char* optarg;
  std::stringstream help_message;
  help_message << argv[0]
               << " alquimia batch chemistry demonstration driver.\n";
  help_message << "  command line options:" << std::endl;
  help_message << "    -d" << std::endl;
  help_message << "         debugging flag for batch driver" << std::endl;
  help_message << "    -i string " << std::endl;
  help_message << "         input file name" << std::endl;
  help_message << std::endl;
  help_message << "    -t string" << std::endl;
  help_message << "         write a template input file" << std::endl;
  help_message << std::endl;

  while ((option = getopt(argc, argv, "di:ht:?")) != -1) {
    switch (option) {
      case 'd': {
        *debug_batch_driver = true;
        break;
      }
      case 'i': {
        /* input file name */
        input_file_name->assign(optarg);
        error = EXIT_SUCCESS;
        break;
      }
      case 't': {
        /* template file name */
        template_file_name->assign(optarg);
        error = EXIT_SUCCESS;
        break;
      }
      case 'v': {
        verbosity_name->assign(optarg);
        break;
      }
      case '?':
      case 'h': {  /* help mode */
        /* print some help stuff and exit without doing anything */
        std::cout << help_message.str();
        error = -1;
        break;
      }
      default: {
        /* no options */
        break;
      }
    }
  }

  if (!input_file_name->c_str() && !template_file_name->c_str()) {
    std::cout << "An input or template file name must be specified." << std::endl;
    std::cout << "Run \"" <<  argv[0] << " -h \" for help." << std::endl;
  }

  if (*debug_batch_driver) {
    std::stringstream message;
    message << "- Command Line Options -----------------------------------------------" << std::endl;
    message << "\tdebug batch driver: " << *debug_batch_driver << std::endl;
    message << "\tinput file name: " << *input_file_name << std::endl;
    message << "\ttemplate file name: " << *template_file_name << std::endl;
    message << "----------------------------------------------- Command Line Options -" << std::endl;
    message << std::endl << std::endl;
    std::cout << message.str();
  }
  if (error == -2) {
    std::cout << help_message.str();
  }
  return error;
}  // end commandLineOptions()


/*******************************************************************************
 **
 **  Output related functions
 **
 *******************************************************************************/
void SetTimeUnits(const std::string& output_time_units,
                  char* time_units,
                  double* time_units_conversion) {
  // default values
  *time_units = 's';
  *time_units_conversion = 1.0;

  // do we want to change the time units for the output?
  if (output_time_units.size() > 0) {
    *time_units = std::tolower(output_time_units.at(0));
    switch (*time_units) {
      case 's':
        break;
      case 'm':
        // 1 min --> 60 sec
        *time_units_conversion = 60.0;
        break;
      case 'h':
        // 1 hour --> 60 min --> 3600.0 sec
        *time_units_conversion = 60.0 * 60.0;
        break;
      case 'd':
        *time_units_conversion = 60.0 * 60.0 * 24.0;
        break;
      case 'y':
        *time_units_conversion = 60.0 * 60.0 * 24.0 * 365.25;
        break;
      default:
        break;
    }
  }
  // base units of driver are seconds, so we actually want 1.0/conversion
  *time_units_conversion = 1.0 / (*time_units_conversion);
}  // end SetTimeUnits()

void WriteOutputHeader(std::fstream* text_output, const char time_units,
                       const AlquimiaProblemMetaData& meta_data,
                       const AlquimiaSizes& sizes,
                       bool* write_pH) {
  if (text_output->is_open()) {
    *text_output << "# \"Time [" << time_units << "]\"";
    int h_index;
    AlquimiaFindIndexFromName("H+", &meta_data.primary_names, &h_index);
    if (h_index >= 0) {
      *write_pH = true;
      *text_output << " , \"pH\"";
    }
    for (int i = 0; i < meta_data.primary_names.size; ++i) {
      *text_output <<  " , \"Total " << meta_data.primary_names.data[i] << " [M]\"";
    }

    for (int i = 0; i < sizes.num_sorbed; ++i) {
      *text_output << " , \"Total Sorbed " << meta_data.primary_names.data[i]
                   << " [mol/m^3]\"";
    }

    for (int i = 0; i < meta_data.mineral_names.size; ++i) {
      *text_output << " , \"" << meta_data.mineral_names.data[i] << " VF\"";
    }
    for (int i = 0; i < meta_data.mineral_names.size; ++i) {
      *text_output << " , \"" << meta_data.mineral_names.data[i] << " Rate [mol/m^3/sec]\"";
    }
    *text_output << std::endl;
  }
}  // end WriteOutputHeader()

void WriteOutput(std::fstream* text_output, const double time,
                 const AlquimiaState& state,
                 const AlquimiaAuxiliaryOutputData& aux_output,
                 const bool write_pH) {
  if (text_output->is_open()) {
    std::string seperator("  ");
    *text_output << std::scientific << std::uppercase
                 << std::setprecision(6);
    *text_output << seperator << time;
    if (write_pH) {
      *text_output << seperator << aux_output.pH;
    }
    for (int i = 0; i < state.total_mobile.size; ++i) {
      *text_output << seperator << state.total_mobile.data[i];
    }
    for (int i = 0; i < state.total_immobile.size; ++i) {
      *text_output << seperator << state.total_immobile.data[i];
    }
    for (int i = 0; i < state.mineral_volume_fraction.size; ++i) {
      *text_output << seperator << state.mineral_volume_fraction.data[i];
    }
    for (int i = 0; i < aux_output.mineral_reaction_rate.size; ++i) {
      *text_output << seperator << aux_output.mineral_reaction_rate.data[i];
    }
    *text_output << std::endl;
  }
}  // end WriteOutput()


/*******************************************************************************
 **
 **  Demo-specific helper rountines for setting up alquimia structs
 **
 *******************************************************************************/
void CopyDemoStateToAlquimiaState(
    const alquimia::drivers::utilities::DemoState& demo_state,
    const AlquimiaProblemMetaData* const alquimia_meta_data,
    AlquimiaState* alquimia_state) {
  alquimia_state->water_density = demo_state.water_density;
  alquimia_state->saturation = demo_state.saturation;
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

void CopyDemoMaterialPropertiesToAlquimiaMaterials(
    const alquimia::drivers::utilities::DemoMaterialProperties& demo_material_props,
    const AlquimiaProblemMetaData& alquimia_meta_data,
    AlquimiaMaterialProperties* alquimia_material_props) {
  alquimia_material_props->volume = demo_material_props.volume;
  
  if (static_cast<size_t>(alquimia_meta_data.isotherm_species_names.size) != 
      demo_material_props.isotherm_species.size()) {
    std::stringstream message;
    message << "ERROR: chemistry engine expects " 
            << alquimia_meta_data.isotherm_species_names.size << " isotherm species. "
            << " but input file only contains " 
            << demo_material_props.isotherm_species.size() << std::endl;
    throw std::runtime_error(message.str());
  }

  char* name;
  name = (char*) calloc(kAlquimiaMaxStringLength, sizeof(char));
  // loop through each species in the *engine's* isotherm list
  for (int i = 0; i < alquimia_meta_data.isotherm_species_names.size; ++i) {
    // save the isotherm species id
    strncpy(name, alquimia_meta_data.isotherm_species_names.data[i],
            kAlquimiaMaxStringLength);
    for (size_t j = 0; j < demo_material_props.isotherm_species.size(); ++j) {
      if (demo_material_props.isotherm_species.at(j).compare(name) == 0) {
        
        alquimia_material_props->isotherm_kd.data[i] = 
            demo_material_props.isotherm_kd.at(j);
        alquimia_material_props->freundlich_n.data[i] = 
            demo_material_props.freundlich_n.at(j);
        alquimia_material_props->langmuir_b.data[i] = 
            demo_material_props.langmuir_b.at(j);
        break;
      }
    }
  }
  free(name);
}  // end CopyDemoMaterialPropertiesToAlquimiaMaterials()


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
    unsigned int max_copy_length = std::min(kAlquimiaMaxStringLength, demo_cond->first.size());
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
                                   demo_aqueous_constraints.at(i).primary_species_name.size());
    std::strncpy(constraint->primary_species_name,
                 demo_aqueous_constraints.at(i).primary_species_name.c_str(),
                 max_copy_length);

    // constraint type
    max_copy_length = std::min(kAlquimiaMaxStringLength,
                               demo_aqueous_constraints.at(i).constraint_type.size());
    std::strncpy(constraint->constraint_type,
                 demo_aqueous_constraints.at(i).constraint_type.c_str(),
                 max_copy_length);

    // associated species
    max_copy_length = std::min(kAlquimiaMaxStringLength,
                               demo_aqueous_constraints.at(i).associated_species.size());
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
                                   demo_mineral_constraints.at(i).mineral_name.size());
    std::strncpy(constraint->mineral_name,
                 demo_mineral_constraints.at(i).mineral_name.c_str(),
                 max_copy_length);
    
    constraint->volume_fraction = demo_mineral_constraints.at(i).volume_fraction;
    constraint->specific_surface_area = demo_mineral_constraints.at(i).specific_surface_area;
  }
}  // end CopyDemoMineralConstraintsToAlquimia()
