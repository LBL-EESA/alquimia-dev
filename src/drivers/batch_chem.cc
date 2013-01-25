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

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

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
  // Create the alquimia structures.  NOTE: These are for a single
  // grid cell. For openmp/threaded code, you'd need one per thead,
  // but we don't know if the engine is thread safe until after the
  // call to GetMetaData()!
  //
  struct AlquimiaInterface chem;
  AlquimiaEngineStatus chem_status;
  AlquimiaData chem_data;
  AlquimiaGeochemicalConditionList alquimia_conditions;

  try {
    // All alquimia function calls require a status object.
    AllocateAlquimiaEngineStatus(&chem_status);
    // Create the chemistry engine
    AllocateAlquimiaInterface(&chem);
    CreateAlquimiaInterface(demo_simulation.engine.c_str(), &chem, &chem_status);
    if (chem_status.error != 0) {
      std::cout << chem_status.message << std::endl;
      return chem_status.error;
    }

    // setup the engine and get the memory requirements
    chem.Setup(demo_simulation.engine_inputfile.c_str(),
               chem.engine_state,
               &chem_data.sizes,
               &chem_status);
    if (chem_status.error != 0) {
      std::cout << chem_status.message << std::endl;
      PrintAlquimiaSizes(&chem_data.sizes);
      return chem_status.error;
    }

    // allocate memory for alquimia data transfer containers.
    AllocateAlquimiaData(&chem_data);

    // allocate the remaining memory in the driver (mesh dependent)

    // get the driver meta data (thread safe, temperature/pressure
    // dependent, species names, etc)
    chem.GetEngineMetaData(chem.engine_state,
                           &chem_data.meta_data,
                           &chem_status);
    if (chem_status.error != 0) {
      std::cout << chem_status.message << std::endl;
      PrintAlquimiaMetaData(&chem_data.meta_data);
      return chem_status.error;
    }
    PrintAlquimiaMetaData(&chem_data.meta_data);

    // finish initializing the driver, e.g. openmp for thread safe
    // engines, verify material properties, etc

    //
    // prepare for constraint processing
    //
    AllocateAlquimiaGeochemicalConditionList(demo_conditions.size(),
                                             &alquimia_conditions);

    // initialize the alquimia state and material properties with
    // appropriate values from the driver's memory.
    CopyDemoStateToAlquimiaState(demo_state, &chem_data.state);
    CopyDemoMaterialPropertiesToAlquimiaMaterials(
        demo_material_props, &chem_data.material_properties);

    // Read the geochemical conditions from the driver's native format
    // and store them in alquimia's format
    CopyDemoConditionsToAlquimiaConditions(demo_conditions, &alquimia_conditions);

    //PrintAlquimiaState(&chem_data.state);
    //PrintAlquimiaGeochemicalConditionList(&alquimia_conditions);

    for (int i = 0; i < alquimia_conditions.num_conditions; ++i) {
      // ask the engine to process the geochemical conditions
      if (demo_simulation.initial_condition.compare(
              alquimia_conditions.conditions[i].name) == 0) {
        // for batch, only care about the IC. If the conditions have
        // spatially dependent values, then the state and material
        // properties need to be updated here!
        chem.ProcessCondition(chem.engine_state,
                              &(alquimia_conditions.conditions[i]),
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
      }
      // store the processed geochemical conditions in driver's memory
    }
    // we are done with the conditions (data is stored in driver state
    // vectors at this point).
    FreeAlquimiaGeochemicalConditionList(&alquimia_conditions);

    char time_units;
    double time_units_conversion;  // [time_units / sec]
    SetTimeUnits(demo_simulation.time_units,
                 &time_units, &time_units_conversion);

    // set delta t: [time_units]/[time_units/sec] = [sec]
    double delta_t = demo_simulation.delta_t / time_units_conversion;
    double time = 0.0;
    // get initial pH
    chem.GetAuxiliaryOutput(chem.engine_state,
                            &chem_data.material_properties,
                            &chem_data.state,
                            &chem_data.aux_data,
                            &chem_data.aux_output,
                            &chem_status);
    // save the IC to our output file
    WriteOutputHeader(&text_output, time_units, chem_data.meta_data);
    WriteOutput(&text_output, time, chem_data.state, chem_data.aux_output);
    std::cout << "Starting reaction stepping with dt = " << delta_t << " [s]\n";
    for (int t = 0; t < demo_simulation.num_time_steps; ++t) {
      time += delta_t;
      std::cout << ".";
      //std::cout << "reaction step : " << t << "  time: " << time << std::endl;
      if (false) {
        PrintAlquimiaState(&chem_data.state);
        PrintAlquimiaAuxiliaryData(&chem_data.aux_data);
      }
      // unpack from driver memory, since this is batch, no unpacking
      chem.ReactionStepOperatorSplit(chem.engine_state,
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
      chem.GetAuxiliaryOutput(chem.engine_state,
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
      WriteOutput(&text_output, out_time, chem_data.state, chem_data.aux_output);
      // repack into driver memory...
    }
    std::cout << std::endl;
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

  if (!error) {
    std::cout << "Success!\n";
  } else {
    std::cout << "Failed!\n";
  }

  text_output.close();

  // cleanup memory
  chem.Shutdown(chem.engine_state, &chem_status);
  FreeAlquimiaInterface(&chem);
  FreeAlquimiaData(&chem_data);
  FreeAlquimiaEngineStatus(&chem_status);

  error = PetscFinalize();

  return error;
}  // end main()

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
                       const AlquimiaMetaData& meta_data) {
  if (text_output->is_open()) {
    *text_output << "# \"Time [" << time_units << "]\"";
    *text_output << " , \"pH\"";
    for (int i = 0; i < meta_data.primary_names.size; ++i) {
      *text_output <<  " , \"Total " << meta_data.primary_names.data[i] << " [M]\"";
    }
    // TODO: sorbed header...
    for (int i = 0; i < meta_data.mineral_names.size; ++i) {
      *text_output << " , \"" << meta_data.mineral_names.data[i] << " VF\"";
    }
    for (int i = 0; i < meta_data.mineral_names.size; ++i) {
      *text_output << " , \"" << meta_data.mineral_names.data[i] << " Rate [moles/sec]\"";
    }
    *text_output << std::endl;
  }
}  // end WriteOutputHeader()

void WriteOutput(std::fstream* text_output, const double time,
                 const AlquimiaState& state,
                 const AlquimiaAuxiliaryOutputData& aux_output) {
  if (text_output->is_open()) {
    std::string seperator(" , ");
    *text_output << std::scientific << std::setprecision(6) << std::setw(15) << time;
    *text_output  << seperator << aux_output.pH;
    for (int i = 0; i < state.total_primary.size; ++i) {
      *text_output << seperator << state.total_primary.data[i];
    }
    for (int i = 0; i < state.total_sorbed.size; ++i) {
      *text_output << seperator << state.total_sorbed.data[i];
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
    AlquimiaState* alquimia_state) {
  alquimia_state->water_density = demo_state.water_density;
  alquimia_state->saturation = demo_state.saturation;
  alquimia_state->porosity = demo_state.porosity;
  alquimia_state->temperature = demo_state.temperature;
  alquimia_state->aqueous_pressure = demo_state.aqueous_pressure;
}  // end CopyDemoStateToAlquimiaState()

void CopyDemoMaterialPropertiesToAlquimiaMaterials(
    const alquimia::drivers::utilities::DemoMaterialProperties& demo_material_props,
    AlquimiaMaterialProperties* alquimia_material_props) {
  alquimia_material_props->volume = demo_material_props.volume;
  // TODO(bja) : loop through and copy vector based material properties!
}  // end CopyDemoMaterialPropertiesToAlquimiaMaterials()


void CopyDemoConditionsToAlquimiaConditions(
    const alquimia::drivers::utilities::DemoConditions& demo_conditions,
    AlquimiaGeochemicalConditionList* alquimia_conditions) {
  namespace util = alquimia::drivers::utilities;

  // copy the geochemical conditions
  util::DemoConditions::const_iterator demo_cond;
  unsigned int i_cond;
  for (demo_cond = demo_conditions.begin(), i_cond = 0;
       i_cond < demo_conditions.size(); ++i_cond, ++demo_cond) {
    std::cout << "    " << demo_cond->first << " : " << i_cond << std::endl;
    AlquimiaGeochemicalCondition* condition =
        &(alquimia_conditions->conditions[i_cond]);
    // std::string.c_str() returns a const char*, so we need to copy
    // it to our own memory.
    char* condition_name = new char [demo_cond->first.size() + 1];
    strcpy(condition_name, demo_cond->first.c_str());
    AllocateAlquimiaGeochemicalCondition(condition_name,
                                         demo_cond->second.size(),
                                         condition);
    delete condition_name;
    for (unsigned int i_const = 0; i_const < demo_cond->second.size(); ++i_const) {
      std::cout << "    " << demo_cond->first << " : " << i_cond << " : "
                << i_const << std::endl;
      AlquimiaGeochemicalConstraint* constraint =
          &(alquimia_conditions->conditions[i_cond].constraints[i_const]);
      AllocateAlquimiaGeochemicalConstraint(constraint);
      // copy demo constraint to alquimia constraint
      int max_copy_length = kAlquimiaMaxStringLength;
      if (strlen(demo_cond->second[i_const].primary_species.c_str()) <
          kAlquimiaMaxStringLength) {
        max_copy_length = strlen(demo_cond->second[i_const].primary_species.c_str());
      }
      std::strncpy(constraint->primary_species,
                   demo_cond->second[i_const].primary_species.c_str(),
                   max_copy_length);

      max_copy_length = kAlquimiaMaxStringLength;
      if (strlen(demo_cond->second[i_const].constraint_type.c_str()) <
          kAlquimiaMaxStringLength) {
        max_copy_length = strlen(demo_cond->second[i_const].constraint_type.c_str());
      }
      std::strncpy(constraint->constraint_type,
                   demo_cond->second[i_const].constraint_type.c_str(),
                   max_copy_length);

      max_copy_length = kAlquimiaMaxStringLength;
      if (strlen(demo_cond->second[i_const].associated_species.c_str()) <
          kAlquimiaMaxStringLength) {
        max_copy_length = strlen(demo_cond->second[i_const].associated_species.c_str());
      }
      std::strncpy(constraint->associated_species,
                   demo_cond->second[i_const].associated_species.c_str(),
                   max_copy_length);
      constraint->value = demo_cond->second[i_const].value;
    }
  }
}  // end CopyDemoConditionsToAlquimiaConditions()
