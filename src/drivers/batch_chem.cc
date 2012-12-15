/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "batch_chem.h"

#ifdef WINDOWS
#include "xgetopt.h"
#else
#include <unistd.h>
#endif

#include <cstdlib>
#include <cctype>

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

#include "alquimia_interface_factory.h"
#include "alquimia_interface.h"
#include "alquimia_containers.h"
#include "alquimia_strings.h"

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
    cfg_reader.ReadInputFile(input_file_name,
                             &demo_simulation, &demo_state, 
                             &demo_material_props, &demo_conditions);
    if (debug_batch_driver) {
      std::cout << "- Input File ---------------------------------------------------------\n";
      demo_simulation.Print();
      demo_state.Print();
      demo_material_props.Print();
      PrintGeochemicalConditions(demo_conditions);
      std::cout << "--------------------------------------------------------- Input File -\n";
    }
  }


  double time_units_conversion = 1.0;
  char time_units = 's';
  std::fstream text_output;
  if (demo_simulation.use_text_output.size() > 0) {
    SetupTextOutput(demo_simulation.use_text_output,
                    demo_simulation.output_time_units,
                    input_file_name,
                    &text_output, &time_units, &time_units_conversion);
  }

  alquimia::AlquimiaInterface* chem = NULL;
  AlquimiaSizes_C alquimia_sizes;
  AlquimiaState_C alquimia_state;
  AlquimiaMaterialProperties_C alquimia_material_props;
  AlquimiaAuxiliaryData_C alquimia_aux_data;
  AlquimiaEngineStatus_C alquimia_status;
  AlquimiaMetaData_C alquimia_meta_data;
  AlquimiaGeochemicalCondition_C alquimia_conditions;
  AlquimiaOutputData_C alquimia_output_data;

  try {
    alquimia::AlquimiaInterfaceFactory aif;
    chem = aif.Create(demo_simulation.engine);

    chem->Setup(demo_simulation.engine_inputfile, &alquimia_sizes);
    if (true) {
      PrintAlquimiaSizes(alquimia_sizes);
    }

    SetupAlquimiaMetaData(alquimia_sizes, &alquimia_meta_data);
    PrintAlquimiaMetaData(alquimia_sizes, alquimia_meta_data);
    chem->GetEngineMetaData(&alquimia_sizes, &alquimia_meta_data);
    if (true) {
      PrintAlquimiaMetaData(alquimia_sizes, alquimia_meta_data);
    }

    std::cout << "-- Setting up alquimia state container...\n";
    SetupAlquimiaState(demo_state, alquimia_sizes, &alquimia_state);
    PrintAlquimiaState(alquimia_sizes, alquimia_state);

    // process the geochemical conditions
    AllocateAlquimiaGeochemicalConditions(demo_conditions.size(),
                                          &alquimia_conditions);
    for (util::DemoConditions::iterator condition = demo_conditions.begin();
         condition != demo_conditions.end(); ++condition) {
      std::cout << "    " << condition->first << " : " << std::endl;
      for (util::DemoGeochemicalCondition::const_iterator g = condition->second.begin();
           g != condition->second.end(); ++g) {
        g->Print();
      }
      
      //chem->ProcessCondition();
    }
    // reaction step loop

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
  FreeAlquimiaMetaData(alquimia_sizes, &alquimia_meta_data);
  FreeAlquimiaState(&alquimia_state);
  FreeAlquimiaGeochemicalConditions(&alquimia_conditions);
  delete chem;
  
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
void SetupTextOutput(const std::string& use_text_output,
                     const std::string& output_time_units,
                     const std::string& input_file_name,
                     std::fstream* text_output, char* time_units,
                     double* time_units_conversion) {
  namespace util = alquimia::drivers::utilities;
  // are we writting to observations to a text file?
  if (util::CaseInsensitiveStringCompare(use_text_output, "true") ||
      util::CaseInsensitiveStringCompare(use_text_output, "yes") ||
      util::CaseInsensitiveStringCompare(use_text_output, "on")) {
    // generate the output file name:
    size_t position = input_file_name.find_last_of('.');
    std::string text_output_name = input_file_name.substr(0, position) + ".txt";

    text_output->open(text_output_name.c_str(), std::fstream::out);

    // do we want to change the time units for the output?
    if (output_time_units.size() > 0) {
      *time_units = std::tolower(output_time_units.at(0));
      switch (*time_units) {
        case 's':
          break;
        case 'm':
          *time_units_conversion = 60.0;
          break;
        case 'h':
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
    *time_units_conversion = 1.0 / (*time_units_conversion);
  }

}  // end SetupTextOutput()

void WriteTextOutputHeader(std::fstream* text_output, const char time_units,
                           const std::vector<std::string>& names,
                           const bool using_sorption) {
  if (text_output->is_open()) {
    *text_output << "# Time(" << time_units << ")";
    for (std::vector<std::string>::const_iterator name = names.begin();
         name != names.end(); ++name) {
      *text_output <<  " , " << *name;
    }
    if (using_sorption) {
      for (std::vector<std::string>::const_iterator name = names.begin();
           name != names.end(); ++name) {
        *text_output <<  " , " << *name << "_sorbed";
      }
    }
    *text_output << std::endl;
  }
}  // end WriteTextOutputHeader()

void WriteTextOutput(std::fstream* text_output, const double time, 
                     const AlquimiaState_C& state) {
  if (text_output->is_open()) {
    // std::string seperator(" , ");
    // *text_output << std::scientific << std::setprecision(6) << std::setw(15) << time;
    // for (int i = 0; i < state.total.size(); ++i) {
    //   *text_output << seperator << state.total.at(i);
    // }
    // for (int i = 0; i < state.total_sorbed.size(); ++i) {
    //   *text_output << seperator << state.total_sorbed.at(i);
    // }
    // *text_output << std::endl;
  }
}  // end WriteTextOutput()


void PrintGeochemicalConditions(
    const alquimia::drivers::utilities::DemoConditions& conditions)
{
  namespace util = alquimia::drivers::utilities;
  std::cout << "  -- Geochemical Conditions :" << std::endl;
  for (util::DemoConditions::const_iterator c = conditions.begin();
       c != conditions.end(); ++c) {
    std::cout << "    " << c->first << " : " << std::endl;
    for (util::DemoGeochemicalCondition::const_iterator g = c->second.begin();
         g != c->second.end(); ++g) {
      g->Print();
    }
  }
}  // end PrintGeochemicalConditions()

/*******************************************************************************
 **
 **  Demo-specific helper rountines for setting up alquimia structs
 **
 *******************************************************************************/
void SetupAlquimiaState(
    const alquimia::drivers::utilities::DemoState& demo_state,
    const AlquimiaSizes_C& alquimia_sizes,
    AlquimiaState_C* alquimia_state) {
  alquimia_state->water_density = demo_state.water_density;
  alquimia_state->saturation = demo_state.saturation;
  alquimia_state->porosity = demo_state.porosity;
  alquimia_state->temperature = demo_state.temperature;
  alquimia_state->aqueous_pressure = demo_state.aqueous_pressure;
  AllocateAlquimiaState(alquimia_sizes, alquimia_state);
  
}  // end SetupAlquimiaState()


void SetupAlquimiaMetaData(
    const AlquimiaSizes_C& alquimia_sizes,
    AlquimiaMetaData_C* alquimia_meta_data) {
  alquimia_meta_data->thread_safe = false;
  alquimia_meta_data->temperature_dependent = false;
  alquimia_meta_data->pressure_dependent = false;
  alquimia_meta_data->porosity_update = false;
  AllocateAlquimiaMetaData(alquimia_sizes, alquimia_meta_data);
  
}  // end SetupAlquimiaState()


/*******************************************************************************
 **
 **  Helper routines that may eventually be part of alquimia library
 **
 *******************************************************************************/

void AllocateAlquimiaState(const AlquimiaSizes_C& sizes,
                           AlquimiaState_C* state) {
  state->total_primary = NULL;
  state->total_sorbed = NULL;
  state->free_ion = NULL;
  state->mineral_volume_fraction = NULL;
  state->mineral_specific_surface_area = NULL;
  state->cation_exchange_capacity = NULL;
  state->surface_site_density = NULL;
  //state-> = NULL;

  if (sizes.num_primary > 0) {
    state->total_primary = (double*) calloc(sizes.num_primary, sizeof(double));
    if (NULL == state->total_primary) {
      // TODO(bja): error handling
    }
    state->free_ion = (double*) calloc(sizes.num_primary, sizeof(double));
    if (NULL == state->free_ion) {
      // TODO(bja): error handling
    }
  }

  if (sizes.num_kinetic_minerals > 0) {
    state->mineral_volume_fraction = (double*) calloc(
        sizes.num_kinetic_minerals, sizeof(double));
    if (NULL == state->mineral_volume_fraction) {
      // TODO(bja): error handling
    }
    state->mineral_specific_surface_area = (double*) calloc(
        sizes.num_kinetic_minerals, sizeof(double));
    if (NULL == state->mineral_specific_surface_area) {
      // TODO(bja): error handling
    }
  }
}  // end AllocateAlquimiaState()

void FreeAlquimiaState(AlquimiaState_C* state) {
  if (state != NULL) {
    free(state->total_primary);
    free(state->free_ion);
    free(state->mineral_volume_fraction);
    free(state->mineral_specific_surface_area);
  }
  state->total_primary = NULL;
  state->free_ion = NULL;
  state->mineral_volume_fraction = NULL;
  state->mineral_specific_surface_area = NULL;
}  // end FreeAlquimiaState()

void AllocateAlquimiaMetaData(const AlquimiaSizes_C& sizes,
                           AlquimiaMetaData_C* meta_data) {
  meta_data->primary_indices = NULL;
  meta_data->primary_names = NULL;
  //meta_data-> = NULL;

  if (sizes.num_primary > 0) {
    meta_data->primary_indices = (int*) calloc(sizes.num_primary, sizeof(int));
    if (NULL == meta_data->primary_indices) {
      // TODO(bja): error handling
    }
    meta_data->primary_names = (char**) calloc(sizes.num_primary, sizeof(char*));
    if (NULL == meta_data->primary_names) {
      // TODO(bja): error handling
    }
    for (int i = 0; i < sizes.num_primary; ++i) {
      meta_data->primary_names[i] = (char*) calloc(ALQUIMIA_MAX_STRING_LENGTH, sizeof(char));
      if (NULL == meta_data->primary_names[i]) {
        // TODO(bja): error handling
      }
    }
  }
}  // end AllocateAlquimiaMetaData()

void FreeAlquimiaMetaData(const AlquimiaSizes_C& sizes,
                          AlquimiaMetaData_C* meta_data) {
  if (meta_data != NULL) {
    free(meta_data->primary_indices);
    for (int i = 0; i < sizes.num_primary; ++i) {
      free(meta_data->primary_names[i]);
    }
    free(meta_data->primary_names);
  }
  meta_data->primary_indices = NULL;
  meta_data->primary_names = NULL;
}  // end FreeAlquimiaMetaData()



void AllocateAlquimiaGeochemicalConditions(
    const int num_conditions,
    AlquimiaGeochemicalCondition_C* conditions) {
  conditions = NULL;
  if (num_conditions > 0) {
  }
}  // end AllocateAlquimiaGeochemicalConditions()

void FreeAlquimiaGeochemicalConditions(
    AlquimiaGeochemicalCondition_C* condition) {
  free(condition);
  condition = NULL;
}  // end FreeAlquimiaState()

void PrintAlquimiaSizes(const AlquimiaSizes_C& sizes) {
  std::cout << "  sizes :\n"
            << "    num primary species : " << sizes.num_primary << "\n"
            << "    num kinetic minerals : " << sizes.num_kinetic_minerals << "\n"
            << "    num aqueous complexes : " << sizes.num_aqueous_complexes << "\n"
            << "    num surface sites : " << sizes.num_surface_sites << "\n"
            << "    num ion exchange sites : " << sizes.num_ion_exchange_sites << "\n";
}  // end PrintAlquimiaSizes()

void PrintAlquimiaMetaData(const AlquimiaSizes_C& sizes,
                           const AlquimiaMetaData_C& meta_data) {
  std::cout << "  meta data :\n"
            << "    thread_safe : " << meta_data.thread_safe << "\n"
            << "    temperature_dependent : " << meta_data.temperature_dependent << "\n"
            << "    pressure_dependent : " << meta_data.pressure_dependent << "\n"
            << "    porosity_update  : " << meta_data.porosity_update << "\n";
  std::cout << "    primary indices (" << meta_data.primary_indices << ") (" 
            << sizes.num_primary << ") :\n      ";
  for (int i = 0; i < sizes.num_primary; ++i) {
    std::cout << meta_data.primary_indices[i] << ", ";
  }
  std::cout << "\n";
  std::cout << "    primary names (" << meta_data.primary_names << ") (" 
            << sizes.num_primary << ") :\n      ";
  for (int i = 0; i < sizes.num_primary; ++i) {
    std::cout << meta_data.primary_names[i] << ", ";
  }
  std::cout << "\n";
}  // end PrintAlquimiaMetaData()

void PrintAlquimiaState(const AlquimiaSizes_C& sizes,
                        const AlquimiaState_C& state) {
  std::cout << "  state:\n"
            << "    water density : " << state.water_density << "\n"
            << "    saturation : " << state.saturation << "\n"
            << "    porosity : " << state.porosity << "\n"
            << "    temperature : " << state.temperature << "\n"
            << "    aqueous_pressure : " << state.aqueous_pressure << "\n";

  std::cout << "    total_primary (" << sizes.num_primary << "):\n      ";
  for (int i = 0; i < sizes.num_primary; ++i) {
    std::cout << state.total_primary[i] << ", ";
  }
  std::cout << "\n";

  std::cout << "    free_ion (" << sizes.num_primary << "):\n      ";
  for (int i = 0; i < sizes.num_primary; ++i) {
    std::cout << state.free_ion[i] << ", ";
  }
  std::cout << "\n";

  std::cout << "    mineral volume fractions (" << sizes.num_kinetic_minerals << "):\n      ";
  for (int i = 0; i < sizes.num_kinetic_minerals; ++i) {
    std::cout << state.mineral_volume_fraction[i] << ", ";
  }
  std::cout << "\n";

  std::cout << "    mineral specific surface area (" << sizes.num_kinetic_minerals << "):\n      ";
  for (int i = 0; i < sizes.num_kinetic_minerals; ++i) {
    std::cout << state.mineral_specific_surface_area[i] << ", ";
  }
  std::cout << "\n";
}
