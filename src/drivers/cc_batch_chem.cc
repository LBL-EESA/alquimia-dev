/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "cc_batch_chem.h"

#ifdef WINDOWS
#include "xgetopt.hh"
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

#include "cc_demo_utils.h"
#include "cfg_reader.h"
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

  util::AlquimiaConfigReader cfg_reader;

  if (!template_file_name.empty()) {
    cfg_reader.WriteTemplateFile(template_file_name);
    exit(EXIT_SUCCESS);
  }

  util::SimulationParameters simulation_params;
  alquimia::AlquimiaState state;
  alquimia::AlquimiaMaterialProperties material_props;
  alquimia::AlquimiaConditions conditions;
  alquimia::AlquimiaMetaData meta_data;
  alquimia::AlquimiaEngineStatus status;

  if (!input_file_name.empty()) {
    cfg_reader.ReadInputFile(input_file_name,
                             &simulation_params, &state, 
                             &material_props, &conditions);
    if (debug_batch_driver) {
      cfg_reader.PrintInput(simulation_params, state, 
                            material_props, conditions);
    }
  }


  double time_units_conversion = 1.0;
  char time_units = 's';
  std::fstream text_output;
  if (simulation_params.use_text_output.size() > 0) {
    SetupTextOutput(simulation_params.use_text_output,
                    simulation_params.output_time_units,
                    input_file_name,
                    &text_output, &time_units, &time_units_conversion);
  }

  alquimia::AlquimiaInterface* chem = NULL;

  try {
    alquimia::AlquimiaInterfaceFactory aif;
    chem = aif.Create(simulation_params.engine);

    chem->Setup(simulation_params.engine_inputfile,
                &meta_data, &status);

    // process the constraints

    // reaction step loop

  } catch (const std::runtime_error& rt_error) {
    std::cout << rt_error.what();
    error = EXIT_FAILURE;
  } catch (const std::logic_error& lg_error) {
    std::cout << lg_error.what();
    error = EXIT_FAILURE;
  }

  if (!error) {
    std::cout << "Success!\n";
  } else {
    std::cout << "Failed!\n";
  }

  text_output.close();
  // cleanup memory
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
                     const alquimia::AlquimiaState& state) {
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


