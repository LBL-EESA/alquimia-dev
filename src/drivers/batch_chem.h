/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_DRIVERS_CC_BATCH_CHEM_H_
#define ALQUIMIA_DRIVERS_CC_BATCH_CHEM_H_

#include <string>
#include <vector>

#include "alquimia_containers.h"

int CommandLineOptions(int argc, char** argv,
                       std::string* verbosity_name,
                       std::string* input_file_name,
                       std::string* template_file_name,
                       bool* debug_batch_driver);

void SetupChemistryOutput(void);

void SetupTextOutput(const std::string& use_text_output,
                     const std::string& output_time_units,
                     const std::string& input_file_name,
                     std::fstream* text_output, char* time_units,
                     double* time_units_conversion);

void WriteTextOutputHeader(std::fstream* text_output,
                           const char time_units,
                           const std::vector<std::string>& names,
                           const bool using_sorption);

void WriteTextOutput(std::fstream* text_output,
                     const double time,
                     const AlquimiaState_C& state);


#endif  /* ALQUIMIA_DRIVERS_CC_BATCH_CHEM_H_ */
