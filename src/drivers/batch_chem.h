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

void PrintGeochemicalConditions(
    const alquimia::drivers::utilities::DemoConditions& conditions);


void SetupAlquimiaState(
    const alquimia::drivers::utilities::DemoState& demo_state,
    const AlquimiaSizes_C& alquimia_sizes,
    AlquimiaState_C* alquimia_state);

void SetupAlquimiaMetaData(
    const AlquimiaSizes_C& alquimia_sizes,
    AlquimiaMetaData_C* alquimia_meta_data);

void AllocateAlquimiaState(const AlquimiaSizes_C& sizes,
                           AlquimiaState_C* state);
void FreeAlquimiaState(AlquimiaState_C* state);

void AllocateAlquimiaMetaData(const AlquimiaSizes_C& sizes,
                           AlquimiaMetaData_C* meta_data);
void FreeAlquimiaMetaData(const AlquimiaSizes_C& sizes,
                          AlquimiaMetaData_C* metda_data);

void AllocateAlquimiaGeochemicalConditionList(int num_conditions,
    AlquimiaGeochemicalConditionList_C* condition_list);

void AllocateAlquimiaGeochemicalCondition(char* name, int num_constraints,
    AlquimiaGeochemicalCondition_C* condition);

void AllocateAlquimiaGeochemicalConstraint(
    AlquimiaGeochemicalConstraint_C* constraint);

void FreeAlquimiaGeochemicalConditionList(
    AlquimiaGeochemicalConditionList_C* conditions);

void FreeAlquimiaGeochemicalCondition(
    AlquimiaGeochemicalCondition_C* condition);

void FreeAlquimiaGeochemicalConstraint(
    AlquimiaGeochemicalConstraint_C* constraint);

void PrintAlquimiaSizes(const AlquimiaSizes_C& sizes);
void PrintAlquimiaMetaData(const AlquimiaSizes_C& sizes,
                           const AlquimiaMetaData_C& meta_data);
void PrintAlquimiaState(const AlquimiaSizes_C& sizes,
                        const AlquimiaState_C& state);


void PrintAlquimiaGeochemicalConditionList(
    AlquimiaGeochemicalConditionList_C* condition_list);
void PrintAlquimiaGeochemicalCondition(
    AlquimiaGeochemicalCondition_C* condition);
void PrintAlquimiaGeochemicalConstraint(
    AlquimiaGeochemicalConstraint_C* constraint);

#endif  /* ALQUIMIA_DRIVERS_CC_BATCH_CHEM_H_ */
