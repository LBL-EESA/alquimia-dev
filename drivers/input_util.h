/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
** Alquimia Copyright (c) 2013-2015, The Regents of the University of California, 
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
*/

#ifndef ALQUIMIA_INPUT_UTIL_H
#define ALQUIMIA_INPUT_UTIL_H

#include "alquimia/alquimia_interface.h"

// This file contains functions that can be helpful in parsing input from an 
// ini_handler.

//------------------------------------------------------------------------
//                      High-level parsing utilities
//------------------------------------------------------------------------
// These can be used to perform sweeps of the input files to gather 
// information, and rely on the low-level parsing utilities below.
//------------------------------------------------------------------------

// Creates an Alquimia interface using the information in the [chemistry]
// section of the input file, initializing the interface and the other given 
// containers. Note that engine_status must be allocated before this call.
void Input_CreateAlquimiaInterface(const char* input_file,
                                   AlquimiaInterface* engine_interface,
                                   AlquimiaSizes* engine_sizes,
                                   AlquimiaEngineFunctionality* engine_functionality,
                                   AlquimiaEngineStatus* engine_status);

// Reads the names of regions from the input file.
void Input_GetRegions(const char* input_file,
                      AlquimiaVectorString* region_names);

// Reads the data (state, properties, initial condition, list of cells) for 
// the region with the given name from the input file.
void Input_GetRegionData(const char* input_file,
                         const char* region_name,
                         AlquimiaProblemMetaData* problem_metadata,
                         AlquimiaState* region_state,
                         AlquimiaProperties* region_properties,
                         AlquimiaGeochemicalCondition* region_initial_condition,
                         AlquimiaVectorInt* region_cells);

// Reads any geochemical conditions specified in the file that are not 
// initial conditions for regions.
void Input_GetGeochemicalConditions(const char* input_file,
                                    AlquimiaGeochemicalConditionVector* conditions);

// Reads any output options/parameters in the file, storing them in the arguments.
// output_type is currently set to "python" or "gnuplot".
void Input_GetOutputParameters(const char* input_file,
                               char* output_type,
                               char* output_file,
                               bool* verbose);

//------------------------------------------------------------------------
//                      Low-level parsing utilities
//------------------------------------------------------------------------
// These can be used within the INIH parser to retrieve data from specific
// sections.
//------------------------------------------------------------------------

// Returns true if the given section name corresponds to a state, and 
// if so, fills state_name with the name of the state.
bool Input_IsStateSection(const char* section, char* state_name);

// Parses the given input value into the given state container. Problem 
// metadata is used to properly sort values associated with species, etc.
void Input_ParseState(const char* section, 
                      const char* name,
                      const char* value,
                      AlquimiaProblemMetaData* metadata,
                      AlquimiaState* state);

// Returns true if the given section name corresponds to a (material) property,
// and if so, fills properties_name with the name of the set of properties.
bool Input_IsPropertiesSection(const char* section, char* properties_name);

// Parses the given input value into the given properties container. Problem
// metadata is used to properly sort values associated with species, etc.
void Input_ParseProperty(const char* section,
                         const char* name,
                         const char* value,
                         AlquimiaProblemMetaData* metadata,
                         AlquimiaProperties* properties);

// Returns true if the given section name corresponds to a geochemical condition,
// and if so, fills condition_name with the name of the geochemical condition.
bool Input_IsGeochemicalConditionSection(const char* section, char* condtion_name);

// Parses the given input value into the given properties container.
void Input_ParseGeochemicalCondition(const char* name,
                                     const char* value,
                                     AlquimiaGeochemicalCondition* condition);

// Returns true if the given section name corresponds to an aqueous constraint,
// and if so, fills constraint_name with the name of the constraint.
bool Input_IsAqueousConstraintSection(const char* section, char* constraint_name);

// Parses the given input value into the given aqueous constraint.
void Input_ParseAqueousConstraint(const char* primary_species,
                                  const char* text,
                                  AlquimiaAqueousConstraint* constraint);

// Returns true if the given section name corresponds to a mineral constraint, and 
// if so, fills constraint_name with the name of the constraint.
bool Input_IsMineralConstraintSection(const char* section, char* constraint_name);

// Parses the given input value into the given aqueous constraint.
void Input_ParseMineralConstraint(const char* mineral_name,
                                  const char* text,
                                  AlquimiaMineralConstraint* constraint);

#endif 
