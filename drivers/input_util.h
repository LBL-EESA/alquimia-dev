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

#ifndef ALQUIMIA_INPUT_UTILS_H
#define ALQUIMIA_INPUT_UTILS_H

#include "alquimia/alquimia_containers.h"

// This file contains functions that can be helpful in parsing input from an 
// ini_handler.

// Returns true if the given section name corresponds to a state, and 
// if so, fills state_name with the name of the state.
bool IsStateSection(const char* section, char* state_name);

// Parses the given input value into the given state container. Problem 
// metadata is used to properly sort values associated with species, etc.
void ParseStateInput(const char* section, 
                     const char* name,
                     const char* value,
                     AlquimiaProblemMetaData* metadata,
                     AlquimiaState* state);

// Returns true if the given section name corresponds to a (material) property,
// and if so, fills properties_name with the name of the set of properties.
bool IsPropertiesSection(const char* section, char* properties_name);

// Parses the given input value into the given properties container. Problem
// metadata is used to properly sort values associated with species, etc.
void ParsePropertyInput(const char* section,
                        const char* name,
                        const char* value,
                        AlquimiaProblemMetaData* metadata,
                        AlquimiaProperties* properties);

// Returns true if the given section name corresponds to a geochemical condition,
// and if so, fills condition_name with the name of the geochemical condition.
bool IsGeochemicalConditionSection(const char* section, char* condtion_name);

// Parses the given input value into the given properties container.
void ParseGeochemicalConditionInput(const char* name,
                                    const char* value,
                                    AlquimiaGeochemicalCondition* condition);

// Returns true if the given section name corresponds to an aqueous constraint,
// and if so, fills constraint_name with the name of the constraint.
bool IsAqueousConstraintSection(const char* section, char* constraint_name);

// Parses the given input value into the given aqueous constraint.
void ParseAqueousConstraintInput(const char* primary_species,
                                 const char* text,
                                 AlquimiaAqueousConstraint* constraint);

// Returns true if the given section name corresponds to a mineral constraint, and 
// if so, fills constraint_name with the name of the constraint.
bool IsMineralConstraintSection(const char* section, char* constraint_name);

// Parses the given input value into the given aqueous constraint.
void ParseMineralConstraintInput(const char* mineral_name,
                                 const char* text,
                                 AlquimiaMineralConstraint* constraint);

#endif 
