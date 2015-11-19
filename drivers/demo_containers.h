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
** 
** Authors: Benjamin Andre <bandre@lbl.gov>
*/

#ifndef ALQUIMIA_DRIVERS_DEMO_CONTAINERS_H_
#define ALQUIMIA_DRIVERS_DEMO_CONTAINERS_H_

/*******************************************************************************
 **
 ** Data transfer containers to move info from the input file to the
 ** demo driver.
 **
 ******************************************************************************/

#include <iostream>
#include <vector>
#include <map>
#include <string>

#include "demo_utils.h"

namespace alquimia {
namespace drivers {
namespace utilities {

struct DemoSimulation {
  std::string description;
  std::string engine;
  std::string engine_inputfile;
  std::string initial_condition;
  double delta_t;
  int num_time_steps;
  std::string time_units;
  void Print(void) const;
};


struct DemoState {
  double water_density;  // [kg/m^3]
  double porosity;  // [-]
  double temperature;  // [celsius]
  double aqueous_pressure; // [Pa]
  std::vector<double> cec;
  std::map<std::string, double> site_density;
  void Print(void) const;
};

struct DemoProperties {
  double volume;  // [m^3]
  double saturation;  // [-]
  std::vector<std::string> isotherm_species;
  std::vector<double> isotherm_kd;  // [?]
  std::vector<double> freundlich_n; // [?]
  std::vector<double> langmuir_b;  // [?]
  void Print(void) const;
};

struct DemoAqueousConstraint {
  std::string primary_species_name;
  std::string constraint_type;
  std::string associated_species;
  double value;
  void Print(void) const;
};
  
struct DemoMineralConstraint {
  std::string mineral_name;
  double volume_fraction;
  double specific_surface_area;
  void Print(void) const;
};

struct DemoGeochemicalCondition {
  std::vector<DemoAqueousConstraint> aqueous_constraints;
  std::vector<DemoMineralConstraint> mineral_constraints;
  void Print(void) const;
};
  
typedef std::map<std::string, DemoGeochemicalCondition> DemoConditions;

void PrintGeochemicalConditions(
    const alquimia::drivers::utilities::DemoConditions& conditions);

}  // namespace utilities
}  // namespace drivers
}  // namespace alquimia
#endif  // ALQUIMIA_DRIVERS_DEMO_CONTAINERS_H_
