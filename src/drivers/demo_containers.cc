/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

//
// Alquimia Copyright (c) 2013, The Regents of the University of California, 
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
// Authors: Benjamin Andre <bandre@lbl.gov>
//

#include "demo_containers.h"

#include <iostream>

namespace alquimia {
namespace drivers {
namespace utilities {

/*******************************************************************************
 **
 **  Simulation
 **
 *******************************************************************************/
void DemoSimulation::Print(void) const {
  std::cout << "  -- Simulation parameters :" << std::endl;
  std::cout << "    description : " << this->description << std::endl;
  std::cout << "    engine : " << this->engine << std::endl;
  std::cout << "    engine inputfile : " << this->engine_inputfile << std::endl;
  std::cout << "    delta t : " << this->delta_t 
            << " [" << this->time_units << "]" << std::endl;
  std::cout << "    num times steps : " << this->num_time_steps << std::endl;
  std::cout << "    initial condition : " << this->initial_condition << std::endl;
  std::cout << std::endl;

}  // end Print()

/*******************************************************************************
 **
 **  State
 **
 *******************************************************************************/

void DemoState::Print(void) const {
  namespace util = alquimia::drivers::utilities;
  std::cout << "  -- State :" << std::endl;
  std::cout << "    density : " << this->water_density << std::endl;
  std::cout << "    porosity : " << this->porosity << std::endl;
  std::cout << "    temperature : " << this->temperature << std::endl;
  std::cout << "    pressure : " << this->aqueous_pressure << std::endl;
  util::PrintVector("    cec", this->cec);
  util::PrintMap("    site density", this->site_density);
  std::cout << std::endl;
}  // end Print()


/*******************************************************************************
 **
 **  Material Properties
 **
 *******************************************************************************/
void DemoProperties::Print(void) const {
  namespace util = alquimia::drivers::utilities;
  std::cout << "  -- Material Properties :" << std::endl;
  std::cout << "    volume : " << this->volume << std::endl;
  std::cout << "    saturation : " << this->saturation << std::endl;
  util::PrintVector("    isotherm_species", this->isotherm_species);
  util::PrintVector("    isotherm_kd", this->isotherm_kd);
  util::PrintVector("    freundlich_n", this->freundlich_n);
  util::PrintVector("    langmuir_b", this->langmuir_b);
  std::cout << std::endl;
}  // end Print()

/*******************************************************************************
 **
 **  Geochemical Constraints
 **
 *******************************************************************************/
void DemoAqueousConstraint::Print(void) const {
  std::cout << "        " << this->primary_species_name << std::endl;
  std::cout << "            type : " << this->constraint_type << std::endl;
  std::cout << "            associated : " << this->associated_species << std::endl;
  std::cout << "            value : " << this->value << std::endl;
}  // end Print()

void DemoMineralConstraint::Print(void) const {
  std::cout << "        " << this->mineral_name << std::endl;
  std::cout << "            volume fraction : " << this->volume_fraction << std::endl;
  std::cout << "            specific surface area : " << this->specific_surface_area << std::endl;
}  // end Print()

/*******************************************************************************
 **
 **  Geochemical Conditions
 **
 *******************************************************************************/

void DemoGeochemicalCondition::Print(void) const {
  
  std::vector<DemoAqueousConstraint>::const_iterator ac;
  for (ac = this->aqueous_constraints.begin(); ac != this->aqueous_constraints.end(); ++ac) {
    ac->Print();
  }

  std::vector<DemoMineralConstraint>::const_iterator mc;
  for (mc = this->mineral_constraints.begin(); mc != this->mineral_constraints.end(); ++mc) {
    mc->Print();
  }

}  // end Print()

void PrintGeochemicalConditions(
    const alquimia::drivers::utilities::DemoConditions& conditions)
{
  namespace util = alquimia::drivers::utilities;
  std::cout << "  -- Geochemical Conditions :" << std::endl;
  for (util::DemoConditions::const_iterator c = conditions.begin();
       c != conditions.end(); ++c) {
    std::cout << "    " << c->first << " : " << std::endl;
    c->second.Print();
  }
}  // end PrintGeochemicalConditions()

}  // namespace utilities
}  // namespace drivers
}  // namespace alquimia
