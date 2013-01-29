/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
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
  std::cout << "  -- State :" << std::endl;
  std::cout << "    density : " << this->water_density << std::endl;
  std::cout << "    saturation : " << this->saturation << std::endl;
  std::cout << "    porosity : " << this->porosity << std::endl;
  std::cout << "    temperature : " << this->temperature << std::endl;
  std::cout << "    pressure : " << this->aqueous_pressure << std::endl;
  std::cout << std::endl;
}  // end Print()


/*******************************************************************************
 **
 **  Material Properties
 **
 *******************************************************************************/
void DemoMaterialProperties::Print(void) const {
  namespace util = alquimia::drivers::utilities;
  std::cout << "  -- Material Properties :" << std::endl;
  std::cout << "    volume : " << this->volume << std::endl;
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
  for (ac = this->aqueous_constraints.begin(); ac != aqueous_constraints.end(); ++ac) {
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
