/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
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
  double saturation;  // [-]
  double porosity;  // [-]
  double temperature;  // [celsius]
  double aqueous_pressure; // [Pa]

  void Print(void) const;
};

struct DemoMaterialProperties {
  double volume;  // [m^3]
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
