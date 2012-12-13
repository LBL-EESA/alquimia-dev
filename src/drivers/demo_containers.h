/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_DRIVERS_DEMO_CONTAINERS_H_
#define ALQUIMIA_DRIVERS_DEMO_CONTAINERS_H_

/*******************************************************************************
 **
 ** Data transfer containers to move info from the input file to the
 ** demo driver.
 **
 ******************************************************************************/

#include <vector>
#include <map>
#include <string>

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
  std::string use_text_output;
  std::string output_time_units;
};


struct DemoState {
  double water_density;  // [kg/m^3]
  double saturation;  // [-]
  double porosity;  // [-]
  double temperature;  // [celsius]
  double aqueous_pressure; // [Pa]
};

struct DemoMaterialProperties {
  double volume;  // [m^3]
  std::vector<double> isotherm_kd;  // [?]
  std::vector<double> freundlich_n; // [?]
  std::vector<double> langmuir_b;  // [?]
};

struct DemoGeochemicalConstraint {
  std::string primary_species;
  std::string constraint_type;
  std::string associated_species;
  std::string value;
};
  
typedef std::vector<DemoGeochemicalConstraint> DemoGeochemicalCondition;

typedef std::map<std::string, DemoGeochemicalCondition> DemoConditions;


}  // namespace utilities
}  // namespace drivers
}  // namespace alquimia
#endif  // ALQUIMIA_DRIVERS_DEMO_CONTAINERS_H_
