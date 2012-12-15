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
  std::string use_text_output;
  std::string output_time_units;

  void Print(void) const {
    std::cout << "  -- Simulation parameters :" << std::endl;
    std::cout << "    description : " << this->description << std::endl;
    std::cout << "    engine : " << this->engine << std::endl;
    std::cout << "    engine inputfile : " << this->engine_inputfile << std::endl;
    std::cout << "    delta t : " << this->delta_t << std::endl;
    std::cout << "    num times steps : " << this->num_time_steps << std::endl;
    std::cout << "    initial condition : " << this->initial_condition << std::endl;
    std::cout << "    text output : " << this->use_text_output << std::endl;
    std::cout << "    output time units : " << this->output_time_units << std::endl;
    std::cout << std::endl;

  };
};


struct DemoState {
  double water_density;  // [kg/m^3]
  double saturation;  // [-]
  double porosity;  // [-]
  double temperature;  // [celsius]
  double aqueous_pressure; // [Pa]

  void Print(void) const {
    std::cout << "  -- State :" << std::endl;
    std::cout << "    density : " << this->water_density << std::endl;
    std::cout << "    saturation : " << this->saturation << std::endl;
    std::cout << "    porosity : " << this->porosity << std::endl;
    std::cout << "    temperature : " << this->temperature << std::endl;
    std::cout << "    pressure : " << this->aqueous_pressure << std::endl;
    std::cout << std::endl;
  };
};

struct DemoMaterialProperties {
  double volume;  // [m^3]
  std::vector<double> isotherm_kd;  // [?]
  std::vector<double> freundlich_n; // [?]
  std::vector<double> langmuir_b;  // [?]
  void Print(void) const {
    namespace util = alquimia::drivers::utilities;
    std::cout << "  -- Material Properties :" << std::endl;
    std::cout << "    volume : " << this->volume << std::endl;
    util::PrintVector("    isotherm_kd", this->isotherm_kd);
    util::PrintVector("    freundlich_n", this->freundlich_n);
    util::PrintVector("    langmuir_b", this->langmuir_b);
    std::cout << std::endl;
  };
};

struct DemoGeochemicalConstraint {
  std::string primary_species;
  std::string constraint_type;
  std::string associated_species;
  double value;
  void Print(void) const {
    std::cout << "        " << this->primary_species << std::endl;
    std::cout << "            type : " << this->constraint_type << std::endl;
    std::cout << "            associated : " << this->associated_species << std::endl;
    std::cout << "            value : " << this->value << std::endl;
  };
};
  
typedef std::vector<DemoGeochemicalConstraint> DemoGeochemicalCondition;

typedef std::map<std::string, DemoGeochemicalCondition> DemoConditions;


}  // namespace utilities
}  // namespace drivers
}  // namespace alquimia
#endif  // ALQUIMIA_DRIVERS_DEMO_CONTAINERS_H_
