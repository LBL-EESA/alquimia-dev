/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_CXX_ALQUIMIA_CXX_CONTAINERS_H_
#define ALQUIMIA_CXX_ALQUIMIA_CXX_CONTAINERS_H_

/*******************************************************************************
 **
 ** C++ implementation of the alquimia containers
 **
 ******************************************************************************/

#include <vector>
#include <map>
#include <string>

namespace alquimia {

struct AlquimiaState {
  double water_density;  // [kg/m^3]
  double saturation;  // [-]
  double porosity;  // [-]
  double temperature;  // [celsius]
  double aqueous_pressure; // [Pa]
  std::vector<double> total_primary;  // [molarity]
  std::vector<double> total_sorbed;  // [moles/m^3 bulk]
  std::vector<double> free_ion;  // [molality]
  std::vector<double> mineral_volume_fraction;  // [-]
  std::vector<double> mineral_specific_surface_area; // [m^2 mineral/m^3 bulk]
  std::vector<double> cation_exchange_capacity;  // [moles/m^3 bulk]
  std::vector<double> surface_site_density;  // [moles/m^3 bulk]
};

struct AlquimiaMaterialProperties {
  double volume;  // [m^3]
  std::vector<double> isotherm_kd;  // [?]
  std::vector<double> freundlich_n; // [?]
  std::vector<double> langmuir_b;  // [?]
};

struct AlquimiaAuxiliaryData {
  std::vector<double> primary_activity_coeff;  // [-]
  std::vector<double> secondary_activity_coeff;  // [-]
  std::vector<double> ion_exchange_ref_cation_conc;  // [?]
  std::vector<double> surface_complex_free_site_conc;  // [?]
};

struct AlquimiaEngineStatus {
  unsigned int num_rhs_evaluations;
  unsigned int num_jacobian_evaluations;
  unsigned int num_newton_iterations;
  bool converged;
};

struct AlquimiaMetaData {
  bool thread_safe;
  bool temperature_dependent;
  bool pressure_dependent;
  bool porosity_update;
  std::vector<std::string> auxiliary_output_names;
};

struct AlquimiaGeochemicalConstraint {
  std::string primary_species;
  std::string constraint_type;
  std::string associated_species;
  std::string value;
};
  
typedef std::vector<AlquimiaGeochemicalConstraint> AlquimiaGeochemicalCondition;

typedef std::map<std::string, AlquimiaGeochemicalCondition> AlquimiaConditions;

struct OutputData {
  double pH;
  std::vector<double> mineral_saturation_index;  // [-]
  std::vector<double> mineral_reaction_rate;  // [?]
};

}  // namespace alquimia
#endif  // ALQUIMIA_CXX_ALQUIMIA_CXX_CONTAINERS_H_
