/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_DRIVERS_CFG_READER_H__H_
#define ALQUIMIA_DRIVERS_CFG_READER_H_

/*******************************************************************************
 **
 ** C++ definition of the alquimia interface as an abstract base class
 **
 ******************************************************************************/

#include "alquimia_containers.h"

// Base class defining the alquimia interface
namespace alquimia {
namespace drivers {
namespace utilities {

struct SimulationParameters {
  std::string description;
  std::string engine;
  std::string engine_inputfile;
  std::string initial_condition;
  double delta_t;
  int num_time_steps;
  std::string use_text_output;
  std::string output_time_units;
};

class AlquimiaConfigReader {
 public:
  AlquimiaConfigReader() {};
  virtual ~AlquimiaConfigReader() {};

  void ReadInputFile(const std::string& file_name,
                     SimulationParameters* simulation_params,
                     alquimia::AlquimiaState* state,
                     alquimia::AlquimiaMaterialProperties* material_props,
                     alquimia::AlquimiaConditions* conditions);

  void PrintInput(const SimulationParameters& params,
                  const alquimia::AlquimiaState& state,
                  const alquimia::AlquimiaMaterialProperties& material_props,
                  const alquimia::AlquimiaConditions& conditions);
  
  void WriteTemplateFile(const std::string& file_name);

 protected:

 private:
  void GetLineCleaned(std::ifstream* input_file,
                      std::string* line);

  void ParseSimulationSection(std::ifstream* input_file,
                              SimulationParameters* params);
  
  void ParseStateSection(std::ifstream* input_file,
                         alquimia::AlquimiaState* state);
  void ParseMaterialPropertySection(
      std::ifstream* input_file,
      alquimia::AlquimiaMaterialProperties* material_props);
  
  void ParseConditionSection(
      std::ifstream* input_file,
      alquimia::AlquimiaConditions* geochemical_conditions);
  
  void PrintSimulationParameters(const SimulationParameters& params);
  void PrintStateParameters(const alquimia::AlquimiaState& state);
  void PrintMaterialPropertyParameters(
      const alquimia::AlquimiaMaterialProperties& material_props);
  void PrintGeochemicalConditions(
      const alquimia::AlquimiaConditions conditions);

  static const std::string kEqual;
  static const std::string kSpaces;

  static const std::string kSimulationSection;
  static const std::string kDescriptionString;
  static const std::string kEngineString;
  static const std::string kEngineInputfileString;
  static const std::string kICString;
  static const std::string kDeltaTimeString;
  static const std::string kNumTimeStepsString;
  static const std::string kUseTextOutputString;
  static const std::string kOutputTimeUnitsString;
  
  static const std::string kStateSection;
  static const std::string kDensityString;
  static const std::string kSaturationString;
  static const std::string kPorosityString;
  static const std::string kTemperatureString;
  static const std::string kPressureString;
  
  static const std::string kMaterialPropertiesSection;
  static const std::string kVolumeString;
  static const std::string kIsothermKdString;
  static const std::string kFreundlichNString;
  static const std::string kLangmuirBString;
  
  static const std::string kGeochemicalConditionsSection;
  static const std::string kNamedConditionSection;
  static const std::string kNameString;
  
};

}  // namespace utilities
}  // namespace drivers
}  // namespace alquimia
#endif  // ALQUIMIA_DRIVERS_CFG_READER_H_
