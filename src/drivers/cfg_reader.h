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

  void ParseSimulationSection(std::ifstream* input_file,
                              SimulationParameters* params);
  
  void ParseStateSection(std::ifstream* input_file,
                         alquimia::AlquimiaState* state);
  void ParseMaterialPropertySection(
      std::ifstream* input_file,
      alquimia::AlquimiaMaterialProperties* material_props);
  
  void ParseGeochemicalConditionsSection(
      std::ifstream* input_file,
      alquimia::AlquimiaConditions* geochemical_conditions);
  
  void ParseConditionSection(
      std::ifstream* input_file,
      alquimia::AlquimiaConditions* geochemical_conditions);
  
  void PrintSimulationParameters(const SimulationParameters& params);
  void PrintStateParameters(const alquimia::AlquimiaState& state);
  void PrintMaterialPropertyParameters(
      const alquimia::AlquimiaMaterialProperties& material_props);
  void PrintGeochemicalConditions(
      const alquimia::AlquimiaConditions conditions);

  static const std::string kSimulationSection;
  static const std::string kDescriptionParam;
  static const std::string kEngineParam;
  static const std::string kEngineInputfileParam;
  static const std::string kUseTextOutputParam;
  static const std::string kOutputTimeUnitsParam;
  
  static const std::string kStateSection;
  static const std::string kDensityParam;
  static const std::string kSaturationParam;
  static const std::string kPorosityParam;
  static const std::string kTemperatureParam;
  static const std::string kPressureParam;
  
  static const std::string kMaterialPropertiesSection;
  static const std::string kVolumeParam;
  static const std::string kIsothermKdParam;
  static const std::string kFreundlichNParam;
  static const std::string kLangmuirBParam;
  
  static const std::string kGeochemicalConditionsSection;
  static const std::string kNamedConditionSection;
  static const std::string kNameParam;
  //static const std::string kICParam;
  
};

}  // namespace utilities
}  // namespace drivers
}  // namespace alquimia
#endif  // ALQUIMIA_DRIVERS_CFG_READER_H_
