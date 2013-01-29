/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_DRIVERS_CFG_READER_H_
#define ALQUIMIA_DRIVERS_CFG_READER_H_

/*******************************************************************************
 **
 ** Simple config style input file reader
 **
 ******************************************************************************/

#include "demo_containers.h"

namespace alquimia {
namespace drivers {
namespace utilities {

class DemoConfigReader {
 public:
  DemoConfigReader() {};
  virtual ~DemoConfigReader() {};

  void ReadInputFile(const std::string& file_name,
                     DemoSimulation* simulation_params,
                     DemoState* state,
                     DemoMaterialProperties* material_props,
                     DemoConditions* conditions);

  void WriteTemplateFile(const std::string& file_name);

  void set_debug(const bool debug) {
    debug_ = debug;
  }

  bool debug(void) const {
    return debug_;
  }

 protected:

 private:
  void GetLineCleaned(std::ifstream* input_file,
                      std::string* line);

  void ParseSimulationSection(std::ifstream* input_file,
                              DemoSimulation* params);
  
  void ParseStateSection(std::ifstream* input_file,
                         DemoState* state);
  void ParseMaterialPropertySection(
      std::ifstream* input_file,
      DemoMaterialProperties* material_props);
  
  void ParseConditionSection(
      std::ifstream* input_file,
      DemoConditions* geochemical_conditions);
  
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
  static const std::string kAqueousString;
  static const std::string kMineralString;

  bool debug_;
  
};

}  // namespace utilities
}  // namespace drivers
}  // namespace alquimia
#endif  // ALQUIMIA_DRIVERS_CFG_READER_H_
