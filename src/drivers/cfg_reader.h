/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
** Alquimia Copyright (c) 2013, The Regents of the University of California, 
** through Lawrence Berkeley National Laboratory (subject to receipt of any 
** required approvals from the U.S. Dept. of Energy).  All rights reserved.
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
  static const std::string kComma;

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
  static const std::string kCECString;
  static const std::string kSiteDensityString;
  
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
