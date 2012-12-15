/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "cfg_reader.h"

#include <cstdlib>

#include <string>
#include <fstream>
#include <iostream>

#include "alquimia_containers.h"

#include "demo_utils.h"
#include "string_tokenizer.h"
namespace alquimia {
namespace drivers {
namespace utilities {

const  std::string DemoConfigReader::kEqual("=:");
const  std::string DemoConfigReader::kSpaces(" \t");

const std::string DemoConfigReader::kSimulationSection("simulation_parameters");
const std::string DemoConfigReader::kDescriptionString("description");
const std::string DemoConfigReader::kEngineString("engine");
const std::string DemoConfigReader::kEngineInputfileString("engine_inputfile");
const std::string DemoConfigReader::kICString("initial_condition");
const std::string DemoConfigReader::kDeltaTimeString("delta_t");
const std::string DemoConfigReader::kNumTimeStepsString("num_time_steps");
const std::string DemoConfigReader::kUseTextOutputString("use_text_output");
const std::string DemoConfigReader::kOutputTimeUnitsString("output_time_units");

const std::string DemoConfigReader::kStateSection("state");
const std::string DemoConfigReader::kDensityString("density");
const std::string DemoConfigReader::kSaturationString("saturation");
const std::string DemoConfigReader::kPorosityString("porosity");
const std::string DemoConfigReader::kTemperatureString("temperature");
const std::string DemoConfigReader::kPressureString("pressure");

const std::string DemoConfigReader::kMaterialPropertiesSection("material_properties");
const std::string DemoConfigReader::kVolumeString("volume");
const std::string DemoConfigReader::kIsothermKdString("isotherm_kd");
const std::string DemoConfigReader::kFreundlichNString("freundlich_n");
const std::string DemoConfigReader::kLangmuirBString("langmuir_b");

const std::string DemoConfigReader::kGeochemicalConditionsSection("geochemical_conditions");
const std::string DemoConfigReader::kNamedConditionSection("condition");
const std::string DemoConfigReader::kNameString("name");

void DemoConfigReader::ReadInputFile(
    const std::string& file_name,
    DemoSimulation* simulation_params,
    DemoState* state,
    DemoMaterialProperties* material_props,
    DemoConditions* conditions)
{
  namespace util = alquimia::drivers::utilities;
  std::ifstream input_file(file_name.c_str());
  if (!input_file) {
    std::cout << "ConfigReader() : \n";
    std::cout << "input file \'" << file_name
              << "\' could not be opened." << std::endl;
    abort();
  }

  enum LineType {
    kCommentLine,
    kSection,
    kParameter
  } line_type;

  int count = 0;
  int max_count = 50;
  while (!input_file.eof() && count < max_count) {
    ++count;
    char next = input_file.peek();
    //std::cout << "ReadInputFile() : next = " << next << std::endl;

    // char first = '\0';
    // if (raw_line.length() > 0) {
    //   first = raw_line[0];
    // }

    if (next == '#' || next == ' ' || next == '\t') {
      // NOTE: treating first character = whitespace as a comment!
      line_type = kCommentLine;
    } else if (next == '\n' || next == '\0') {
      // treat blank lines as comments.
      line_type = kCommentLine;
    } else if (next == '[') {
      line_type = kSection;
    } else {
      line_type = kParameter;
    }

    if (line_type == kSection) {
      std::string raw_line;
      GetLineCleaned(&input_file, &raw_line);
      //std::cout << raw_line << std::endl;
      size_t first = raw_line.find_first_not_of('[');
      size_t last = raw_line.find_last_of(']');
      --last;
      std::string section_name = raw_line.substr(first, last);
      if (util::CaseInsensitiveStringCompare(section_name,
                                             kSimulationSection)) {
        ParseSimulationSection(&input_file, simulation_params);
      } else if (util::CaseInsensitiveStringCompare(section_name,
                                                    kStateSection)) {
        ParseStateSection(&input_file, state);
      } else if (util::CaseInsensitiveStringCompare(section_name,
                                                    kMaterialPropertiesSection)) {
        ParseMaterialPropertySection(&input_file, material_props);
      } else if (util::CaseInsensitiveStringCompare(section_name,
                                                    kNamedConditionSection)) {
        ParseConditionSection(&input_file, conditions);
      } else {
        std::cout << "Unknown section name '" << section_name << "'n";
      }
    }
  }  // end while(input_file)

  input_file.close();
}  // end ReadInputFile()

void DemoConfigReader::GetLineCleaned(
    std::ifstream* input_file,
    std::string* line) {
  namespace util = alquimia::drivers::utilities;
  getline(*input_file, *line);
  if ((line->size() > 0) && (line->at(line->size() - 1) == '\r')) {
    // getline only searches for \n line ends. windows files use \r\n
    // check for a hanging \r and remove it if it is there
    line->resize(line->size() - 1);
  }
  util::RemoveLeadingAndTrailingWhitespace(line);
}  // end GetLineCleaned()


void DemoConfigReader::ParseSimulationSection(
    std::ifstream* input_file,
    DemoSimulation* simulation)
{
  namespace util = alquimia::drivers::utilities;
  std::string raw_line;
  bool new_section(false);
  while (!input_file->eof() && !new_section) {
    char next = input_file->peek();
    if (next == '[') {
      new_section = true;
    } else if (next == '#' || next == ' ') {
      // comment line, do nothing
      GetLineCleaned(input_file, &raw_line);
    } else {
      // parameter line
      GetLineCleaned(input_file, &raw_line);
      util::StringTokenizer params(raw_line, kEqual);
      //std::cout << "\'" << raw_line << "\'" << std::endl;
      // if params.size() == 0 then we have a blank line
      if (params.size() == 2) {
        std::string key = params.at(0);
        util::RemoveLeadingAndTrailingWhitespace(&key);
        std::string value = params.at(1);
        util::RemoveLeadingAndTrailingWhitespace(&value);
        //std::cout << "Parsing ----->  '" << params.at(0) << "'" << std::endl;
        if (util::CaseInsensitiveStringCompare(key, kDescriptionString)) {
          // the description probably has spaces in it, so we want to use
          // the raw parameter value from params.at(1)
          simulation->description.assign(value);
        } else if (util::CaseInsensitiveStringCompare(key, kEngineString)) {
          simulation->engine.assign(value);
        } else if (util::CaseInsensitiveStringCompare(key, kEngineInputfileString)) {
          simulation->engine_inputfile.assign(value);
        } else if (util::CaseInsensitiveStringCompare(key, kDeltaTimeString)) {
          simulation->delta_t = std::atof(value.c_str());
        } else if (util::CaseInsensitiveStringCompare(key, kNumTimeStepsString)) {
          simulation->num_time_steps = std::atoi(value.c_str());
        } else if (util::CaseInsensitiveStringCompare(key, kUseTextOutputString)) {
          simulation->use_text_output.assign(value);
        } else if (util::CaseInsensitiveStringCompare(key, kOutputTimeUnitsString)) {
          simulation->output_time_units.assign(value);
        } else if (util::CaseInsensitiveStringCompare(key, kICString)) {
          simulation->initial_condition.assign(value);
        }
      }
    }  // end else(parameter line)
  }  // end while(input_file)
}  // end ParseSimulationSection()


void DemoConfigReader::ParseStateSection(
    std::ifstream* input_file,
    DemoState* state)
{
  namespace util = alquimia::drivers::utilities;
  std::string raw_line;
  bool new_section(false);
  while (!input_file->eof() && !new_section) {
    char next = input_file->peek();
    if (next == '[') {
      new_section = true;
    } else if (next == '#' || next == ' ') {
      // comment line, do nothing
      GetLineCleaned(input_file, &raw_line);
    } else {
      // parameter line
      GetLineCleaned(input_file, &raw_line);
      util::StringTokenizer params(raw_line, kEqual);
      //std::cout << "\'" << raw_line << "\'" << std::endl;
      // if param.size() == 0 then we have a blank line
      if (params.size() == 2) {
        std::string key(params.at(0));
        util::RemoveLeadingAndTrailingWhitespace(&key);
        std::string value(params.at(1));
        util::RemoveLeadingAndTrailingWhitespace(&value);
        //std::cout << "Parsing ----->  '" << param.at(0) << "'" << std::endl;
        if (util::CaseInsensitiveStringCompare(key, kDensityString)) {
          state->water_density = std::atof(value.c_str());
        } else if (util::CaseInsensitiveStringCompare(key, kSaturationString)) {
          state->saturation = std::atof(value.c_str());
        } else if (util::CaseInsensitiveStringCompare(key, kPorosityString)) {
          state->porosity = std::atof(value.c_str());
        } else if (util::CaseInsensitiveStringCompare(key, kTemperatureString)) {
          state->temperature = std::atof(value.c_str());
        } else if (util::CaseInsensitiveStringCompare(key, kPressureString)) {
          state->aqueous_pressure = std::atof(value.c_str());
        }
      }
    }  // end else(parameter line)
  }  // end while(input_file)
}  // end ParseStateSection();

void DemoConfigReader::ParseMaterialPropertySection(
    std::ifstream* input_file,
    DemoMaterialProperties* material_props)
{
  namespace util = alquimia::drivers::utilities;
  std::string raw_line;
  bool new_section(false);
  while (!input_file->eof() && !new_section) {
    char next = input_file->peek();
    if (next == '[') {
      new_section = true;
    } else if (next == '#' || next == ' ') {
      // comment line, do nothing
      GetLineCleaned(input_file, &raw_line);
    } else {
      // parameter line
      GetLineCleaned(input_file, &raw_line);
      util::StringTokenizer params(raw_line, kEqual);
      //std::cout << "\'" << raw_line << "\'" << std::endl;
      // if params.size() == 0 then we have a blank line
      if (params.size() == 2) {
        std::string key = params.at(0);
        util::RemoveLeadingAndTrailingWhitespace(&key);
        std::string value = params.at(1);
        util::RemoveLeadingAndTrailingWhitespace(&value);
        // single valued parameters go first
        if (util::CaseInsensitiveStringCompare(key, kVolumeString)) {
          material_props->volume = std::atof(value.c_str());
        } else {
          // now we deal with vector parameters.

          // TODO(bja): can we check or deal with the order/species names
          // somehow....
          util::StringTokenizer vec_values(value, ",");
          for (util::StringTokenizer::iterator s = vec_values.begin();
               s != vec_values.end(); ++s) {
            util::StringTokenizer data(*s, kSpaces);
            std::string species(data.at(0));
            util::RemoveLeadingAndTrailingWhitespace(&key);
            std::string value_str(data.at(1));
            util::RemoveLeadingAndTrailingWhitespace(&value_str);
            double species_value(std::atof(value_str.c_str()));
            //std::cout << "Parsing ----->  " << key << " : " << species << std::endl;
            if (util::CaseInsensitiveStringCompare(key, kIsothermKdString)) {
              material_props->isotherm_kd.push_back(species_value);
            } else if (util::CaseInsensitiveStringCompare(key, kFreundlichNString)) {
              material_props->freundlich_n.push_back(species_value);
            } else if (util::CaseInsensitiveStringCompare(key, kLangmuirBString)) {
              material_props->langmuir_b.push_back(species_value);
            }
          }  // end for(s)
        }  // end else (vector parameters)
      }  // end if(params.size())
    }  // end else(parameter line)
  }  // end while(input_file)
}  // end ParseMaterialPropertiesSection();

void DemoConfigReader::ParseConditionSection(
    std::ifstream* input_file,
    DemoConditions* geochemical_conditions)
{
  namespace util = alquimia::drivers::utilities;
  std::string condition_name("");
  std::string raw_line;
  bool new_section(false);
  while (!input_file->eof() && !new_section) {
    char next = input_file->peek();
    if (next == '[') {
      new_section = true;
    } else if (next == '#' || next == ' ') {
      // comment line, do nothing
      GetLineCleaned(input_file, &raw_line);
    } else {
      // parameter line
      GetLineCleaned(input_file, &raw_line);
      util::StringTokenizer params(raw_line, kEqual);
      //std::cout << "\'" << raw_line << "\'" << std::endl;
      // if params.size() == 0 then we have a blank line
      if (params.size() == 2) {
        std::string key = params.at(0);
        util::RemoveLeadingAndTrailingWhitespace(&key);
        std::string value = params.at(1);
        util::RemoveLeadingAndTrailingWhitespace(&value);
        // single valued parameters go first
        if (util::CaseInsensitiveStringCompare(key, kNameString)) {
          condition_name = value;
          (*geochemical_conditions)[condition_name] = DemoGeochemicalCondition();
        } else {
          // now we deal with vector parameters.

          // TODO(bja): can we check or deal with the order/species names
          // somehow....
          DemoGeochemicalConstraint constraint;
          constraint.primary_species = key;
          util::StringTokenizer constraint_data(value, kSpaces);
          constraint.value = std::atof(constraint_data.at(0).c_str());
          constraint.constraint_type = constraint_data.at(1);
          if (constraint_data.size() > 2) {
            constraint.associated_species = constraint_data.at(2);
          }
          (*geochemical_conditions)[condition_name].push_back(constraint);
        }  // end else (vector parameters)
      }  // end if(params.size())
    }  // end else(parameter line)
  }  // end while(input_file)
}  // end ParseConditionSection();


void DemoConfigReader::WriteTemplateFile(const std::string& file_name)
{
  std::ofstream template_file(file_name.c_str());
  if (!template_file) {
    std::cout << "ConfigReader : \n";
    std::cout << "template file \'" << file_name
              << "\' could not be opened." << std::endl;
    abort();
  }
  template_file << "[" << kSimulationSection << "]" << std::endl;
  template_file << kDescriptionString << " = " << std::endl;
  template_file << kEngineString << " = " << std::endl;
  template_file << kEngineInputfileString << " = " << std::endl;
  template_file << kICString << " = " << std::endl;
  template_file << kDeltaTimeString << " = " << std::endl;
  template_file << kNumTimeStepsString << " = " << std::endl;
  template_file << std::endl;

  template_file << "[" << kStateSection << "]" << std::endl;
  template_file << kDensityString << " = " << std::endl;
  template_file << kSaturationString << " = " << std::endl;
  template_file << kPorosityString << " = " << std::endl;
  template_file << kTemperatureString << " = " << std::endl;
  template_file << kPressureString << " = " << std::endl;
  template_file << std::endl;

  template_file << "[" << kMaterialPropertiesSection << "]" << std::endl;
  template_file << kVolumeString << " = " << std::endl;
  template_file << kIsothermKdString << " = [species_name value, ...]" << std::endl;
  template_file << kFreundlichNString << " = [species_name value, ...]" << std::endl;
  template_file << kLangmuirBString << " = [species_name value, ...]" << std::endl;
  template_file << std::endl;

  template_file << "[condition]" << std::endl;
  template_file << "name = condition_name" << std::endl;
  template_file << "species_name = value type association"<< std::endl;
  template_file << std::endl;
  template_file.close();
}  // end WriteTemplateFile()



}  //  namespace utilities
}  //  namespace drivers
}  //  namespace alquimia

