/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "cfg_reader.h"

#include <string>
#include <fstream>
#include <iostream>

#include "alquimia_containers.h"

#include "cc_demo_utils.h"
#include "string_tokenizer.h"
namespace alquimia {
namespace drivers {
namespace utilities {

const std::string AlquimiaConfigReader::kSimulationSection("simulation_parameters");
const std::string AlquimiaConfigReader::kDescriptionParam("description");
const std::string AlquimiaConfigReader::kEngineParam("engine");
const std::string AlquimiaConfigReader::kEngineInputfileParam("engine_inputfile");
const std::string AlquimiaConfigReader::kUseTextOutputParam("use_text_output");
const std::string AlquimiaConfigReader::kOutputTimeUnitsParam("output_time_units");

const std::string AlquimiaConfigReader::kStateSection("state");
const std::string AlquimiaConfigReader::kDensityParam("density");
const std::string AlquimiaConfigReader::kSaturationParam("saturation");
const std::string AlquimiaConfigReader::kPorosityParam("porosity");
const std::string AlquimiaConfigReader::kTemperatureParam("temperature");
const std::string AlquimiaConfigReader::kPressureParam("pressure");

const std::string AlquimiaConfigReader::kMaterialPropertiesSection("material_properties");
const std::string AlquimiaConfigReader::kVolumeParam("volume");
const std::string AlquimiaConfigReader::kIsothermKdParam("isotherm_kd");
const std::string AlquimiaConfigReader::kFreundlichNParam("freundlich_n");
const std::string AlquimiaConfigReader::kLangmuirBParam("langmuir_b");

const std::string AlquimiaConfigReader::kGeochemicalConditionsSection("geochemical_conditions");
const std::string AlquimiaConfigReader::kNamedConditionSection("condition");
const std::string AlquimiaConfigReader::kNameParam("name");
//const std::string AlquimiaConfigReader::kICParam("ic");

void AlquimiaConfigReader::ReadInputFile(
    const std::string& file_name,
    SimulationParameters* simulation_params,
    alquimia::AlquimiaState* state,
    alquimia::AlquimiaMaterialProperties* material_props,
    alquimia::AlquimiaConditions* conditions)
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

    // TODO(bja): need to add windows line end check to all getline() calls!
    // if ((raw_line.size() > 0) && (raw_line.at(raw_line.size() - 1) == '\r')) {
    //   // getline only searches for \n line ends. windows files use \r\n
    //   // check for a hanging \r and remove it if it is there
    //   raw_line.resize(raw_line.size() - 1);
    // }

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
      getline(input_file, raw_line);
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
                                                    kGeochemicalConditionsSection)) {
        ParseGeochemicalConditionsSection(&input_file, conditions);
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


void AlquimiaConfigReader::ParseSimulationSection(
    std::ifstream* input_file,
    SimulationParameters* simulation)
{
  namespace util = alquimia::drivers::utilities;
  std::string equal("=:");
  std::string spaces(" \t");
  std::string raw_line;
  bool new_section(false);
  while (!input_file->eof() && !new_section) {
    char next = input_file->peek();
    if (next == '[') {
      new_section = true;
    } else if (next == '#' || next == ' ') {
      // comment line, do nothing
      getline(*input_file, raw_line);
    } else {
      // parameter line
      getline(*input_file, raw_line);
      util::StringTokenizer params(raw_line, equal);
      //std::cout << "\'" << raw_line << "\'" << std::endl;
      // if params.size() == 0 then we have a blank line
      if (params.size() == 2) {
        std::string key = params.at(0);
        util::RemoveLeadingAndTrailingWhitespace(&key);
        std::string value = params.at(1);
        util::RemoveLeadingAndTrailingWhitespace(&value);
        //std::cout << "Parsing ----->  '" << params.at(0) << "'" << std::endl;
        if (util::CaseInsensitiveStringCompare(key, kDescriptionParam)) {
          // the description probably has spaces in it, so we want to use
          // the raw parameter value from params.at(1)
          simulation->description.assign(value);
        } else if (util::CaseInsensitiveStringCompare(key, kEngineParam)) {
          simulation->engine.assign(value);
        } else if (util::CaseInsensitiveStringCompare(key, kEngineInputfileParam)) {
          simulation->engine_inputfile.assign(value);
        } else if (util::CaseInsensitiveStringCompare(key, kUseTextOutputParam)) {
          simulation->use_text_output.assign(value);
        } else if (util::CaseInsensitiveStringCompare(key, kOutputTimeUnitsParam)) {
          simulation->output_time_units.assign(value);
        }
      }
    }  // end else(parameter line)
  }  // end while(input_file)
}  // end ParseSimulationSection()


void AlquimiaConfigReader::ParseStateSection(
    std::ifstream* input_file,
    alquimia::AlquimiaState* state)
{
  namespace util = alquimia::drivers::utilities;
  std::string equal("=:");
  std::string spaces(" \t");
  std::string raw_line;
  bool new_section(false);
  while (!input_file->eof() && !new_section) {
    char next = input_file->peek();
    if (next == '[') {
      new_section = true;
    } else if (next == '#' || next == ' ') {
      // comment line, do nothing
      getline(*input_file, raw_line);
    } else {
      // parameter line
      getline(*input_file, raw_line);
      util::StringTokenizer params(raw_line, equal);
      //std::cout << "\'" << raw_line << "\'" << std::endl;
      // if param.size() == 0 then we have a blank line
      if (params.size() == 2) {
        std::string key(params.at(0));
        util::RemoveLeadingAndTrailingWhitespace(&key);
        std::string value(params.at(1));
        util::RemoveLeadingAndTrailingWhitespace(&value);
        //std::cout << "Parsing ----->  '" << param.at(0) << "'" << std::endl;
        if (util::CaseInsensitiveStringCompare(key, kDensityParam)) {
          state->water_density = std::atof(value.c_str());
        } else if (util::CaseInsensitiveStringCompare(key, kSaturationParam)) {
          state->saturation = std::atof(value.c_str());
        } else if (util::CaseInsensitiveStringCompare(key, kPorosityParam)) {
          state->porosity = std::atof(value.c_str());
        } else if (util::CaseInsensitiveStringCompare(key, kTemperatureParam)) {
          state->temperature = std::atof(value.c_str());
        } else if (util::CaseInsensitiveStringCompare(key, kPressureParam)) {
          state->aqueous_pressure = std::atof(value.c_str());
        }
      }
    }  // end else(parameter line)
  }  // end while(input_file)
}  // end ParseStateSection();

void AlquimiaConfigReader::ParseMaterialPropertySection(
    std::ifstream* input_file,
    alquimia::AlquimiaMaterialProperties* material_props)
{
  namespace util = alquimia::drivers::utilities;
  std::string equal("=:");
  std::string spaces(" \t");
  std::string raw_line;
  bool new_section(false);
  while (!input_file->eof() && !new_section) {
    char next = input_file->peek();
    if (next == '[') {
      new_section = true;
    } else if (next == '#' || next == ' ') {
      // comment line, do nothing
      getline(*input_file, raw_line);
    } else {
      // parameter line
      getline(*input_file, raw_line);
      util::StringTokenizer params(raw_line, equal);
      //std::cout << "\'" << raw_line << "\'" << std::endl;
      // if params.size() == 0 then we have a blank line
      if (params.size() == 2) {
        std::string key = params.at(0);
        util::RemoveLeadingAndTrailingWhitespace(&key);
        std::string value = params.at(1);
        util::RemoveLeadingAndTrailingWhitespace(&value);
        // single valued parameters go first
        if (util::CaseInsensitiveStringCompare(key, kVolumeParam)) {
          material_props->volume = std::atof(value.c_str());
        } else {
          // now we deal with vector parameters.

          // TODO(bja): can we check or deal with the order/species names
          // somehow....
          util::StringTokenizer vec_values(value, ",");
          for (util::StringTokenizer::iterator s = vec_values.begin();
               s != vec_values.end(); ++s) {
            util::StringTokenizer data(*s, spaces);
            std::string species(data.at(0));
            util::RemoveLeadingAndTrailingWhitespace(&key);
            std::string value_str(data.at(1));
            util::RemoveLeadingAndTrailingWhitespace(&value_str);
            double species_value(std::atof(value_str.c_str()));
            //std::cout << "Parsing ----->  " << key << " : " << species << std::endl;
            if (util::CaseInsensitiveStringCompare(key, kIsothermKdParam)) {
              material_props->isotherm_kd.push_back(species_value);
            } else if (util::CaseInsensitiveStringCompare(key, kFreundlichNParam)) {
              material_props->freundlich_n.push_back(species_value);
            } else if (util::CaseInsensitiveStringCompare(key, kLangmuirBParam)) {
              material_props->langmuir_b.push_back(species_value);
            }
          }  // end for(s)
        }  // end else (vector parameters)
      }  // end if(params.size())
    }  // end else(parameter line)
  }  // end while(input_file)
}  // end ParseMaterialPropertiesSection();

void AlquimiaConfigReader::ParseGeochemicalConditionsSection(
    std::ifstream* input_file,
    alquimia::AlquimiaConditions* geochemical_conditions)
{
  namespace util = alquimia::drivers::utilities;
  std::string equal("=:");
  std::string spaces(" \t");
  std::string raw_line;
  bool new_section(false);
  while (!input_file->eof() && !new_section) {
    char next = input_file->peek();
    if (next == '[') {
      new_section = true;
    } else if (next == '#' || next == ' ') {
      // comment line, do nothing
      getline(*input_file, raw_line);
    } else {
      // parameter line
      getline(*input_file, raw_line);
      util::StringTokenizer params(raw_line, equal);
      //std::cout << "\'" << raw_line << "\'" << std::endl;
      // if params.size() == 0 then we have a blank line
      if (params.size() > 0) {
        std::string key = params.at(0);
        util::RemoveLeadingAndTrailingWhitespace(&key);
        // we ignore the values here
        //std::string value = params.at(1);
        //util::RemoveLeadingAndTrailingWhitespace(&value);
        (*geochemical_conditions)[key] = alquimia::AlquimiaGeochemicalCondition();
      }
    }  // end else(parameter line)
  }  // end while(input_file)
}  // end ParseGeochemicalConditionsSection();

void AlquimiaConfigReader::ParseConditionSection(
    std::ifstream* input_file,
    alquimia::AlquimiaConditions* geochemical_conditions)
{
  namespace util = alquimia::drivers::utilities;
  std::string equal("=:");
  std::string spaces(" \t");
  std::string condition_name("");
  std::string raw_line;
  bool new_section(false);
  while (!input_file->eof() && !new_section) {
    char next = input_file->peek();
    if (next == '[') {
      new_section = true;
    } else if (next == '#' || next == ' ') {
      // comment line, do nothing
      getline(*input_file, raw_line);
    } else {
      // parameter line
      getline(*input_file, raw_line);
      util::StringTokenizer params(raw_line, equal);
      //std::cout << "\'" << raw_line << "\'" << std::endl;
      // if params.size() == 0 then we have a blank line
      if (params.size() == 2) {
        std::string key = params.at(0);
        util::RemoveLeadingAndTrailingWhitespace(&key);
        std::string value = params.at(1);
        util::RemoveLeadingAndTrailingWhitespace(&value);
        // single valued parameters go first
        if (util::CaseInsensitiveStringCompare(key, kNameParam)) {
          condition_name = value;
        } else {
          // now we deal with vector parameters.

          // TODO(bja): can we check or deal with the order/species names
          // somehow....
          alquimia::AlquimiaGeochemicalConstraint constraint;
          constraint.primary_species = key;
          util::StringTokenizer constraint_data(value, spaces);
          constraint.value = constraint_data.at(0);
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


void AlquimiaConfigReader::PrintInput(
    const SimulationParameters& params,
    const alquimia::AlquimiaState& state,
    const alquimia::AlquimiaMaterialProperties& material_props,
    const alquimia::AlquimiaConditions& conditions)
{
  std::cout << "- Input File ---------------------------------------------------------\n";
  PrintSimulationParameters(params);
  PrintStateParameters(state);
  PrintMaterialPropertyParameters(material_props);
  PrintGeochemicalConditions(conditions);
  std::cout << "--------------------------------------------------------- Input File -\n";
}  // end PrintInput()


void AlquimiaConfigReader::PrintSimulationParameters(
    const SimulationParameters& sim_params)
{
  std::cout << "  -- Simulation parameters :" << std::endl;
  std::cout << "    description : " << sim_params.description << std::endl;
  std::cout << "    engine : " << sim_params.engine << std::endl;
  std::cout << "    engine inputfile : " << sim_params.engine_inputfile << std::endl;
  std::cout << "    text output : " << sim_params.use_text_output << std::endl;
  std::cout << "    output time units : " << sim_params.output_time_units << std::endl;
  std::cout << std::endl;
}  // end PrintSimulationParameters()

void AlquimiaConfigReader::PrintStateParameters(
    const alquimia::AlquimiaState& state)
{
  std::cout << "  -- State :" << std::endl;
  std::cout << "    density : " << state.water_density << std::endl;
  std::cout << "    saturation : " << state.saturation << std::endl;
  std::cout << "    porosity : " << state.porosity << std::endl;
  std::cout << "    temperature : " << state.temperature << std::endl;
  std::cout << "    pressure : " << state.aqueous_pressure << std::endl;
  std::cout << std::endl;
}  // end PrintStateParameters()

void AlquimiaConfigReader::PrintMaterialPropertyParameters(
    const alquimia::AlquimiaMaterialProperties& material_props)
{
  namespace util = alquimia::drivers::utilities;
  std::cout << "  -- Material Properties :" << std::endl;
  std::cout << "    volume : " << material_props.volume << std::endl;
  util::PrintVector("    isotherm_kd", material_props.isotherm_kd);
  util::PrintVector("    freundlich_n", material_props.freundlich_n);
  util::PrintVector("    langmuir_b", material_props.langmuir_b);
  std::cout << std::endl;
}  // end PrintMaterialPropertyParameters()

void AlquimiaConfigReader::PrintGeochemicalConditions(
    const alquimia::AlquimiaConditions conditions)
{
  std::cout << "  -- Geochemical Conditions :" << std::endl;
  for (alquimia::AlquimiaConditions::const_iterator c = conditions.begin();
       c != conditions.end(); ++c) {
    std::cout << "    " << c->first << " : " << std::endl;
    for (alquimia::AlquimiaGeochemicalCondition::const_iterator g = c->second.begin();
         g != c->second.end(); ++g) {
      std::cout << "        " << g->primary_species << std::endl;
      std::cout << "            type : " << g->constraint_type << std::endl;
      std::cout << "            associated : " << g->associated_species << std::endl;
      std::cout << "            value : " << g->value << std::endl;
    }
  }
}  // end PrintGeochemicalConditions()

void AlquimiaConfigReader::WriteTemplateFile(const std::string& file_name)
{
  std::ofstream template_file(file_name.c_str());
  if (!template_file) {
    std::cout << "ConfigReader : \n";
    std::cout << "template file \'" << file_name
              << "\' could not be opened." << std::endl;
    abort();
  }
  template_file << "[" << kSimulationSection << "]" << std::endl;
  template_file << kDescriptionParam << " = " << std::endl;
  template_file << kEngineParam << " = " << std::endl;
  template_file << kEngineInputfileParam << " = " << std::endl;
  template_file << std::endl;

  template_file << "[" << kStateSection << "]" << std::endl;
  template_file << kDensityParam << " = " << std::endl;
  template_file << kSaturationParam << " = " << std::endl;
  template_file << kPorosityParam << " = " << std::endl;
  template_file << kTemperatureParam << " = " << std::endl;
  template_file << kPressureParam << " = " << std::endl;
  template_file << std::endl;

  template_file << "[" << kMaterialPropertiesSection << "]" << std::endl;
  template_file << kVolumeParam << " = " << std::endl;
  template_file << kIsothermKdParam << " = [species_name value, ...]" << std::endl;
  template_file << kFreundlichNParam << " = [species_name value, ...]" << std::endl;
  template_file << kLangmuirBParam << " = [species_name value, ...]" << std::endl;
  template_file << std::endl;

  template_file << "[" << kGeochemicalConditionsSection << "]" << std::endl;
  template_file << "condition_name = " << std::endl;
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

