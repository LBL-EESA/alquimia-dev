/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

//
// Alquimia Copyright (c) 2013, The Regents of the University of California, 
// through Lawrence Berkeley National Laboratory (subject to receipt of any 
// required approvals from the U.S. Dept. of Energy).  All rights reserved.
// 
// If you have questions about your rights to use or distribute this software, 
// please contact Berkeley Lab's Technology Transfer and Intellectual Property 
// Management at TTD@lbl.gov referring to Alquimia (LBNL Ref. 2013-119).
// 
// NOTICE.  This software was developed under funding from the U.S. Department 
// of Energy.  As such, the U.S. Government has been granted for itself and 
// others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide 
// license in the Software to reproduce, prepare derivative works, and perform 
// publicly and display publicly.  Beginning five (5) years after the date 
// permission to assert copyright is obtained from the U.S. Department of Energy, 
// and subject to any subsequent five (5) year renewals, the U.S. Government is 
// granted for itself and others acting on its behalf a paid-up, nonexclusive, 
// irrevocable, worldwide license in the Software to reproduce, prepare derivative
// works, distribute copies to the public, perform publicly and display publicly, 
// and to permit others to do so.
// 
// Authors: Benjamin Andre <bandre@lbl.gov>
//

#include "demo_output.h"

#include <iostream>
#include <fstream>

#include "alquimia_containers.h"
#include "alquimia_util.h"

namespace alquimia {
namespace drivers {
namespace utilities {

/*******************************************************************************
 **
 **  DemoOutput base class
 **
 *******************************************************************************/
DemoOutput::DemoOutput(void) {
}  // end DemoOutput()

DemoOutput::~DemoOutput(void) {
  this->output.close();
}  // end ~DemoOutput()

void DemoOutput::Setup(const std::string& filename) {
  size_t position = filename.find_last_of('.');
  std::string text_output_name = filename.substr(0, position) + ".txt";
  this->output.open(text_output_name.c_str(), std::fstream::out);
}  // end DemoOutput()

/*******************************************************************************
 **
 **  PFloTran
 **
 *******************************************************************************/
DemoOutputPFloTran::DemoOutputPFloTran(void)
    : DemoOutput() {
}  // end DemoOutputPFloTran()

DemoOutputPFloTran::~DemoOutputPFloTran(void) {
}  // end ~DemoOutputPFloTran()

void DemoOutputPFloTran::WriteHeader(const char time_units,
                                     const AlquimiaProblemMetaData& meta_data,
                                     const AlquimiaSizes& sizes) {
  if (this->output.is_open()) {
    this->output << "# \"Time [" << time_units << "]\"";
    int h_index;
    AlquimiaFindIndexFromName("H+", &meta_data.primary_names, &h_index);
    if (h_index >= 0) {
      this->write_pH = true;
      this->output << " , \"pH\"";
    }
    for (int i = 0; i < meta_data.primary_names.size; ++i) {
      this->output <<  " , \"Total " << meta_data.primary_names.data[i] << " [M]\"";
    }

    for (int i = 0; i < sizes.num_sorbed; ++i) {
      this->output << " , \"Total Sorbed " << meta_data.primary_names.data[i]
                   << " [mol/m^3]\"";
    }

    for (int i = 0; i < meta_data.mineral_names.size; ++i) {
      this->output << " , \"" << meta_data.mineral_names.data[i] << " VF\"";
    }
    for (int i = 0; i < meta_data.mineral_names.size; ++i) {
      this->output << " , \"" << meta_data.mineral_names.data[i] << " Rate [mol/m^3/sec]\"";
    }
    this->output << std::endl;
  }
}  // end DemoOutputPFloTran::WriteHeader()

void DemoOutputPFloTran::Write(const double time,
                               const AlquimiaState& state,
                               const AlquimiaAuxiliaryOutputData& aux_output) {
  if (this->output.is_open()) {
    std::string seperator("  ");
    this->output << std::scientific << std::uppercase
                 << std::setprecision(6);
    this->output << seperator << time;
    if (this->write_pH) {
      this->output << seperator << aux_output.pH;
    }
    for (int i = 0; i < state.total_mobile.size; ++i) {
      this->output << seperator << state.total_mobile.data[i];
    }
    for (int i = 0; i < state.total_immobile.size; ++i) {
      this->output << seperator << state.total_immobile.data[i];
    }
    for (int i = 0; i < state.mineral_volume_fraction.size; ++i) {
      this->output << seperator << state.mineral_volume_fraction.data[i];
    }
    for (int i = 0; i < aux_output.mineral_reaction_rate.size; ++i) {
      this->output << seperator << aux_output.mineral_reaction_rate.data[i];
    }
    this->output << std::endl;
  }
}  // end DemoOutputPFloTran::Write()


/*******************************************************************************
 **
 **  CrunchFLow
 **
 *******************************************************************************/


}  // namespace utilities
}  // namespace drivers
}  // namespace alquimia
