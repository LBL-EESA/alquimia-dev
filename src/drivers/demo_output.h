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

#ifndef ALQUIMIA_DRIVERS_DEMO_OUTPUT_H_
#define ALQUIMIA_DRIVERS_DEMO_OUTPUT_H_

/*******************************************************************************
 **
 ** Output classes so the batch chem driver can write reaction
 ** stepping output in a way that is easily comparable to the engines
 ** observation file.
 **
 ******************************************************************************/

#include <iostream>
#include <fstream>
#include <string>

#include "alquimia_containers.h"

#include "demo_utils.h"

namespace alquimia {
namespace drivers {
namespace utilities {

class DemoOutput {
 public:
  DemoOutput(void);
  virtual ~DemoOutput();

  virtual void Setup(const std::string& filename);

  virtual void WriteHeader(const char time_units,
                           const AlquimiaProblemMetaData& meta_data,
                           const AlquimiaSizes& sizes) = 0;

  virtual void Write(const double time,
                     const AlquimiaState& state,
                     const AlquimiaAuxiliaryOutputData& aux_output) = 0;

 protected:
  std::fstream output;

};

class DemoOutputPFloTran : public DemoOutput {
 public:
  DemoOutputPFloTran(void);
  virtual ~DemoOutputPFloTran();

  virtual void WriteHeader(const char time_units,
                           const AlquimiaProblemMetaData& meta_data,
                           const AlquimiaSizes& sizes);

  virtual void Write(const double time,
                     const AlquimiaState& state,
                     const AlquimiaAuxiliaryOutputData& aux_output);

 protected:

 private:
  bool write_pH;
};



}  // namespace utilities
}  // namespace drivers
}  // namespace alquimia
#endif  // ALQUIMIA_DRIVERS_DEMO_OUTPUT_H_
