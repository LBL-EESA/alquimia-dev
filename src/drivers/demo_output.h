/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
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
