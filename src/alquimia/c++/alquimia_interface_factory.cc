/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "alquimia_interface_factory.h"

#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>

#include "pflotran_alquimia_interface.h"
#include "crunch_alquimia_interface.h"
#include "alquimia_interface.h"
#include "alquimia_strings.h"

namespace alquimia {

AlquimiaInterface* AlquimiaInterfaceFactory::Create(
    const std::string& engine_name) {

  AlquimiaInterface* interface = NULL;

  if (engine_name == alquimia::strings::kPFloTran) {
#ifdef HAVE_PFLOTRAN
    interface = new PFloTranAlquimiaInterface();
#else
    std::stringstream message;
    message << "\nERROR : AlquimiaInterfaceFactory::Create() : "
            << "PFloTran interface requested, but alquimia was not compiled with PFloTran!\n";
    throw std::runtime_error(message.str());
#endif

  } else if (engine_name == alquimia::strings::kCrunchFlow) {
#ifdef HAVE_CRUNCH
    interface = new CrunchFlowAlquimiaInterface();
#else
    std::stringstream message;
    message << "\nERROR : AlquimiaInterfaceFactory::Create() : "
            << "CrunchFlow interface requested, but alquimia was not compiled with CrunchFlow!\n";
    throw std::runtime_error(message.str());
#endif

  } else if  (engine_name == alquimia::strings::kCppChem) {
    std::stringstream message;
    message << "\nERROR : AlquimiaInterfaceFactory::Create() : "
            << "C++ chemistry interface has not been implemented!\n";
    throw std::runtime_error(message.str());

  } else {
    std::stringstream message;
    message << "\nERROR : AlquimiaInterfaceFactory::Create() : "
            << "The requested interface '" << engine_name << "' is invalid."
            << "\n  Valid engines names are:\n"
            << "    '" << alquimia::strings::kPFloTran << "'\n"
            << "    '" << alquimia::strings::kCrunchFlow << "'\n"
            << "    '" << alquimia::strings::kCppChem << "'\n";
    throw std::runtime_error(message.str());
  }

  if (interface != NULL) {
    std::cout << "AlquimiaInterfaceFactory::Create() : successfully created "
              << engine_name << " interface.\n";
  }

  return interface;
}  //  end Create()

}  //  namespace alquimia
