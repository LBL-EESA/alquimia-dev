/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_CXX_ALQUIMIA_INTERFACE_FACTORY_H_
#define ALQUIMIA_CXX_ALQUIMIA_INTERFACE_FACTORY_H_

#include "alquimia_interface.h"

#include <string>

namespace alquimia {

class AlquimiaInterfaceFactory {
 public:
  AlquimiaInterfaceFactory() {};
  ~AlquimiaInterfaceFactory() {};

  AlquimiaInterface* Create(const std::string& engine);

 protected:

 private:
};

}  // namespace alquimia
#endif  // ALQUIMIA_CXX_ALQUIMIA_INTERFACE_FACTORY_H_
