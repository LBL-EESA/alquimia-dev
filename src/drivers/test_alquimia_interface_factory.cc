/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

#include "alquimia_interface_factory.h"
#include "alquimia_interface.h"
#include "alquimia_containers.h"
#include "alquimia_strings.h"

int main(int argc, char* argv[])
{
  alquimia::AlquimiaInterfaceFactory aif;

  try {
    std::string engine_name(alquimia::strings::kAlquimiaStringPFloTran);
    alquimia::AlquimiaInterface* interface = aif.Create(engine_name);
    std::cout << "Successfully created '" << engine_name << "' interface!" << std::endl;
    delete interface;
  } catch (std::exception& e) {
    std::cout << e.what() << std::endl;

  }

  try {
    std::string engine_name(alquimia::strings::kAlquimiaStringCrunchFlow);
    alquimia::AlquimiaInterface* interface = aif.Create(engine_name);
    std::cout << "Successfully created '" << engine_name << "' interface!" << std::endl;
    delete interface;
  } catch (std::exception& e) {
    std::cout << e.what() << std::endl;

  }

  try {
    std::string engine_name("foo");
    alquimia::AlquimiaInterface* interface = aif.Create(engine_name);
    std::cout << "Successfully created '" << engine_name << "' interface!" << std::endl;
    delete interface;
  } catch (std::exception& e) {
    std::cout << e.what() << std::endl;

  }
  
  
  return EXIT_SUCCESS;
}
