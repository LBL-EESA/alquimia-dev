include_directories(${ALQUIMIA_INCLUDE_DIRS})
link_libraries(alquimia;${ALQUIMIA_TPLS})

# Batch chemistry (0D) driver.
add_executable(batch_chem batch_chem.cc cfg_reader.cc demo_containers.cc 
                          demo_output.cc demo_utils.cc string_tokenizer.cc
                          fpe.c)

# Reactive transport (1D) driver.
add_executable(transport TransportDriver.c DriverOutput.c transport.c ini.c)

# Install executable.
install(TARGETS batch_chem transport DESTINATION bin)

