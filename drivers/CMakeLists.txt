include_directories(${ALQUIMIA_INCLUDE_DIRS})

# Batch chemistry (0D) driver.
add_executable(batch_chem BatchChemDriver.c DriverOutput.c batch_chem.c ini.c input_util.c)
target_link_libraries(batch_chem alquimia;${ALQUIMIA_TPLS})
add_dependencies(batch_chem alquimia)

# Reactive transport (1D) driver.
add_executable(transport TransportDriver.c DriverOutput.c transport.c ini.c)
target_link_libraries(transport alquimia;${ALQUIMIA_TPLS})
add_dependencies(batch_chem alquimia)

# Install executables.
install(TARGETS batch_chem transport DESTINATION bin)

