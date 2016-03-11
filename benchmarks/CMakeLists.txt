# NOTE: Currently all the benchmarks are just smoke tests.

# Batch chemistry benchmarks.
add_subdirectory(batch_chem)

# Transport benchmarks.
add_subdirectory(transport)

# Sanity checks.
if (ALQUIMIA_HAVE_PFLOTRAN)
  add_custom_target(pflotran_batch_chem_sanity_check ${CMAKE_INSTALL_PREFIX}/bin/batch_chem calcite-short-pc-pflotran.cfg
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/batch_chem)
  add_custom_target(pflotran_transport_sanity_check ${CMAKE_INSTALL_PREFIX}/bin/transport tracer-1d-pflotran.cfg
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/transport)
  set(sanity_checks pflotran_batch_chem_sanity_check pflotran_transport_sanity_check)
endif()
if (ALQUIMIA_HAVE_CRUNCHFLOW)
  add_custom_target(crunchflow_batch_chem_sanity_check ${CMAKE_INSTALL_PREFIX}/bin/batch_chem calcite-short-cc-crunch.cfg
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/batch_chem)
  add_custom_target(crunchflow_transport_sanity_check ${CMAKE_INSTALL_PREFIX}/bin/transport tracer-1d-crunch.cfg
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/transport)
  set(sanity_checks ${sanity_checks} crunchflow_batch_chem_sanity_check crunchflow_transport_sanity_check)
endif()
add_custom_target(test_install ${CMAKE_COMMAND} -E echo "All sanity checks PASSED." DEPENDS ${sanity_checks})
