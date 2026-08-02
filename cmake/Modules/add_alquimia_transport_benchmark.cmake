# This function adds a reactive transport benchmark test for Alquimia.
function(add_alquimia_transport_benchmark benchmark input_file)
  set(exe ${PROJECT_BINARY_DIR}/drivers/transport)
  set(test_name transport_${benchmark}__smoke)
  add_test(${test_name} ${exe} ${input_file})
  set_tests_properties(${test_name} PROPERTIES
    LABELS "smoke;transport"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endfunction()
