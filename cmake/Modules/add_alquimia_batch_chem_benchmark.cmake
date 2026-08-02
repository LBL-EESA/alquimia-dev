# This function adds a batch chemistry benchmark test for Alquimia.
function(add_alquimia_batch_chem_benchmark benchmark input_file)
  set(exe ${PROJECT_BINARY_DIR}/drivers/batch_chem)
  set(test_name batch_chem_${benchmark}__smoke)
  add_test(${test_name} ${exe} ${input_file})
  set_tests_properties(${test_name} PROPERTIES
    LABELS "batch;smoke"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endfunction()

function(add_alquimia_batch_chem_reference benchmark input_file)
  set(exe ${PROJECT_BINARY_DIR}/drivers/batch_chem)
  set(tools_dir ${PROJECT_SOURCE_DIR}/benchmarks/tools)
  set(reference_dir ${CMAKE_CURRENT_SOURCE_DIR}/references)
  set(results_dir ${CMAKE_CURRENT_BINARY_DIR}/results/${benchmark})
  set(test_name batch_chem_${benchmark}__engine_reference)
  add_test(
    ${test_name}
    ${PYTHON_EXECUTABLE} ${tools_dir}/run_batch_reference_test.py
      --batch-chem ${exe}
      --benchmark-dir ${CMAKE_CURRENT_SOURCE_DIR}
      --manifest ${reference_dir}/manifest.json
      --case ${benchmark}
      --work-dir ${results_dir})
  set_tests_properties(${test_name} PROPERTIES
    LABELS "batch;engine-reference"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endfunction()
