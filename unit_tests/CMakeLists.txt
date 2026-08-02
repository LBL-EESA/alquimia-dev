include(add_alquimia_unit_test)
add_alquimia_unit_test(test_alquimia_c_utils_exe test_alquimia_c_utils.c)
add_alquimia_unit_test(test_input_util_exe
  test_input_util.c
  ../drivers/input_util.c
  ../drivers/ini.c)

add_test(
  test_batch_results_py
  ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/python/test_batch_results.py)
set_tests_properties(test_batch_results_py PROPERTIES
  ENVIRONMENT "PYTHONPATH=${PROJECT_SOURCE_DIR}/benchmarks/tools"
  LABELS "unit"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
