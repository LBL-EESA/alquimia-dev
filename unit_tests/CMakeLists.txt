include(add_alquimia_unit_test)
add_alquimia_unit_test(test_alquimia_c_utils_exe test_alquimia_c_utils.c)
add_alquimia_unit_test(test_input_util_exe
  test_input_util.c
  ../drivers/input_util.c
  ../drivers/ini.c)
