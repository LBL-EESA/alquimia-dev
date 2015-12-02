# This function adds a unit test executable for Alquimia.
function(add_alquimia_unit_test exe)
  add_executable(${exe} ${ARGN})
  target_link_libraries(${exe} alquimia;${ALQUIMIA_TPLS})
  set_target_properties(${exe} PROPERTIES COMPILE_FLAGS "-DCMAKE_CURRENT_SOURCE_DIR=\\\"${CMAKE_CURRENT_SOURCE_DIR}\\\"")
  add_test(${exe} ${exe})
  set_tests_properties(${exe} PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endfunction()
