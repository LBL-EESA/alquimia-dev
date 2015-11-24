include_directories(${ALQUIMIA_INCLUDE_DIRS})
link_libraries(alquimia;${ALQUIMIA_TPLS})
add_executable(batch_chem batch_chem.cc cfg_reader.cc demo_containers.cc 
                          demo_output.cc demo_utils.cc string_tokenizer.cc)

# Install executable.
install(TARGETS batch_chem DESTINATION bin)

