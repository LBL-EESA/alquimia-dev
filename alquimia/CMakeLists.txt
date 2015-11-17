# Build a configuration header file from our options.
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/alquimia.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/alquimia.h"
)

add_library(alquimia alquimia.c alquimia_constants.c alquimia_interface.c 
            alquimia_memory.c alquimia_util.c alquimia_containers.F90 
            c_f_interface_module.F90)

# Install library.
install(TARGETS alquimia DESTINATION lib)

# Install headers, including generated files.
file(GLOB HEADERS *.h)
LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/alquimia.h)
LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/alquimia_version.h)
install(FILES ${HEADERS} DESTINATION include/alquimia)
